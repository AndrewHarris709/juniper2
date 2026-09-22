#include "GPUOperators.h"

#include <cmath>
#include <iostream>

#include "SimConfig.h"
#include "Tree.h"

#pragma omp begin declare target
namespace GPU {
    constexpr int MAX_STACK_DEPTH = 64;

    float valueAt(const float q)
    {
        float aq = std::abs(q);

        if (aq >= 2)
        {
            return 0;
        }

        if (aq >= 1)
        {
            return NORM * 0.25 * (2 - aq) * (2 - aq) * (2 - aq);
        }

        return NORM * (1 - 1.5 * aq * aq + 0.75 * aq * aq * aq);
    }

    float gradientAt(const float q) {
        float aq = std::abs(q);

        if (aq >= 2) {
            return 0;
        }

        if (aq >= 1) {
            return NORM * -0.75 * (aq - 2) * (aq - 2);
        }

        return NORM * aq * (2.25 * aq - 3);
    }

    float dWdhAt(const float q) {
        return -1.0 * (3 * valueAt(q) + q * gradientAt(q));
    }

    NodeRange getPartsFromNode(TreeData* tree, int nodeIndex) {
        TreeNode *node = &tree->contents[nodeIndex];

        NodeRange range;
        range.data = tree->mapping + node->mappingStart;
        range.size = node->mappingSize;

        return range;
    }

    float distBetween(SimConfigDevice* config, float x1, float x2, float y1, float y2, float z1, float z2) {
        // Periodic Boundary-aware distance calculation
        float dx = std::abs(x1 - x2);
        float dy = std::abs(y1 - y2);
        float dz = std::abs(z1 - z2);

        if (dx > (config->boundingBox.x2 - config->boundingBox.x1) / 2) {
            dx = (config->boundingBox.x2 - config->boundingBox.x1) - dx;
        }
        if (dy > (config->boundingBox.y2 - config->boundingBox.y1) / 2) {
            dy = (config->boundingBox.y2 - config->boundingBox.y1) - dy;
        }
        if (dz > (config->boundingBox.z2 - config->boundingBox.z1) / 2) {
            dz = (config->boundingBox.z2 - config->boundingBox.z1) - dz;
        }

        if (x1 > x2)
            dx *= -1;
        if (y1 > y2)
            dy *= -1;
        if (z1 > z2)
            dz *= -1;

        return sqrtf(dx * dx + dy * dy + dz * dz);
    }

    float distBetween(SimConfigDevice* config, float* xyzh, int part1, int part2) {
        float x1 = xyzh[4 * part1], y1 = xyzh[4 * part1 + 1], z1 = xyzh[4 * part1 + 2];
        float x2 = xyzh[4 * part2], y2 = xyzh[4 * part2 + 1], z2 = xyzh[4 * part2 + 2];

        return distBetween(config, x1, x2, y1, y2, z1, z2);
    }

    float distBetweenNodes(SimConfigDevice* config, TreeData* tree, int nodeIndex1, int nodeIndex2) {
        TreeNode* node1 = &tree->contents[nodeIndex1];
        TreeNode* node2 = &tree->contents[nodeIndex2];

        return distBetween(config, node1->x, node2->x, node1->y, node2->y, node1->z, node2->z);
    }

    bool atEndCondition(float newH, float oldH, float origH) {
        return std::abs(newH - oldH) / origH <= 1e-4;
    }

    void getNeighbours(TreeData* tree, int nodeIndex, SimConfigDevice* config, float partMax, NeighbourList* result) {
        int stackScratch[MAX_STACK_DEPTH]{};
        result->count = 0;

        int stackSize = 0;
        stackScratch[(stackSize++)] = 0;
        TreeNode* targetNode = &(tree->contents[nodeIndex]);

        while (stackSize > 0) {
            int nextNodeIndex = stackScratch[(--stackSize)];
            TreeNode* nextNode = &tree->contents[nextNodeIndex];

            float distance = distBetweenNodes(config, tree, nextNodeIndex, nodeIndex);
            float targetBounds = nextNode->size + targetNode->size + (2 * std::max(targetNode->hmax, partMax));

            if (distance * distance < targetBounds * targetBounds) {
                if (nextNode->leftChild == -1) {
                    NodeRange nodeContents = getPartsFromNode(tree, nextNodeIndex);
                    for (int i = 0; i < nodeContents.size; i++) {
                        int realIndex = nodeContents.data[i];
                        result->indices[(result->count++)] = realIndex;
                    }
                } else {
                    stackScratch[(stackSize++)] = nextNode->leftChild;
                    stackScratch[(stackSize++)] = nextNode->rightChild;
                }
            }
        }
    }

    float densityIterateAtParticle(float* xyzh, TreeData* tree, int partIndex, int nodeIndex, SimConfigDevice* config) {
        int partA = partIndex;
        float oldH = std::numeric_limits<float>::max();
        float newH = xyzh[partA * 4 + 3];
        float origH = newH;
        int iterationCount = 0;

        while (!atEndCondition(newH, oldH, origH)) {
            NeighbourList neighbours;
            getNeighbours(tree, nodeIndex, config, newH, &neighbours);

            float density = config->mass * std::pow(1.2 / newH, 3);
            float grad = -newH / (3 * density);

            float density_sum = 0;
            float omega = 0;
            for (int i = 0; i < neighbours.count; i++) {
                int partB = neighbours.indices[i];

                float dist = distBetween(config, xyzh, partA, partB);

                density_sum += config->mass * GPU::valueAt(dist / newH) / std::pow(newH, 3);
                omega += config->mass * GPU::dWdhAt(dist / newH);
            }
            omega = 1 - grad * omega / std::pow(newH, 4);

            oldH = newH;
            newH -= (density - density_sum) / ((-3 * density * omega) / newH);
            iterationCount++;

            if (newH > 1.4 * oldH) {
                newH = 1.4 * oldH;
            } else if (newH < 0.7 * oldH) {
                newH = 0.7 * oldH;
            }

            if (iterationCount > 10) {
                break;
            }
        }

        return newH;
    }

}
#pragma omp end declare target
