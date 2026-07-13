#include <cmath>
#include <iostream>

#include "kernel.h"
#include "Tree.h"

//#pragma declare target start
namespace GPU {
    constexpr float NORM = 1 / M_PI;

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
        return -1.0 * NORM * (3 * valueAt(q) + q * gradientAt(q));
    }

    NodeRange getPartsFromNode(TreeNode* treeContents, int* mapping, int nodeIndex) {
        TreeNode *node = &treeContents[nodeIndex];

        NodeRange range;
        range.data = mapping + node->mappingStart;
        range.size = node->mappingSize;

        return range;
    }

    float distBetween(float x1, float x2, float y1, float y2, float z1, float z2) {
        // Periodic Boundary-aware distance calculation
        float dx = std::abs(x1 - x2);
        float dy = std::abs(y1 - y2);
        float dz = std::abs(z1 - z2);

        float temp_bounds[6]{0, 1, 0, 1, 0, 1};

        if (dx > (temp_bounds[1] - temp_bounds[0]) / 2) {
            dx = (temp_bounds[1] - temp_bounds[0]) - dx;
        }
        if (dy > (temp_bounds[3] - temp_bounds[2]) / 2) {
            dy = (temp_bounds[3] - temp_bounds[2]) - dy;
        }
        if (dz > (temp_bounds[5] - temp_bounds[4]) / 2) {
            dz = (temp_bounds[5] - temp_bounds[4]) - dz;
        }

        if (x1 > x2)
            dx *= -1;
        if (y1 > y2)
            dy *= -1;
        if (z1 > z2)
            dz *= -1;

        return sqrtf(dx * dx + dy * dy + dz * dz);
    }

    float distBetween(float* xyzh, int part1, int part2) {
        float x1 = xyzh[4 * part1], y1 = xyzh[4 * part1 + 1], z1 = xyzh[4 * part1 + 2];
        float x2 = xyzh[4 * part2], y2 = xyzh[4 * part2 + 1], z2 = xyzh[4 * part2 + 2];

        return distBetween(x1, x2, y1, y2, z1, z2);
    }

    float distBetweenNodes(TreeNode* treeContents, int nodeIndex1, int nodeIndex2) {
        TreeNode* node1 = &treeContents[nodeIndex1];
        TreeNode* node2 = &treeContents[nodeIndex2];

        return distBetween(node1->x, node2->x, node1->y, node2->y, node1->z, node2->z);
    }

    bool atEndCondition(float newH, float oldH, float origH) {
        return std::abs(newH - oldH) / origH <= 10e-4;
    }

    bool* getNeighbours(TreeNode* treeContents, int* treeMapping, int nodeCount, int partCount, int nodeIndex, float partMax) {
        bool* isNeighbour = new bool[partCount]();

        int* nodeStack = new int[nodeCount];
        int stackSize = 0;

        nodeStack[stackSize++] = 0;
        TreeNode* targetNode = &treeContents[nodeIndex];

        auto *kernel = new Kernel();

        while (stackSize > 0) {
            int nextNodeIndex = nodeStack[--stackSize];
            TreeNode* nextNode = &treeContents[nextNodeIndex];

            float distance = distBetweenNodes(treeContents, nextNodeIndex, nodeIndex);
            float targetBounds = nextNode->size + targetNode->size + (kernel->getRadius() * std::max(targetNode->hmax, partMax));

            if (distance * distance < targetBounds * targetBounds) {
                if (nextNode->leftChild == -1) {
                    NodeRange nodeContents = getPartsFromNode(treeContents, treeMapping, nextNodeIndex);
                    for (int i = 0; i < nodeContents.size; i++) {
                        int partIndex = nodeContents.data[i];
                        isNeighbour[partIndex] = true;
                    }
                } else {
                    nodeStack[stackSize++] = nextNode->leftChild;
                    nodeStack[stackSize++] = nextNode->rightChild;
                }
            }
        }

        return isNeighbour;
    }

    float densityIterateAtParticle(float* xyzh, TreeNode* treeContents, int* treeMapping, int nodeCount, int partCount, int partIndex, int nodeIndex, float m) {
        int partA = partIndex;
        float oldH = std::numeric_limits<float>::max();
        float newH = xyzh[partA * 4 + 3];
        float origH = newH;
        int iterationCount = 0;

        bool* isPartANeighbour = new bool[partCount];

        while (!atEndCondition(newH, oldH, origH)) {
            isPartANeighbour = getNeighbours(treeContents, treeMapping, nodeCount, partCount, nodeIndex, newH);

            float density = m * std::pow(1.2 / newH, 3);
            float grad = -newH / (3 * density);

            float density_sum = 0;
            float omega = 0;
            for (int i = 0; i < partCount; i++) {
                int partB = treeMapping[i];

                if (!isPartANeighbour[partB]) {
                    continue;
                }

                float dist = distBetween(xyzh, partA, partB);

                density_sum += m * GPU::valueAt(dist / newH) / std::pow(newH, 3);
                omega += m * GPU::dWdhAt(dist / newH);
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

        delete isPartANeighbour;

        return newH;
    }


}

//#pragma declare target end
