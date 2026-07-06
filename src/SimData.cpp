//
// Created by Andrew on 2025-09-16.
//

#include "SimData.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ranges>
#include <vector>
#include <stack>

#include "Kernel.h"
#include "Tree.h"
#include "commons.h"

constexpr float HFACT = 1.2;
constexpr int MAX_DENSITY_ITERATIONS = 10;

static inline const std::vector<std::string> posCols{"x", "y", "z", "h"};
static inline const std::vector<std::string> velCols{"vx", "vy", "vz", "u"};
static inline const std::vector<std::string> varCols{"fx", "fy", "fz"};

//#pragma omp declare target
float SimData::distBetween(int part1, int part2) {
    float x1 = this->xyzh[4 * part1], y1 = this->xyzh[4 * part1 + 1], z1 = this->xyzh[4 * part1 + 2];
    float x2 = this->xyzh[4 * part2], y2 = this->xyzh[4 * part2 + 1], z2 = this->xyzh[4 * part2 + 2];

    return junipermath::distBetween(x1, x2, y1, y2, z1, z2);
}

bool* getNeighbours(SimData& simData, Tree& tree, int nodeIndex) {
    bool* isNeighbour = new bool[simData.getParticleCount()]();

    std::vector<int> neighbours;
    std::stack<int> nodeStack;
    nodeStack.push(0);
    TreeNode* targetNode = tree.getNode(nodeIndex);

    auto *kernel = new Kernel();

    while (!nodeStack.empty()) {
        int nextNodeIndex = nodeStack.top();
        TreeNode* nextNode = tree.getNode(nextNodeIndex);
        nodeStack.pop();

        float distance = tree.distBetweenNodes(nextNodeIndex, nodeIndex);
        float targetBounds = nextNode->size + targetNode->size + (kernel->getRadius() * std::max(targetNode->hmax, nextNode->hmax));

        if (distance * distance < targetBounds * targetBounds) {
            if (nextNode->leftChild == -1) {
                NodeRange nodeContents = tree.getPartsFromNode(nextNodeIndex);
                for (int i = 0; i < nodeContents.size; i++) {
                    int partIndex = nodeContents.data[i];
                    isNeighbour[partIndex] = true;
                }
            } else {
                nodeStack.push(nextNode->leftChild);
                nodeStack.push(nextNode->rightChild);
            }
        }
    }

    return isNeighbour;
}

bool SimData::atEndCondition(float newH, float oldH, float origH) {
    return std::abs(newH - oldH) / origH <= 10e-4;
}

float SimData::densityIterateAtParticle(SimData& simData, Tree& tree, int partIndex, int nodeIndex) {
    int partA = partIndex;
    float oldH = std::numeric_limits<float>::max();
    float newH = simData.xyzh[partA * 4 + 3];
    float origH = newH;
    int iterationCount = 0;

    bool* isPartANeighbour = new bool[simData.getParticleCount()];

    while (!atEndCondition(newH, oldH, origH)) {
        isPartANeighbour = getNeighbours(simData, tree, nodeIndex);
        Kernel *kernel = new Kernel();

        float density = simData.m * std::pow(HFACT / newH, 3);
        float grad = -newH / (3 * density);

        float density_sum = 0;
        float omega = 0;
        for (int i = 0; i < simData.getParticleCount(); i++) {
            int partB = tree.data->mapping[i];

            if (!isPartANeighbour[partB]) {
                continue;
            }

            float dist = simData.distBetween(partA, partB);
            density_sum += simData.m * kernel->valueAt(dist / newH) / std::pow(newH, 3);
            omega += simData.m * kernel->dWdhAt(dist / newH);
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

        if (iterationCount > MAX_DENSITY_ITERATIONS) {
            break;
        }
    }

    delete isPartANeighbour;

    return newH;
}
//#pragma omp end declare target

void SimData::densityIterate(SimData& simData) {
    Tree* tree = new Tree(simData.getParticleCount());
    tree->build(simData);
    std::cout << "Tree built with " << tree->getNodeCount() << "nodes." << std::endl;

    float m = simData.m;
    float* xyzh = simData.xyzh;
    TreeNode* treeContents = tree->data->contents;
    int* mapping = tree->data->mapping;
    int nodeCount = tree->data->nodeCount;
    int partCount = tree->data->partCount;

    //#pragma omp target enter data map(to: xyzh[0:pCount*4])
    //#pragma omp target teams distribute parallel for map(tofrom: xyzh[0:pCount*4]) \
    //                                                 map(to: treeContents[0:nodeCount], mapping[0:partCount]) \
    //                                                 map(to: nodeCount, partCount) \
    //                                                 defaultmap(none)
    for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
        if (treeContents[nodeIndex].leftChild != -1) {
            // non-leaf node
            continue;
        }

        NodeRange partIndices = tree->getPartsFromNode(nodeIndex);
        for (int i = 0; i < partIndices.size; i++) {
            int partIndex = partIndices.data[i];
            xyzh[4 * partIndex + 3] = densityIterateAtParticle(simData, *tree, partIndex, nodeIndex);
        }
    }

    delete tree;
}

SimData::SimData(const std::string& filename) {
    time = 0.0;

    if (filename.empty()) {
        return;
    }

    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << std::filesystem::current_path().string() << std::endl;
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    std::vector<std::vector<std::string>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            // Linux systems will sometimes return a '/r' character from getline.
            std::erase(cell, '\r' );
            row.push_back(cell);
        }

        data.push_back(row);
    }

    particleCount = data.size() - 1;
    xyzh    = new float[particleCount * 4];
    vxyzu   = new float[particleCount * 4];
    fxyz    = new float[particleCount * 3];
    int rowIndex = 0;
    for (const auto& row : data | std::views::drop(1)) {

        int i = 0;
        for (const auto& cell : row) {
            std::string column = data[0][i];

            if (auto posOffset = std::ranges::find(posCols, column); posOffset != std::end(posCols)) {
                xyzh[rowIndex * 4 + std::distance(posCols.begin(), posOffset)] = std::stof(cell);
            }
            if (auto velOffset = std::ranges::find(velCols, column); velOffset != std::end(velCols)) {
                vxyzu[rowIndex * 4 + std::distance(velCols.begin(), velOffset)] = std::stof(cell);
            }
            if (auto varOffset = std::ranges::find(varCols, column); varOffset != std::end(varCols)) {
                fxyz[rowIndex * 3 + std::distance(varCols.begin(), varOffset)] = std::stof(cell);
            }
            i++;
        }
        rowIndex++;
    }

    density = new float[particleCount]();
    omega = new float[particleCount]();
    accs = new float[particleCount]();
    energies = new float[particleCount]();
    this->particleCount = particleCount;

    file.close();
}

SimData::~SimData() {
    delete[] xyzh;
    delete[] vxyzu;
    delete[] fxyz;
    delete[] density;
    delete[] omega;
    delete[] accs;
    delete[] energies;
}

void SimData::toCSV(const std::string &filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cout << std::filesystem::current_path().string() << std::endl;
        std::cerr << "Error writing to file " << filename << std::endl;
    }

    file << "i,x,y,z,h,vx,vy,vz,u,fx,fy,fz" << std::endl;

    for (int i = 0; i < getParticleCount(); i++) {
        file << i << ",";
        file << xyzh[4*i] << "," << xyzh[4*i+1] << "," << xyzh[4*i+2] << "," << xyzh[4*i+3] << ",";
        file << vxyzu[4*i] << "," << vxyzu[4*i+1] << "," << vxyzu[4*i+2] << "," << vxyzu[4*i+3] << ",";
        if (getParticleCount() > 0) {
            file << fxyz[3*i] << "," << fxyz[3*i+1] << "," << fxyz[3*i+2];
        }
        file << std::endl;
    }

    file.close();
}

int SimData::getParticleCount() const {
    return this->particleCount;
}

