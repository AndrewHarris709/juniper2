//
// Created by Andrew on 2026-07-17.
//

#include "ReportBuilders.h"

#include <map>
#include <sstream>
#include <fstream>

#include "GPUOperators.h"
#include "Tree.h"

bool ReportBuilders::output(std::stringstream& ss, std::string &filename) {
    std::ofstream outFile(filename);

    if (outFile.is_open()) {
        outFile << ss.rdbuf();
        outFile.close();
        return true;
    }

    return false;
}

bool ReportBuilders::buildNeighboursReport(SimData& data, std::string& name) {
    Tree* tree = new Tree(data.getParticleCount());
    tree->build(data);

    std::map<int, std::vector<int>> treeNeighbours;

    for (int nodeIndex = 0; nodeIndex < tree->getNodeCount(); nodeIndex++) {
        if (!tree->isLeaf(nodeIndex)) {
            continue;
        }

        NodeRange partIndices = GPU::getPartsFromNode(tree->data, nodeIndex);
        for (int i = 0; i < partIndices.size; i++) {
            int partIndex = partIndices.data[i];

            GPU::NeighbourList neighbours;
            GPU::getNeighbours(tree->data, nodeIndex, data.getOffloadConfig(), data.h(partIndex), &neighbours);

            std::vector<int> neighbourVector(neighbours.indices, neighbours.indices + neighbours.count);
            treeNeighbours.emplace(partIndex, neighbourVector);
        }
    }

    std::stringstream ss;
    for (int i = 0; i < data.getParticleCount(); i++) {
        ss << i << ",tree";
        for (const int neighbour: treeNeighbours[i]) {
            ss << "," << neighbour;
        }
        ss << std::endl;

        ss << i << ",true";
        for (int j = 0; j < data.getParticleCount(); j++) {
            SimConfigDevice config = data.getOffloadConfig();
            if (GPU::distBetween(&config, data.x(i), data.x(j), data.y(i), data.y(j), data.z(i), data.z(j)) / data.h(i) < 2) {
                ss << "," << j;
            }
        }
        ss << std::endl;
    }

    std::string outputName = ".//output//" + name + "_neighbours.jrep";

    return output(ss, outputName);
}
