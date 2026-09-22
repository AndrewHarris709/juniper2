//
// Created by Andrew on 2026-07-17.
//

#include "ReportBuilders.h"

#include <map>
#include <sstream>
#include <fstream>

#include "LeapfrogIntegrator.h"
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
            GPU::getNeighbours(tree->data, nodeIndex, &data.config.g_config, data.h(partIndex), &neighbours);

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
            SimConfigDevice* config = data.getOffloadConfig();
            if (GPU::distBetween(config, data.x(i), data.x(j), data.y(i), data.y(j), data.z(i), data.z(j)) / data.h(i) < 2) {
                ss << "," << j;
            }
        }
        ss << std::endl;
    }

    std::string outputName = ".//output//" + name + "_neighbours.jrep";

    return output(ss, outputName);
}

bool ReportBuilders::buildAccelerationEnergyReport(SimData& data, std::string& name) {
    Tree* tree = new Tree(data.getParticleCount());
    tree->build(data);
    TreeData* td = tree->data;
    float* xyzh = data.xyzh;
    float* vxyzu = data.vxyzu;
    SimConfigDevice* config = data.getOffloadConfig();

    auto* accs = new junipermath::Point3f[data.getParticleCount()]{{0}};
    auto* energy = new float[data.getParticleCount()]{0};
    auto* density = data.density;
    auto* omega = data.omega;

    // transfer omega & density so they're accessible after this block too
    for (int nodeIndex = 0; nodeIndex < tree->getNodeCount(); nodeIndex++) {
        if (td->contents[nodeIndex].leftChild != -1) {
            // non-leaf node
            continue;
        }

        NodeRange partIndices = GPU::getPartsFromNode(td, nodeIndex);
        for (int i = 0; i < partIndices.size; i++) {
            int partIndex = partIndices.data[i];
            density[partIndex] = LeapfrogIntegrator::densityForParticle(partIndex, nodeIndex, xyzh, td, config);
            omega[partIndex] = LeapfrogIntegrator::omegaForParticle(partIndex, nodeIndex, density[partIndex], xyzh, td, config);
        }
    }

    #pragma omp target teams distribute parallel for  map(to: vxyzu[0:data.getParticleCount() * 4]) \
                                                      map(tofrom: accs[0:data.getParticleCount()], energy[0:data.getParticleCount()]) \
                                                      map(to: density[0:data.getParticleCount()], omega[0:data.getParticleCount()]) \
                                                      map(to: td->nodeCount, td->partCount) \
                                                      map(to: td->contents[0:td->nodeCount]) \
                                                      map(to: td->mapping[0:td->partCount])
    for (int nodeIndex = 0; nodeIndex < td->nodeCount; nodeIndex++) {
        if (td->contents[nodeIndex].leftChild != -1) {
            // non-leaf node
            continue;
        }

        NodeRange partIndices = GPU::getPartsFromNode(td, nodeIndex);
        for (int i = 0; i < partIndices.size; i++) {
            int partIndex = partIndices.data[i];
            accs[partIndex] = LeapfrogIntegrator::accForParticle(partIndex, nodeIndex, xyzh, vxyzu, density, omega, td, config);
            energy[partIndex] = LeapfrogIntegrator::energyChangeForParticle(partIndex, nodeIndex, xyzh, vxyzu, density, omega, td, config);
        }
    }

    std::stringstream ss;
    for (int i = 0; i < data.getParticleCount(); i++) {
        ss << i << "," << energy[i] << "," << accs[i].x << "," << accs[i].y << "," << accs[i].z << std::endl;
    }

    std::string outputName = ".//output//" + name + "_accenergy.jrep";

    return output(ss, outputName);
}
