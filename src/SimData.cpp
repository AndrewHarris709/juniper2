//
// Created by Andrew on 2025-09-16.
//

#include "SimData.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <ranges>
#include <vector>

#include "Tree.h"
#include "GPUOperators.h"
#include "ReportBuilders.h"

constexpr float HFACT = 1.2;
constexpr int MAX_DENSITY_ITERATIONS = 10;
constexpr int THREAD_LIMIT = 512;

static inline const std::vector<std::string> posCols{"x", "y", "z", "h"};
static inline const std::vector<std::string> velCols{"vx", "vy", "vz", "u"};
static inline const std::vector<std::string> varCols{"fx", "fy", "fz"};

void SimData::densityIterate() {
    Tree* tree = new Tree(this->getParticleCount());
    tree->build(*this);
    std::cout << "Tree built with " << tree->getNodeCount() << "nodes." << std::endl;

    float* xyzh = this->xyzh;
    TreeData* td = tree->data;
    SimConfigDevice* deviceConfig = this->deviceConfig;

    assert(omp_target_is_present(xyzh, omp_get_default_device()) && "xyzh not mapped on device");
    assert(omp_target_is_present(deviceConfig, omp_get_default_device()) && "config not mapped on device");

    #pragma omp target teams distribute parallel for  map(to: td->nodeCount, td->partCount) \
                                                      map(to: td->contents[0:td->nodeCount]) \
                                                      map(to: td->mapping[0:td->partCount]) \
                                                      thread_limit(THREAD_LIMIT)
    for (int nodeIndex = 0; nodeIndex < td->nodeCount; nodeIndex++) {
        if (td->contents[nodeIndex].leftChild != -1) {
            // non-leaf node
            continue;
        }

        NodeRange partIndices = GPU::getPartsFromNode(td, nodeIndex);
        for (int i = 0; i < partIndices.size; i++) {
            int partIndex = partIndices.data[i];
            float result = GPU::densityIterateAtParticle(xyzh, td, partIndex, nodeIndex, deviceConfig);
            xyzh[4 * partIndex + 3] = result;
        }
    }

    pullGPU();
}

SimData::SimData(const std::string& filename, const std::string& configname): config(configname) {
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

    this->deviceConfig = getOffloadConfig();
    #pragma omp target enter data map(to: xyzh[0:this->getParticleCount()*4], *this->deviceConfig)
}

SimData::~SimData() {
    #pragma omp target exit data map(release: xyzh[0:this->getParticleCount()*4], *this->deviceConfig)

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

float SimData::getMass() {
    return this->config.getMass();
}

SimConfigDevice* SimData::getOffloadConfig() {
    return &config.g_config;
}

void SimData::generateReports() {
    if (this->config.shouldReport("neighbours")) {
        std::cout << "Generating Particle Neighbours Report...";
        ReportBuilders::buildNeighboursReport(*this, this->config.getName());
        std::cout << "Done" << std::endl;
    }

    if (this->config.shouldReport("accenergy")) {
        std::cout << "Generating Acceleration and Energy Report...";
        ReportBuilders::buildAccelerationEnergyReport(*this, this->config.getName());
        std::cout << "Done" << std::endl;
    }
}

void SimData::updateGPU() {
    #pragma omp target update to(xyzh[0:this->getParticleCount()*4])
}

void SimData::pullGPU() {
    #pragma omp target update from(xyzh[0:this->getParticleCount()*4])
}

