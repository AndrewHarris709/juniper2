//
// Created by Andrew on 2025-09-16.
//

#ifndef SIMDATA_H
#define SIMDATA_H

#include <string>

#include "SimConfig.h"


class SimData {

public:
    explicit SimData(const std::string& filename, const std::string& configname);
    ~SimData();

    inline float& x(int i) { return xyzh[4*i+0]; }
    inline float& y(int i) { return xyzh[4*i+1]; }
    inline float& z(int i) { return xyzh[4*i+2]; }
    inline float& h(int i) { return xyzh[4*i+3]; }

    // const overloads, for use on const SimData& / read-only access
    inline float x(int i) const { return xyzh[4*i+0]; }
    inline float y(int i) const { return xyzh[4*i+1]; }
    inline float z(int i) const { return xyzh[4*i+2]; }
    inline float h(int i) const { return xyzh[4*i+3]; }

    int particleCount;
    double time;
    float* xyzh = nullptr;
    float* vxyzu = nullptr;
    float* fxyz = nullptr;

    float* density = nullptr;
    float* omega = nullptr;
    float* accs = nullptr;
    float* energies = nullptr;

    SimConfig config;
    SimConfigDevice* deviceConfig;

    int getParticleCount() const;
    float getMass();

    void toCSV(const std::string& filename);
    void densityIterate();

    SimConfigDevice* getOffloadConfig();
    void generateReports();

    void updateGPU();
    void pullGPU();
};


#endif //SIMDATA_H