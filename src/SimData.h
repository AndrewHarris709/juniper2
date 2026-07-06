//
// Created by Andrew on 2025-09-16.
//

#ifndef SIMDATA_H
#define SIMDATA_H

#include <string>

#include "Tree.h"


class SimData {
    float densityIterateAtParticle(SimData& simData, Tree& tree, int partIndex, int nodeIndex);
    bool atEndCondition(float newH, float oldH, float origH);
    float distBetween(int part1, int part2);

public:
    explicit SimData(const std::string& filename);
    explicit SimData() : SimData("") {};
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
    float m;
    float* xyzh = nullptr;
    float* vxyzu = nullptr;
    float* fxyz = nullptr;

    float* density = nullptr;
    float* omega = nullptr;
    float* accs = nullptr;
    float* energies = nullptr;

    int getParticleCount() const;

    void toCSV(const std::string& filename);
    void densityIterate(SimData& simData);
};


#endif //SIMDATA_H