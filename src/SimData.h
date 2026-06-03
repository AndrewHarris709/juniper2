//
// Created by Andrew on 2025-09-16.
//

#ifndef SIMDATA_H
#define SIMDATA_H
#include <string>


class SimData {
public:
    explicit SimData(const std::string& filename);
    explicit SimData() : SimData("") {};
    ~SimData();

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
};


#endif //SIMDATA_H