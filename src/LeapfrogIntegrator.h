//
// Created by Andrew on 2026-07-31.
//

#ifndef LEAPFROGINTEGRATOR_H
#define LEAPFROGINTEGRATOR_H

#include "SimConfig.h"
#include "Tree.h"

#pragma omp begin declare target
namespace LeapfrogIntegrator {

    constexpr float GAMMA = 5.0 / 3;

    float densityForParticle(int partIndex, float *xyzh, TreeData *tree, SimConfigDevice *simConfig);
    float omegaForParticle(int partIndex, int nodeIndex, float *xyzh, TreeData *tree, SimConfigDevice *simConfig);
    float pressureforParticle(int partIndex, float *xyzh, TreeData *tree, SimConfigDevice *simConfig);
    junipermath::Point3f accForParticle(int partIndex, int nodeIndex, float *xyzh, float *vxyzu, TreeData* tree, SimConfigDevice* config);
    float energyChangeForParticle(int partIndex, int nodeIndex, float *xyzh, float *vxyzu, TreeData* tree, SimConfigDevice* config);
    float qabForParticle(int partIndex, int otherPart, float *xyzh, TreeData* tree, SimConfigDevice* config);
    junipermath::Point3f displacementBetween(int partIndex, int otherPart, float* xyzh, SimConfigDevice* config);
    junipermath::Point3f displacementBetween(float x1, float x2, float y1, float y2, float z1, float z2, SimConfigDevice* config);
    junipermath::Point3f velocityDiffBetween(int partIndex, int otherPart, float *vxyzu);
    junipermath::Point3f norm(junipermath::Point3f p);
    float dot(const junipermath::Point3f p1, const junipermath::Point3f p2);

};
#pragma omp end declare target

#endif //LEAPFROGINTEGRATOR_H
