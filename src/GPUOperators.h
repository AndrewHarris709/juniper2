//
// Created by Andrew on 2026-07-13.
//

#ifndef GPUOPERATORS_H
#define GPUOPERATORS_H

#include "SimConfig.h"
#include "Tree.h"

#pragma omp begin declare target
namespace GPU {
    constexpr float NORM = 1 / M_PI;
    constexpr int MAX_NEIGHBOURS = 30000;

    struct NeighbourList {
        int indices[MAX_NEIGHBOURS];
        int count;
    };

    float valueAt(float q);
    float gradientAt(float q);
    float dWdhAt(float q);
    NodeRange getPartsFromNode(TreeData* tree, int nodeIndex);
    float distBetween(SimConfigDevice* config, float x1, float x2, float y1, float y2, float z1, float z2);
    float distBetweenNodes(TreeData* tree, int nodeIndex1, int nodeIndex2);
    bool atEndCondition(float newH, float oldH, float origH);
    void getNeighbours(TreeData* tree, int nodeIndex, SimConfigDevice* config, float partMax, NeighbourList* result);
    float densityIterateAtParticle(float* xyzh, TreeData* tree, int partIndex, int nodeIndex, SimConfigDevice* config);

}
#pragma omp end declare target

#endif //GPUOPERATORS_H
