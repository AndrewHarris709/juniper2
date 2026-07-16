//
// Created by Andrew on 2026-07-13.
//

#ifndef GPUOPERATORS_H
#define GPUOPERATORS_H
#include <math.h>

#include "Tree.h"

#pragma omp declare target
namespace GPU {
    constexpr float NORM = 1 / M_PI;

    float valueAt(float q);
    float gradientAt(float q);
    float dWdhAt(float q);
    NodeRange getPartsFromNode(TreeNode* treeContents, int* mapping, int nodeIndex);
    float distBetween(float x1, float x2, float y1, float y2, float z1, float z2);
    float distBetweenNodes(TreeNode* treeContents, int nodeIndex1, int nodeIndex2);
    bool atEndCondition(float newH, float oldH, float origH);
    void getNeighbours(TreeNode* treeContents, int* treeMapping, int nodeCount, int partCount, int nodeIndex, int partIndex, SimConfigDevice config, float partMax, bool* neighbourScratch, int* stackScratch);
    float densityIterateAtParticle(float* xyzh, TreeNode* treeContents, int* treeMapping, int nodeCount, int partCount, int partIndex, int nodeIndex, SimConfigDevice config, bool* neighbourScratch, int* stackScratch);

}
#pragma omp end declare target

#endif //GPUOPERATORS_H
