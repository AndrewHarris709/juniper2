//
// Created by Andrew on 2025-09-12.
//

#ifndef JUNIPER2_TREEOPERATORS_H
#define JUNIPER2_TREEOPERATORS_H
#include <numeric>
#include <cmath>

#include "commons.h"
#include "SimData.h"
#include "Tree.h"

namespace PartOps {
    Point3f getCentreOfMass(SimData& data, Tree& tree, int index) {
        Point3f centreOfMass(0, 0, 0);
        std::span<int> contents = tree.getNodeIndices(index);

        int partCount = 0;
        for (const int i : contents) {
            centreOfMass.x += data.xyzh[4 * i];
            centreOfMass.y += data.xyzh[4 * i + 1];
            centreOfMass.z += data.xyzh[4 * i + 2];
            partCount++;
        }

        centreOfMass.x /= partCount;
        centreOfMass.y /= partCount;
        centreOfMass.z /= partCount;

        return centreOfMass;
    }

    Box3f getBoundingBox(SimData& data, Tree& tree, int index) {
        std::span<int> contents = tree.getNodeIndices(index);
        int fi = contents[0];

        Box3f box(data.xyzh[4 * fi], data.xyzh[4 * fi + 1], data.xyzh[4 * fi + 2],
                 data.xyzh[4 * fi], data.xyzh[4 * fi + 1], data.xyzh[4 * fi + 2]);
        for (int i : contents) {
            float x = data.xyzh[4 * i];
            float y = data.xyzh[4 * i + 1];
            float z = data.xyzh[4 * i + 2];

            if (x < box.x1) box.x1 = x;
            if (y < box.y1) box.y1 = y;
            if (z < box.z1) box.z1 = z;
            if (x > box.x2) box.x2 = x;
            if (y > box.y2) box.y2 = y;
            if (z > box.z2) box.z2 = z;
        }

        return box;
    }

    float getBoundingRadius(SimData& data, Tree& tree, int index) {
        std::span<int> contents = tree.getNodeIndices(index);
        Point3f centre = getCentreOfMass(data, tree, index);
        float boundingRadius = 0;

        for (int i : contents) {
            float x = data.xyzh[4 * i];
            float y = data.xyzh[4 * i + 1];
            float z = data.xyzh[4 * i + 2];
            float centreDist = std::sqrt(((centre.x - x) * (centre.x - x) +
                                (centre.y - y) * (centre.y - y) +
                                (centre.z - z) * (centre.z - z)));

            if (centreDist > boundingRadius) boundingRadius = centreDist;
        }

        return boundingRadius;
    }

    float getMaxSmoothingLength(SimData& data, Tree& tree, int index) {
        std::span<int> contents = tree.getNodeIndices(index);
        float maxSmoothingLength = 0;

        for (int i: contents) {
            if (data.xyzh[4 * i + 3] > maxSmoothingLength) {
                maxSmoothingLength = data.xyzh[4 * i + 3];
            }
        }

        return maxSmoothingLength;
    }
}

#endif JUNIPER2_TREEOPERATORS_H
