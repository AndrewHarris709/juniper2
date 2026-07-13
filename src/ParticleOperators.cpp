//
// Created by Andrew on 2025-09-12.
//

#include <cmath>

#include "commons.h"
#include "SimData.h"
#include "Tree.h"

namespace PartOps {
    junipermath::Point3f getCentreOfMass(SimData& data, Tree& tree, int index) {
        junipermath::Point3f centreOfMass(0, 0, 0);
        NodeRange contents = tree.getPartsFromNode(index);

        int partCount;
        for (partCount = 0; partCount < contents.size; partCount++) {
            int partIndex = contents.data[partCount];
            centreOfMass.x += data.x(partIndex);
            centreOfMass.y += data.y(partIndex);
            centreOfMass.z += data.z(partIndex);
        }

        centreOfMass.x /= partCount;
        centreOfMass.y /= partCount;
        centreOfMass.z /= partCount;

        return centreOfMass;
    }

    junipermath::Box3f getBoundingBox(SimData& data, Tree& tree, int index) {
        NodeRange contents = tree.getPartsFromNode(index);
        int fi = contents.data[0];

        junipermath::Box3f box(data.x(fi), data.y(fi), data.z(fi),
                              data.x(fi), data.y(fi), data.z(fi));
        for (int i = 0; i < contents.size; i++) {
            int p = contents.data[i];

            float x = data.x(p);
            float y = data.y(p);
            float z = data.z(p);

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
        NodeRange contents = tree.getPartsFromNode(index);
        junipermath::Point3f centre = getCentreOfMass(data, tree, index);
        float boundingRadius = 0;

        for (int i = 0; i < contents.size; i++) {
            int p = contents.data[i];
            float x = data.x(p);
            float y = data.y(p);
            float z = data.z(p);
            float centreDist = std::sqrt(((centre.x - x) * (centre.x - x) +
                                (centre.y - y) * (centre.y - y) +
                                (centre.z - z) * (centre.z - z)));

            if (centreDist > boundingRadius) boundingRadius = centreDist;
        }

        return boundingRadius;
    }

    float getMaxSmoothingLength(SimData& data, Tree& tree, int index) {
        NodeRange contents = tree.getPartsFromNode(index);
        float maxSmoothingLength = 0;

        for (int i = 0; i < contents.size; i++) {
            int p = contents.data[i];
            if (data.h(p) > maxSmoothingLength) {
                maxSmoothingLength = data.h(i);
            }
        }

        return maxSmoothingLength;
    }
}
