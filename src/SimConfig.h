//
// Created by Andrew on 2026-07-15.
//

#ifndef SIMCONFIG_H
#define SIMCONFIG_H
#include <string>

#include "commons.h"

struct SimConfigDevice {
    float mass;
    junipermath::Box3f boundingBox;
};

class SimConfig {
    junipermath::Box3f boundingBox;
    float mass;
    std::string name;

public:
    explicit SimConfig(const std::string& file);

    junipermath::Box3f getBoundingBox();
    float getMass();
    std::string& getName();
    SimConfigDevice getOffloadConfig();
};

#endif //SIMCONFIG_H
