//
// Created by Andrew on 2026-07-15.
//

#ifndef SIMCONFIG_H
#define SIMCONFIG_H
#include <string>
#include <vector>

#include "commons.h"

struct SimConfigDevice {
    float mass;
    junipermath::Box3f boundingBox;
};

class SimConfig {
    junipermath::Box3f boundingBox;
    float mass;
    std::string name;
    std::vector<std::string> reports;

public:
    explicit SimConfig(const std::string& file);

    junipermath::Box3f getBoundingBox();
    float getMass();
    std::string& getName();
    SimConfigDevice getOffloadConfig();
    bool shouldReport(const std::string& report);
};

#endif //SIMCONFIG_H
