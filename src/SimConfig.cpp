//
// Created by Andrew on 2026-07-15.
//

#include "toml.hpp"
#include "SimConfig.h"

#include <iostream>

SimConfig::SimConfig(const std::string &file) {
    boundingBox = junipermath::Box3f();

    toml::table tbl;
    try {
        tbl = toml::parse_file(file);
        this->mass = tbl["simconfig"]["mass"].value<float>().value();
        this->name = tbl["simconfig"]["name"].value<std::string>().value();
        this->boundingBox.x1 = tbl["simconfig"]["limits"][0][0].value<float>().value();
        this->boundingBox.x2 = tbl["simconfig"]["limits"][0][1].value<float>().value();
        this->boundingBox.y1 = tbl["simconfig"]["limits"][1][0].value<float>().value();
        this->boundingBox.y2 = tbl["simconfig"]["limits"][1][1].value<float>().value();
        this->boundingBox.z1 = tbl["simconfig"]["limits"][2][0].value<float>().value();
        this->boundingBox.z2 = tbl["simconfig"]["limits"][2][1].value<float>().value();
    }
    catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

junipermath::Box3f SimConfig::getBoundingBox() {
    return boundingBox;
}

float SimConfig::getMass() {
    return mass;
}

std::string& SimConfig::getName() {
    return name;
}

SimConfigDevice SimConfig::getOffloadConfig() {
    return SimConfigDevice{mass, boundingBox};
}
