//
// Created by Andrew on 2026-07-15.
//

#include "toml.hpp"
#include "SimConfig.h"

#include <algorithm>
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

        if (toml::array* arr = tbl["simconfig"]["reports"].as_array()) {
            for (auto&& el : *arr) {
                if (auto str = el.as_string()) {
                    this->reports.push_back(str->get());
                }
            }
        }
    }
    catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }

    g_config = SimConfigDevice{mass, boundingBox};
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

bool SimConfig::shouldReport(const std::string& report) {
    auto offset = std::ranges::find(this->reports, report);
    return offset != std::end(this->reports);
}
