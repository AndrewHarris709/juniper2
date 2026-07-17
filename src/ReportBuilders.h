//
// Created by Andrew on 2026-07-17.
//

#ifndef REPORTBUILDERS_H
#define REPORTBUILDERS_H
#include <string>
#include <SimData.h>

class ReportBuilders {
    static bool output(std::stringstream& contents, std::string& filename);

public:
    static bool buildNeighboursReport(SimData& data, std::string& name);
};

#endif //REPORTBUILDERS_H
