// This file defines and calculate the energy terminology based on the Gromacs computation tool.

#ifndef ENERGY_H
#define ENERGY_H

#include <string>
// Use system command here to stop wasting time.
// TODO: add args to the main functions in Python main.

class energy {
    float energy;
    std::string filePath;

    float computeEnergy();

    static std::string getDirectory(std::string fileName);
};


#endif //MCSOFTWARE_ENERGY_H
