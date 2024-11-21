#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>

#include "parameters.h"

namespace ccp {

class Simulation {
public:
    explicit Simulation(const Parameters& parameters, const std::string& data_path)
        : parameters_(parameters), data_path_(data_path) {};
    void run();

private:
    Parameters parameters_;
    std::string data_path_;
};

} // ccp

#endif //SIMULATION_H
