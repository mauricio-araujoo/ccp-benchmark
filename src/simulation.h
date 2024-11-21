#ifndef SIMULATION_H
#define SIMULATION_H

#include "parameters.h"

namespace ccp {

class Simulation {
public:
    explicit Simulation(const Parameters& parameters) : parameters_(parameters) {};
    void run();

private:
    Parameters parameters_;
};

} // ccp

#endif //SIMULATION_H
