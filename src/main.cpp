#include <cstdio>

#include "simulation.h"
#include "kn/random/random.h"


int main() {
    kn::random::initialize(42);

    printf("starting benchmark simulation\n");

    ccp::Simulation sim(ccp::Parameters::case_1());
    sim.run();

    return 0;
}

