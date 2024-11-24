#ifndef SIMULATION_H
#define SIMULATION_H

#include <kn/collisions/mcc.h>
#include <kn/particle/species.h>
#include <kn/spatial/grid.h>

#include <string>

#include "events.h"
#include "parameters.h"

namespace ccp {

class Simulation {
public:
    class StateInterface {
    public:
        StateInterface(Simulation& sim) : sim_(sim) {}
        const kn::spatial::UniformGrid& electron_density() const { return sim_.electron_density_; }
        const kn::spatial::UniformGrid& ion_density() const { return sim_.ion_density_; }
        const kn::particle::ChargedSpecies<1, 3>& ions() const { return sim_.ions_; }
        const kn::particle::ChargedSpecies<1, 3>& electrons() const { return sim_.electrons_; }

        const Parameters& parameters() const { return sim_.parameters_; }

        size_t step() const { return sim_.step; }

    private:
        Simulation& sim_;
    };

    friend StateInterface;

    explicit Simulation(const Parameters& parameters, const std::string& data_path);

    void run();

    enum class Event { Start, Step, End };

    struct EventAction {
        virtual void notify(const StateInterface&) = 0;
        virtual ~EventAction() {}
    };

    Events<Event, EventAction>& events();
    StateInterface& state() { return state_; };

private:
    Parameters parameters_;
    std::string data_path_;
    StateInterface state_;

    size_t step = 0;
    kn::particle::ChargedSpecies<1, 3> ions_;
    kn::particle::ChargedSpecies<1, 3> electrons_;

    kn::spatial::UniformGrid electron_density_;
    kn::spatial::UniformGrid ion_density_;

    kn::spatial::UniformGrid rho_field_;
    kn::spatial::UniformGrid phi_field_;
    kn::spatial::UniformGrid electric_field_;

    Events<Event, EventAction> events_;

    void set_initial_conditions();
    kn::collisions::MCCReactionSet<1, 3> load_electron_collisions();
    kn::collisions::MCCReactionSet<1, 3> load_ion_collisions();
};
}  // namespace ccp

#endif  // SIMULATION_H
