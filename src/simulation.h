#ifndef SIMULATION_H
#define SIMULATION_H

#include <spark/collisions/mcc.h>
#include <spark/particle/species.h>
#include <spark/spatial/grid.h>

#include <string>

#include "events.h"
#include "parameters.h"

namespace ccp {

class Simulation {
public:
    class StateInterface {
    public:
        StateInterface(Simulation& sim) : sim_(sim) {}
        const spark::spatial::UniformGrid<1>& electron_density() const {
            return sim_.electron_density_;
        }
        const spark::spatial::UniformGrid<1>& ion_density() const { return sim_.ion_density_; }
        const spark::particle::ChargedSpecies<1, 3>& ions() const { return sim_.ions_; }
        const spark::particle::ChargedSpecies<1, 3>& electrons() const { return sim_.electrons_; }

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
    spark::particle::ChargedSpecies<1, 3> ions_;
    spark::particle::ChargedSpecies<1, 3> electrons_;

    spark::spatial::UniformGrid<1> electron_density_;
    spark::spatial::UniformGrid<1> ion_density_;

    spark::spatial::UniformGrid<1> rho_field_;
    spark::spatial::UniformGrid<1> phi_field_;
    spark::spatial::TUniformGrid<spark::core::Vec<1>, 1> electric_field_;

    Events<Event, EventAction> events_;

    void set_initial_conditions();
    spark::collisions::MCCReactionSet<1, 3> load_electron_collisions();
    spark::collisions::MCCReactionSet<1, 3> load_ion_collisions();
};
}  // namespace ccp

#endif  // SIMULATION_H
