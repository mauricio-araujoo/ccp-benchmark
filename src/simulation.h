#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <kn/collisions/mcc.h>
#include <kn/particle/species.h>
#include <kn/spatial/grid.h>

#include "parameters.h"
#include "events.h"

namespace ccp {
    class Simulation {
    public:

        explicit Simulation(const Parameters &parameters, const std::string &data_path);

        void run();

        struct State {
            size_t step = 0;
            kn::particle::ChargedSpecies<1, 3> ions_;
            kn::particle::ChargedSpecies<1, 3> electrons_;

            kn::spatial::UniformGrid electron_density_;
            kn::spatial::UniformGrid ion_density_;

            kn::spatial::UniformGrid rho_field_;
            kn::spatial::UniformGrid phi_field_;
            kn::spatial::UniformGrid electric_field_;
        };

        enum class Event {
            Start,
            Step,
            End
        };

        struct EventAction {
            virtual void notify(const State&) = 0;
            virtual ~EventAction() {}
        };

        const State& state() const;
        const Parameters& parameters() const;
        Events<Event, EventAction>& events();

    private:
        Parameters parameters_;
        std::string data_path_;

        State state_;

        Events<Event, EventAction> events_;

        void set_initial_conditions();
        kn::collisions::MCCReactionSet<1, 3> load_electron_collisions();
        kn::collisions::MCCReactionSet<1, 3> load_ion_collisions();
    };
} // ccp

#endif //SIMULATION_H
