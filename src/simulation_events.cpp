
#include "simulation_events.h"

namespace ccp {
void setup_events(Simulation& simulation) {
    struct PrintStartAction : public Simulation::EventAction {
        int i = 0;

        void notify(const Simulation::State&) override { printf("Starting simulation (%d)\n", i); }
    };

    simulation.events().add_action<PrintStartAction>(Simulation::Event::Start);

    struct PrintStepAction : public Simulation::EventAction {
        void notify(const Simulation::State& s) override {
            if (s.step % 1000 == 0)
                printf("step: %lu\n", s.step);
        }
    };

    simulation.events().add_action<PrintStepAction>(Simulation::Event::Step);

    struct AverageFieldAction : public Simulation::EventAction {
        kn::spatial::AverageGrid av_electron_density;
        kn::spatial::AverageGrid av_ion_density;
        Parameters parameters_;

        AverageFieldAction(const Simulation::State& s, const Parameters& parameters)
            : parameters_(parameters) {
            av_electron_density = kn::spatial::AverageGrid(s.electron_density_);
            av_ion_density = kn::spatial::AverageGrid(s.ion_density_);
        }

        void notify(const Simulation::State& s) override {
            if (s.step > parameters_.n_steps - parameters_.n_steps_avg) {
                av_electron_density.add(s.electron_density_);
                av_ion_density.add(s.ion_density_);
            }
        }
    };

    auto avg_field_action = simulation.events().add_action(Simulation::Event::Step,
                                   AverageFieldAction(simulation.state(), simulation.parameters()));

    struct SaveDataAction : public Simulation::EventAction {
        std::weak_ptr<AverageFieldAction> avg_field_action_;
        explicit SaveDataAction(const std::weak_ptr<AverageFieldAction>& avg_field_action)
            : avg_field_action_(avg_field_action) {}

        void notify(const Simulation::State& s) override {
            if (!avg_field_action_.expired()) {
                //save
            }
        }
    };

    simulation.events().add_action(Simulation::Event::End, SaveDataAction(avg_field_action));
}
}  // namespace ccp