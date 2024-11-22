
#include "simulation_events.h"

#include <fstream>
#include <span>

namespace {
void save_vec(const char* filename, const std::span<double>& v) {
    std::ofstream out_file(filename);

    for (size_t i = 0; i < v.size(); i++) {
        if (i != 0) {
            out_file << "\n";
        }
        out_file << v[i];
    }
}
}  // namespace

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

        AverageFieldAction(const Parameters& parameters) : parameters_(parameters) {
            av_electron_density = kn::spatial::AverageGrid(parameters_.l, parameters_.nx);
            av_ion_density = kn::spatial::AverageGrid(parameters_.l, parameters_.nx);
        }

        void notify(const Simulation::State& s) override {
            if (s.step > parameters_.n_steps - parameters_.n_steps_avg) {
                av_electron_density.add(s.electron_density_);
                av_ion_density.add(s.ion_density_);
            }
        }
    };

    auto avg_field_action = simulation.events().add_action(
        Simulation::Event::Step, AverageFieldAction(simulation.parameters()));

    struct SaveDataAction : public Simulation::EventAction {
        std::weak_ptr<AverageFieldAction> avg_field_action_;
        explicit SaveDataAction(const std::weak_ptr<AverageFieldAction>& avg_field_action)
            : avg_field_action_(avg_field_action) {}

        void notify(const Simulation::State& s) override {
            if (!avg_field_action_.expired()) {
                auto avg_field_action_ptr = avg_field_action_.lock();
                save_vec("density_e.txt", avg_field_action_ptr->av_electron_density.get());
                save_vec("density_i.txt", avg_field_action_ptr->av_ion_density.get());
            }
        }
    };

    simulation.events().add_action(Simulation::Event::End, SaveDataAction(avg_field_action));
}
}  // namespace ccp