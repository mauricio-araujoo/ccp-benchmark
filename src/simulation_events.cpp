
#include "simulation_events.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <span>

namespace {
template <class It>
void save_vec(const char* filename, const It& vec) {
    std::ofstream out_file(filename);

    size_t i = 0;
    for (const auto& v : vec) {
        if (i != 0) {
            out_file << "\n";
        }
        out_file << v;
        i++;
    }
}

std::vector<double> count_to_density(double particle_weight,
                                     double dx,
                                     const std::vector<double>& count) {
    auto d = std::vector<double>(count.size());
    std::ranges::transform(count, d.begin(), [particle_weight, dx](const double val) {
        return val * particle_weight / dx;
    });
    return d;
}
}  // namespace

namespace ccp {
void setup_events(Simulation& simulation) {
    struct PrintStartAction : public Simulation::EventAction {
        void notify(const Simulation::State&) override { printf("Starting simulation\n"); }
    };

    simulation.events().add_action<PrintStartAction>(Simulation::Event::Start);

    struct MeasureTimeStepAction : public Simulation::EventAction {
        typedef std::chrono::high_resolution_clock clk;
        typedef std::chrono::milliseconds ms;

        std::chrono::time_point<std::chrono::steady_clock> t_last;
        std::chrono::duration<double> mean_duration{0.0};
        size_t initial_step = 0;

        void notify(const Simulation::State& s) override {
            if (initial_step == 0)
                t_last = clk::now();

            if ((s.step % 1000 == 0) && s.step > 0) {
                const auto now = clk::now();
                const double dur = std::chrono::duration_cast<ms>(now - t_last).count() /
                                   static_cast<double>(s.step - initial_step);
                printf("Average step duration: %.2f ms\n", dur);
                t_last = now;
                initial_step = s.step;
            }
        }
    };

    simulation.events().add_action<MeasureTimeStepAction>(Simulation::Event::Step);

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
        Parameters parameters_;

        explicit SaveDataAction(const std::weak_ptr<AverageFieldAction>& avg_field_action,
                                const Parameters& parameters)
            : avg_field_action_(avg_field_action), parameters_(parameters) {}

        void notify(const Simulation::State& s) override {
            if (!avg_field_action_.expired()) {
                const auto avg_field_action_ptr = avg_field_action_.lock();
                const auto& avg_e = avg_field_action_ptr->av_electron_density.get();
                const auto& avg_i = avg_field_action_ptr->av_ion_density.get();

                save_vec("density_e.txt",
                         count_to_density(parameters_.particle_weight, parameters_.dx, avg_e));
                save_vec("density_i.txt",
                         count_to_density(parameters_.particle_weight, parameters_.dx, avg_i));
            }
        }
    };

    simulation.events().add_action(Simulation::Event::End,
                                   SaveDataAction(avg_field_action, simulation.parameters()));
}
}  // namespace ccp