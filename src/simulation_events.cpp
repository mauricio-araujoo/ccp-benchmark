
#include "simulation_events.h"

#include <chrono>
#include <fstream>
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
    constexpr size_t print_step_interval = 1000;

    struct PrintStartAction : public Simulation::EventAction {
        void notify(const Simulation::StateInterface&) override { printf("Starting simulation\n"); }
    };

    simulation.events().add_action<PrintStartAction>(Simulation::Event::Start);

    struct PrintEvolutionAction : public Simulation::EventAction {
        typedef std::chrono::high_resolution_clock clk;
        typedef std::chrono::duration<double, std::milli> ms;
        std::chrono::time_point<std::chrono::steady_clock> t_last;
        size_t initial_step = 0;

        void notify(const Simulation::StateInterface& s) override {
            auto step = s.step();

            if (step == 0)
                t_last = clk::now();

            if ((step % (print_step_interval / 10) == 0) && (step > 0)) {
                printf("-");
            }

            if ((step % print_step_interval == 0) && (step > 0)) {
                printf("\n");

                const auto now = clk::now();
                const double dur = std::chrono::duration_cast<ms>(now - t_last).count() /
                                   static_cast<double>(s.step() - initial_step);
                t_last = now;
                initial_step = step;

                const float progress =
                    static_cast<float>(step) /
                    static_cast<float>(std::max(1ul, s.parameters().n_steps - 1));

                const double dur_per_particle =
                    dur / (static_cast<double>(s.electrons().n() + s.ions().n()));

                printf("Info (Step: %zu/%zu, %.2f%%):\n", step, s.parameters().n_steps,
                       progress * 100.0);
                printf("    Avg step duration: %.2fms (%.2eus/p)\n", dur, dur_per_particle * 1e3);
                printf("    Sim electrons: %zu\n", s.electrons().n());
                printf("    Sim ions: %zu\n", s.ions().n());
                printf("\n");
            }
        }
    };

    simulation.events().add_action<PrintEvolutionAction>(Simulation::Event::Step);

    struct AverageFieldAction : public Simulation::EventAction {
        spark::spatial::AverageGrid av_electron_density;
        spark::spatial::AverageGrid av_ion_density;
        Parameters parameters_;

        explicit AverageFieldAction(const Parameters& parameters) : parameters_(parameters) {
            av_electron_density = spark::spatial::AverageGrid(parameters_.l, parameters_.nx);
            av_ion_density = spark::spatial::AverageGrid(parameters_.l, parameters_.nx);
        }

        void notify(const Simulation::StateInterface& s) override {
            if (s.step() > parameters_.n_steps - parameters_.n_steps_avg) {
                av_electron_density.add(s.electron_density());
                av_ion_density.add(s.ion_density());
            }
        }
    };

    auto avg_field_action = simulation.events().add_action(
        Simulation::Event::Step, AverageFieldAction(simulation.state().parameters()));

    struct SaveDataAction : public Simulation::EventAction {
        std::weak_ptr<AverageFieldAction> avg_field_action_;
        Parameters parameters_;

        explicit SaveDataAction(const std::weak_ptr<AverageFieldAction>& avg_field_action,
                                const Parameters& parameters)
            : avg_field_action_(avg_field_action), parameters_(parameters) {}

        void notify(const Simulation::StateInterface& s) override {
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

    simulation.events().add_action(
        Simulation::Event::End, SaveDataAction(avg_field_action, simulation.state().parameters()));
}
}  // namespace ccp