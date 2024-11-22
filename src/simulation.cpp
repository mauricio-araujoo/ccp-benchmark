#include "simulation.h"

#include <kn/collisions/mcc.h>
#include <kn/constants/constants.h>
#include <kn/electromagnetics/poisson.h>
#include <kn/interpolate/field.h>
#include <kn/interpolate/weight.h>
#include <kn/particle/boundary.h>
#include <kn/particle/pusher.h>
#include <kn/random/random.h>
#include <kn/spatial/grid.h>

#include "reactions.h"

namespace {
    auto maxwellian_emitter(double t, double l, double m) {
        return [l, t, m](kn::core::Vec<3> &v, kn::core::Vec<1> &x) {
            x.x = l * kn::random::uniform();
            double vth = std::sqrt(kn::constants::kb * t / m);
            v = {
                kn::random::normal(0.0, vth), kn::random::normal(0.0, vth),
                kn::random::normal(0.0, vth)
            };
        };
    }
}

namespace ccp {
    void Simulation::run() {
        set_initial_conditions();

        auto electron_collisions = load_electron_collisions();
        auto ion_collisions = load_ion_collisions();

        auto av_electron_density = kn::spatial::AverageGrid(electron_density_);
        auto av_ion_density = kn::spatial::AverageGrid(ion_density_);

        auto poisson_solver = kn::electromagnetics::DirichletPoissonSolver(parameters_.nx, parameters_.dx);

        printf("Starting main loop...\n");
        for (size_t i = 0; i < parameters_.n_steps; ++i) {
            kn::interpolate::weight_to_grid(electrons_, electron_density_);
            kn::interpolate::weight_to_grid(ions_, ion_density_);
            kn::electromagnetics::charge_density(parameters_.particle_weight, ion_density_, electron_density_,
                                                 rho_field_);

            const double vc = parameters_.volt * std::sin(
                                  2.0 * kn::constants::pi * parameters_.f * parameters_.dt * static_cast<double>(i));
            poisson_solver.solve(rho_field_.data(), phi_field_.data(), 0.0, vc);
            poisson_solver.efield(phi_field_.data(), electric_field_.data());

            kn::interpolate::field_at_particles(electric_field_, electrons_);
            kn::interpolate::field_at_particles(electric_field_, ions_);

            kn::particle::move_particles(electrons_, parameters_.dt);
            kn::particle::move_particles(ions_, parameters_.dt);

            kn::particle::apply_absorbing_boundary(electrons_, 0, parameters_.l);
            kn::particle::apply_absorbing_boundary(ions_, 0, parameters_.l);

            electron_collisions.react_all();
            ion_collisions.react_all();

            // Generate average field
            if (i > parameters_.n_steps - parameters_.n_steps_avg) {
                av_electron_density.add(electron_density_);
                av_ion_density.add(ion_density_);
            }
        }
    }

    void Simulation::set_initial_conditions() {
        // Charged species
        electrons_ = kn::particle::ChargedSpecies<1, 3>(-kn::constants::e, kn::constants::m_e);
        electrons_.add(parameters_.n_initial, maxwellian_emitter(parameters_.te, parameters_.l, kn::constants::m_e));

        ions_ = kn::particle::ChargedSpecies<1, 3>(kn::constants::e, parameters_.m_he);
        ions_.add(parameters_.n_initial, maxwellian_emitter(parameters_.ti, parameters_.l, parameters_.m_he));

        // Fields
        electron_density_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
        ion_density_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
        rho_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
        phi_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
        electric_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
    }

    kn::collisions::MCCReactionSet<1, 3> Simulation::load_electron_collisions() {
        // Load electron reactions
        auto electron_reactions = reactions::load_electron_reactions(data_path_, parameters_, ions_);
        kn::collisions::ReactionConfig<1, 3> electron_reaction_config{
            parameters_.dt, parameters_.dx,
            std::make_unique<kn::collisions::StaticUniformTarget<1, 3> >(parameters_.ng, parameters_.tg),
            std::move(electron_reactions), kn::collisions::RelativeDynamics::FastProjectile
        };

        return kn::collisions::MCCReactionSet(electrons_, std::move(electron_reaction_config));
    }

    kn::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions() {
        // Load ion reactions
        auto ion_reactions = reactions::load_ion_reactions(data_path_, parameters_);
        kn::collisions::ReactionConfig<1, 3> ion_reaction_config{
            parameters_.dt, parameters_.dx,
            std::make_unique<kn::collisions::StaticUniformTarget<1, 3> >(parameters_.ng, parameters_.tg),
            std::move(ion_reactions), kn::collisions::RelativeDynamics::SlowProjectile
        };

        return kn::collisions::MCCReactionSet(ions_, std::move(ion_reaction_config));
    }
} // ccp
