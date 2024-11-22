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
    return [l, t, m](kn::core::Vec<3>& v, kn::core::Vec<1>& x) {
        x.x = l * kn::random::uniform();
        double vth = std::sqrt(kn::constants::kb * t / m);
        v = {kn::random::normal(0.0, vth), kn::random::normal(0.0, vth),
             kn::random::normal(0.0, vth)};
    };
}
}  // namespace

namespace ccp {
Simulation::Simulation(const Parameters& parameters, const std::string& data_path)
    : parameters_(parameters), data_path_(data_path) {}

void Simulation::run() {
    set_initial_conditions();

    auto electron_collisions = load_electron_collisions();
    auto ion_collisions = load_ion_collisions();

    auto poisson_solver =
        kn::electromagnetics::DirichletPoissonSolver(parameters_.nx, parameters_.dx);

    events().notify(Event::Start, state_);

    for (state_.step = 0; state_.step < parameters_.n_steps; ++state_.step) {
        kn::interpolate::weight_to_grid(state_.electrons_, state_.electron_density_);
        kn::interpolate::weight_to_grid(state_.ions_, state_.ion_density_);
        kn::electromagnetics::charge_density(parameters_.particle_weight, state_.ion_density_,
                                             state_.electron_density_, state_.rho_field_);

        const double boundary_voltage =
            parameters_.volt * std::sin(2.0 * kn::constants::pi * parameters_.f * parameters_.dt *
                                        static_cast<double>(state_.step));
        poisson_solver.solve(state_.rho_field_.data(), state_.phi_field_.data(), 0.0,
                             boundary_voltage);
        poisson_solver.efield(state_.phi_field_.data(), state_.electric_field_.data());

        kn::interpolate::field_at_particles(state_.electric_field_, state_.electrons_);
        kn::interpolate::field_at_particles(state_.electric_field_, state_.ions_);

        kn::particle::move_particles(state_.electrons_, parameters_.dt);
        kn::particle::move_particles(state_.ions_, parameters_.dt);

        kn::particle::apply_absorbing_boundary(state_.electrons_, 0, parameters_.l);
        kn::particle::apply_absorbing_boundary(state_.ions_, 0, parameters_.l);

        electron_collisions.react_all();
        ion_collisions.react_all();

        events().notify(Event::Step, state_);
    }

    events().notify(Event::End, state_);
}

const Simulation::State& Simulation::state() const {
    return state_;
}

const Parameters& Simulation::parameters() const {
    return parameters_;
}

Events<Simulation::Event, Simulation::EventAction>& Simulation::events() {
    return events_;
}

void Simulation::set_initial_conditions() {
    // Charged species
    state_.electrons_ = kn::particle::ChargedSpecies<1, 3>(-kn::constants::e, kn::constants::m_e);
    state_.electrons_.add(parameters_.n_initial,
                          maxwellian_emitter(parameters_.te, parameters_.l, kn::constants::m_e));

    state_.ions_ = kn::particle::ChargedSpecies<1, 3>(kn::constants::e, parameters_.m_he);
    state_.ions_.add(parameters_.n_initial,
                     maxwellian_emitter(parameters_.ti, parameters_.l, parameters_.m_he));

    // Fields
    state_.electron_density_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
    state_.ion_density_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
    state_.rho_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
    state_.phi_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
    state_.electric_field_ = kn::spatial::UniformGrid(parameters_.l, parameters_.nx);
}

kn::collisions::MCCReactionSet<1, 3> Simulation::load_electron_collisions() {
    // Load electron reactions
    auto electron_reactions =
        reactions::load_electron_reactions(data_path_, parameters_, state_.ions_);
    kn::collisions::ReactionConfig<1, 3> electron_reaction_config{
        parameters_.dt, parameters_.dx,
        std::make_unique<kn::collisions::StaticUniformTarget<1, 3>>(parameters_.ng, parameters_.tg),
        std::move(electron_reactions), kn::collisions::RelativeDynamics::FastProjectile};

    return kn::collisions::MCCReactionSet(state_.electrons_, std::move(electron_reaction_config));
}

kn::collisions::MCCReactionSet<1, 3> Simulation::load_ion_collisions() {
    // Load ion reactions
    auto ion_reactions = reactions::load_ion_reactions(data_path_, parameters_);
    kn::collisions::ReactionConfig<1, 3> ion_reaction_config{
        parameters_.dt, parameters_.dx,
        std::make_unique<kn::collisions::StaticUniformTarget<1, 3>>(parameters_.ng, parameters_.tg),
        std::move(ion_reactions), kn::collisions::RelativeDynamics::SlowProjectile};

    return kn::collisions::MCCReactionSet(state_.ions_, std::move(ion_reaction_config));
}
}  // namespace ccp
