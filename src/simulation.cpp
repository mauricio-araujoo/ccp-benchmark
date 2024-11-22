#include "simulation.h"

#include <kn/collisions/mcc.h>
#include <kn/constants/constants.h>
#include <kn/random/random.h>

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
}

namespace ccp {
    void Simulation::run() {

        auto electrons = kn::particle::ChargedSpecies<1, 3>(-kn::constants::e, kn::constants::m_e);
        electrons.add(parameters_.n_initial, maxwellian_emitter(parameters_.te, parameters_.l, kn::constants::m_e));

        auto ions = kn::particle::ChargedSpecies<1, 3>(kn::constants::e, parameters_.m_he);
        ions.add(parameters_.n_initial, maxwellian_emitter(parameters_.ti, parameters_.l, parameters_.m_he));

        // Load electron reactions
        auto electron_reactions = reactions::load_electron_reactions(data_path_, parameters_, ions);
        kn::collisions::ReactionConfig<1, 3> electron_reaction_config{
            parameters_.dt, parameters_.dx,
            std::make_unique<kn::collisions::StaticUniformTarget<1, 3> >(parameters_.ng, parameters_.tg),
            std::move(electron_reactions), kn::collisions::RelativeDynamics::FastProjectile
        };
        auto electron_collisions =
                kn::collisions::MCCReactionSet(electrons, std::move(electron_reaction_config));

        // Load ion reactions
        auto ion_reactions = reactions::load_ion_reactions(data_path_, parameters_);
        kn::collisions::ReactionConfig<1, 3> ion_reaction_config{
            parameters_.dt,  parameters_.dx, std::make_unique<kn::collisions::StaticUniformTarget<1, 3>>( parameters_.ng,  parameters_.tg),
            std::move(ion_reactions), kn::collisions::RelativeDynamics::SlowProjectile};
        auto ion_collisions = kn::collisions::MCCReactionSet(ions, std::move(ion_reaction_config));
    }
} // ccp
