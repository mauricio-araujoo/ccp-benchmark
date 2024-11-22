
#ifndef REACTIONS_H
#define REACTIONS_H

#include <filesystem>

#include "kn/collisions/reaction.h"
#include "parameters.h"

namespace ccp::reactions {
kn::collisions::Reactions<1, 3> load_electron_reactions(const std::filesystem::path& dir,
                                                        const Parameters& par,
                                                        kn::particle::ChargedSpecies<1, 3>& ions);

kn::collisions::Reactions<1, 3> load_ion_reactions(const std::filesystem::path& dir,
                                                   const Parameters& par);
}  // namespace ccp::reactions

#endif  // REACTIONS_H
