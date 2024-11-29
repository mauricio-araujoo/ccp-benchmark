
#ifndef REACTIONS_H
#define REACTIONS_H

#include <filesystem>

#include "spark/collisions/reaction.h"
#include "parameters.h"

namespace ccp::reactions {
spark::collisions::Reactions<1, 3> load_electron_reactions(const std::filesystem::path& dir,
                                                        const Parameters& par,
                                                        spark::particle::ChargedSpecies<1, 3>& ions);

spark::collisions::Reactions<1, 3> load_ion_reactions(const std::filesystem::path& dir,
                                                   const Parameters& par);
}  // namespace ccp::reactions

#endif  // REACTIONS_H
