
#include "reactions.h"
#include "kn/collisions/reactions/he_reactions.h"
#include "rapidcsv.h"

namespace {
    kn::collisions::CrossSection load_cross_section(const std::string &path, double energy_threshold) {
        kn::collisions::CrossSection cs;
        rapidcsv::Document doc(path, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams(';'));
        cs.energy = doc.GetColumn<double>(0);
        cs.cross_section = doc.GetColumn<double>(1);
        cs.threshold = energy_threshold;
        return cs;
    }
}

kn::collisions::Reactions<1, 3> ccp::reactions::load_electron_reactions(const std::filesystem::path &dir,
                                                                        const Parameters &par,
                                                                        kn::particle::ChargedSpecies<1, 3> &ions) {
    kn::collisions::Reactions<1, 3> electron_reactions;
    electron_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeElectronElasticCollision<1, 3> >(
            kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Elastic_He.csv", 0.0)));

    electron_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeExcitationCollision<1, 3> >(
            kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Excitation1_He.csv", 19.82)));

    electron_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeExcitationCollision<1, 3> >(
            kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Excitation2_He.csv", 20.61)));

    electron_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeIonizationCollision<1, 3> >(
            ions, par.tg, kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Ionization_He.csv", 24.59)));

    return electron_reactions;
}

kn::collisions::Reactions<1, 3> ccp::reactions::load_ion_reactions(const std::filesystem::path &dir,
                                                                   const Parameters &par) {
    kn::collisions::Reactions<1, 3> ion_reactions;
    ion_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeIonElasticCollision<1, 3> >(
            kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Isotropic_He.csv", 0.0)));

    ion_reactions.push_back(
        std::make_unique<kn::collisions::reactions::HeIonChargeExchangeCollision<1, 3> >(
            kn::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Backscattering_He.csv", 0.0)));

    return ion_reactions;
}
