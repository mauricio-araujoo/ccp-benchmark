
#include "parameters.h"

namespace ccp {

    void Parameters::fixed_parameters() {
        tg = 300.0;
        te = 30'000.0;
        ti = 300.0;
        m_he = 6.67e-27;
        m_e = 9.109e-31;
        l = 6.7e-2;
        f = 13.56e6;
    }

    void Parameters::computed_parameters() {
        dx = l / static_cast<double>(nx - 1);
        particle_weight = n0 * l / static_cast<double>(ppc * (nx - 1));
    }

    Parameters Parameters::case_1() {
        Parameters p{};
        p.fixed_parameters();

        p.nx = 129;
        p.dt = 1.0 / (400.0 * p.f);
        p.ng = 9.64e20;
        p.n0 = 2.56e14;
        p.volt = 450.0;
        p.ppc = 512;
        p.n_steps = 512'000;
        p.n_steps_avg = 12'800;

        p.computed_parameters();
        return p;
    }

    Parameters Parameters::case_2() {
        Parameters p{};
        p.fixed_parameters();

        p.nx = 257;
        p.dt = 1.0 / (800.0 * p.f);
        p.ng = 32.1e20;
        p.n0 = 5.12e14;
        p.volt = 200.0;
        p.ppc = 256;
        p.n_steps = 4'096'000;
        p.n_steps_avg = 25'600;

        p.computed_parameters();
        return p;
    }

    Parameters Parameters::case_3() {
        Parameters p{};
        p.fixed_parameters();

        p.nx = 513;
        p.dt = 1.0 / (1600.0 * p.f);
        p.ng = 96.4e20;
        p.n0 = 5.12e14;
        p.volt = 150.0;
        p.ppc = 128;
        p.n_steps = 8'192'000;
        p.n_steps_avg = 51'200;

        p.computed_parameters();
        return p;
    }

    Parameters Parameters::case_4() {
        Parameters p{};
        p.fixed_parameters();

        p.nx = 513;
        p.dt = 1.0 / (3200.0 * p.f);
        p.ng = 321.0e20;
        p.n0 = 3.84e14;
        p.volt = 120.0;
        p.ppc = 64;
        p.n_steps = 49'152'000;
        p.n_steps_avg = 102'400;

        p.computed_parameters();
        return p;
    }
}