#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstddef>

namespace ccp {

struct Parameters {

    size_t nx;
    double f;
    double dt;
    double l;
    double dx;
    double ng;
    double tg;
    double te;
    double ti;
    double n0;
    double m_he;
    double m_e;
    double volt;
    size_t ppc;
    size_t n_steps;
    size_t n_steps_avg;
    double particle_weight;

    static Parameters case_1();
    static Parameters case_2();
    static Parameters case_3();
    static Parameters case_4();

private:
    void fixed_parameters();
    void computed_parameters();
};

}
#endif //PARAMETERS_H
