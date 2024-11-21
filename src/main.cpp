#include <cstdio>
#include <string>

#include <argparse/argparse.hpp>

#include "simulation.h"
#include "kn/random/random.h"

ccp::Parameters get_case_parameters(int case_number) {

    ccp::Parameters p;

    switch (case_number) {
        case 1:
            p = ccp::Parameters::case_1();
        break;
        case 2:
            p = ccp::Parameters::case_2();
        break;
        case 3:
            p = ccp::Parameters::case_3();
        break;
        case 4:
            p = ccp::Parameters::case_4();
        break;
        default:
        break;
    }

    return p;
}

int main(int argc, char *argv[]) {
    kn::random::initialize(42);
    argparse::ArgumentParser args("cpp-benchmark");

    int case_number = 0;
    args.add_argument("case_number")
        .help("Benchmark case to be simulated")
        .scan<'i', int>()
        .default_value(1)
        .choices(1, 2, 3, 4)
        .store_into(case_number);

    std::string data_path{"../data"};
    args.add_argument("-d", "--data")
        .help("Path to folder with cross section data")
        .default_value(data_path)
        .store_into(data_path);

    args.parse_args(argc, argv);

    printf("Starting benchmark case %d simulation\n", case_number);
    printf("Data path set to %s\n", data_path.c_str());

    ccp::Simulation sim(get_case_parameters(case_number), data_path);
    sim.run();

    return 0;
}
