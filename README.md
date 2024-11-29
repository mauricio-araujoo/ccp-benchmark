# Capacitively Coupled Plasma benchmark
Implementation of the PIC capacitively coupled plasma benchmarks by 

> M. M. Turner, et al. *Simulation benchmarks for low-pressure plasmas: Capacitive discharges.* Phys. Plasmas 1 January 2013; 20 (1): 013507. [10.1063/1.4775084](https://doi.org/10.1063/1.4775084).

using the [spark](https://github.com/lase-unb/spark) particle kinetic simulation library.

## Building and running

To build the code it is necessary to have a compiler that supports C++20 and CMake. Then you can simply run from the project's root folder

```sh
mkdir build && cd build
cmake .
cmake --build .
```

To run simulation, you need to pass the number of the benchmark to be executed and, optionally, the path to the folder containing the collision cross sections. This folder is the `data` folder in the project. The usage pattern is shown below:

```sh
Usage: cpp-benchmark [--help] [--version] [--data VAR] case_number

Positional arguments:
  case_number    Benchmark case to be simulated [default: 1]

Optional arguments:
  -h, --help     shows help message and exits
  -v, --version  prints version information and exits
  -d, --data     Path to folder with cross section data [default: "../data"]
```
