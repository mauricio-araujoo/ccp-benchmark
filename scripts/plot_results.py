import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    prog='plot_results',
    description='Plot benchmark and calculation results')

parser.add_argument("out",
                    help="Path to folder containing the calculation data.")

parser.add_argument("-d", "--data_path", help="Path to folder with benchmark data",
                    default="../data")

args = parser.parse_args()

L = 6.7e-2

de = np.genfromtxt(f"{args.out}/density_e.txt")
di = np.genfromtxt(f"{args.out}/density_i.txt")
x = np.linspace(0, L, de.size)

data_benchmark_1 = np.genfromtxt(f"{args.data_path}/Benchmark_A.csv", delimiter=' ')
xb = data_benchmark_1[:, 0]
density_e_b = data_benchmark_1[:, 1]
density_i_b = data_benchmark_1[:, 4]

plt.plot(x, de, label='electrons')
plt.plot(xb, density_e_b, '--')

plt.plot(x, di, label='ions')
plt.plot(xb, density_i_b, '--')

plt.legend()
plt.show()
