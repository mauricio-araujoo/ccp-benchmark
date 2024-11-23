import matplotlib.pyplot as plt
import numpy as np

L = 6.7e-2

de = np.genfromtxt("../cmake-build-release/density_e.txt")
di = np.genfromtxt("../cmake-build-release/density_i.txt")
x = np.linspace(0, L, de.size)

data_benchmark_1 = np.genfromtxt("../data/Benchmark_A.csv", delimiter=' ')
xb = data_benchmark_1[:, 0]
density_e_b = data_benchmark_1[:, 1]
density_i_b = data_benchmark_1[:, 4]

plt.plot(x, de, label='electrons')
plt.plot(xb, density_e_b, '--')

plt.plot(x, di, label='ions')
plt.plot(xb, density_i_b, '--')

plt.legend()
plt.show()
