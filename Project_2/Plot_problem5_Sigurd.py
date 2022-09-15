import numpy as np
import matplotlib.pyplot as plt

main = np.loadtxt("run_times_problem_5_not_dense.txt", delimiter=", ", usecols=range(2))
plt.plot(
    main[:, 0],
    main[:, 1],
    color="midnightblue",
    marker="o",
    markerfacecolor="red",
    markeredgecolor="red",
)
# plt.xscale("log")
plt.grid()
plt.title("Run times for Jacobi's algorithm for non-dense matricies of size N")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.xticks([2**i for i in range(2, 8)])
plt.show()
