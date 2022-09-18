import numpy as np
import matplotlib.pyplot as plt

plt.subplot(1, 2, 1)

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
plt.title("Run times non-dense matricies of size N")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.ylim([10 ** (-6), 10])
plt.xticks([2**i for i in range(2, 8)])

plt.subplot(1, 2, 2)

main = np.loadtxt("run_times_problem_5_dense.txt", delimiter=", ", usecols=range(2))
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
plt.title("Run times dense matricies of size N")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.ylim([10 ** (-6), 10])
plt.xticks([2**i for i in range(2, 8)])
plt.show()
