import numpy as np
import matplotlib.pyplot as plt

"""Script to plot"""

main = np.loadtxt("run_times_problem_5_not_dense.txt", delimiter=",", usecols=range(2))

fig = plt.figure()
fig.set_size_inches(w=5, h=8)
# plt.tick_params(axis="y", which="both", labelleft=False, labelright=True)

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
plt.title("Run times non-dense matrices of size N")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.ylim([10 ** (-6), 10])
plt.xticks([2**i for i in range(2, 8)])
plt.rc("pgf", texsystem="pdflatex")
plt.savefig("plot_problem_5_A.pgf")

plt.clf()

main = np.loadtxt("run_times_problem_5_dense.txt", delimiter=",", usecols=range(2))
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
plt.title("Run times dense matrices of size N")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.ylim([10 ** (-6), 10])
plt.xticks([2**i for i in range(2, 8)])
plt.savefig("plot_problem_5_dense.pgf")
