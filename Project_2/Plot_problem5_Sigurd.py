import numpy as np
import matplotlib.pyplot as plt

main = np.loadtxt("run_times_problem_5_not_dense.txt", delimiter=",", usecols=range(2))

fig = plt.figure()
fig.set_size_inches(w=3, h=4)

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
plt.title("Run times dense matricies of size N")
plt.xlabel("N")
# plt.ylabel("Runtime (s)")
plt.yscale("log")
plt.ylim([10 ** (-6), 10])
plt.xticks([2**i for i in range(2, 8)])
plt.savefig("plot_problem_5_dense.pgf")
