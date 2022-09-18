import numpy as np
import matplotlib.pyplot as plt


n_values = [10, 100]

for n in n_values:
    filename = "values_problem_6_n_" + str(n) + ".txt"
    print(filename)
    main = np.loadtxt(filename, delimiter=",", usecols=range(4))
    colors = ["midnightblue", "red", "green"]
    for clr, i in zip(colors, range(1, 4)):
        plt.plot(main[:, 0], main[:, i], color=clr, label="beam " + str(i))
    plt.grid()
    plt.legend()
    plt.xlim([0, 1])
    plt.xlabel("$x_i$")
    plt.ylabel("$v_i$ (height of beam at point $x_i$)")
    plt.title("Plot of three first eigenvectors, when n = " + str(n))
    plt.rc("pgf", texsystem="pdflatex")  # or luatex, xelatex...
    plt.savefig("plot_problem_6_n_" + str(n) + ".pgf")
    # plt.show()
    plt.clf()
