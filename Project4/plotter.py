import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["figure.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["legend.fontsize"] = "large"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 16


def plot_energy_per_spin(filename, l, runs, random, T):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []
    N = l * l

    for length in simulations:
        cur = main[: length - 1]
        results.append(np.sum(cur / length / N))

    title = f"Approximation of expected energy per spin in a {l} x {l} grid"
    title += " with random initial state" * random + " with ordered inital state" * (
        not random
    )
    title += f" and T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Expected energy per spin")
    plt.plot(
        simulations, results, label="Expected energy per spin", color="midnightblue"
    )
    plt.legend()
    plt.xscale("log")
    plt.show()


def plot_magnetisation_per_spin(filename, l, runs, random, T):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[: length - 1]
        results.append(np.sum(np.abs(cur)) / length / N)

    title = f"Approximation of expected magnetisation per spin in a {l} x {l} grid"
    title += " with random initial state" * random + " with ordered inital state" * (
        not random
    )
    title += f" and T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Expected magnetisation per spin")
    plt.plot(
        simulations,
        results,
        label="Expected magnetisation per spin",
        color="midnightblue",
    )
    plt.legend()
    plt.xscale("log")
    plt.show()


def plot_specific_heat_capacity(filename, l, runs, T, random):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[: length - 1]
        E_squared = np.sum(np.power(cur, 2)) / length
        E = (np.sum(cur) / length) ** 2
        results.append((E_squared - E) / N / T / T)

    title = f"Approximation of specific heat capacity in a {l} x {l} grid"
    title += " with random initial state" * random + " with ordered inital state" * (
        not random
    )
    title += f" and T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Specific heat capacity")
    plt.plot(simulations, results, label="Specific heat capacity", color="midnightblue")
    plt.legend()
    plt.xscale("log")
    plt.show()


def plot_susceptibility(filename, l, runs, T, random):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[: length - 1]
        M_squared = np.sum(np.power(cur, 2)) / length
        M = (np.sum(np.abs(cur)) / length) ** 2
        results.append((M_squared - M) / N / T / T)

    title = f"Approximation of susceptibility in a {l} x {l} grid"
    title += " with random initial state" * random + " with ordered inital state" * (
        not random
    )
    title += f" and T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Susceptibility")
    plt.plot(simulations, results, label="Susceptibility", color="midnightblue")
    plt.legend()
    plt.xscale("log")
    plt.show()


def generate_historgram(filename, l):
    N = l * l
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    main /= N
    plt.hist(main, density=True, bins="auto")
    plt.yscale("log")
    plt.show()


if __name__ == "__main__":

    # Problem 4
    # plot_energy_per_spin("Energy_states_l_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1)
    # plot_magnetisation_per_spin(
    #     "Magnetism_states_l_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # plot_specific_heat_capacity(
    #     "Energy_states_l_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # plot_susceptibility(
    #     "Magnetism_states_l_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )

    # Problem 5
    # plot_energy_per_spin(
    #     "Energy_states_l_20_T_1.000000_1000000_random.txt", 20, 7, 1, 1
    # )
    # plot_energy_per_spin("Energy_states_l_20_T_1_1000000_positiv.txt", 20, 7, 0, 1)
    # plot_energy_per_spin("Energy_states_l_20_T_2_1000000_random.txt", 20, 7, 1, 1)
    # plot_energy_per_spin("Energy_states_l_20_T_2_1000000_positiv.txt", 20, 7, 0, 1)

    # plot_magnetisation_per_spin("Magnestism_states_l_20_T_1_1000000_random.txt", 20, 7, 1, 1)
    # plot_magnetisation_per_spin("Magnestism_states_l_20_T_1_1000000_positiv.txt", 20, 7, 0, 1)
    # plot_magnetisation_per_spin("Magnestism_states_l_20_T_2_1000000_random.txt", 20, 7, 1, 1)
    # plot_magnetisation_per_spin("Magnestism_states_l_20_T_2_1000000_positiv.txt", 20, 7, 0, 1)

    # Problem 6
    generate_historgram("Energy_states_l_20_T_1.000000_1000000_random.txt", 20)
    # generate_historgram("Energy_states_l_20_T_2.400000_1000000_random.txt", 20)
