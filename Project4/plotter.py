import numpy as np
import matplotlib.pyplot as plt


def plot_energy_per_spin(filename, N, runs):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([10**i for i in range(1, runs + 1)])
    results = []

    for length in simulations:
        cur = main[: length - 1]
        results.append(np.sum(cur) / length / N)

    plt.plot(simulations, results)
    plt.xscale("log")
    plt.show()


def plot_magnetisation_per_spin(filename, N, runs):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([10**i for i in range(1, runs + 1)])
    results = []

    for length in simulations:
        cur = main[: length - 1]
        results.append(np.sum(np.abs(cur)) / length / N)

    plt.plot(simulations, results)
    plt.xscale("log")
    plt.show()


def plot_specific_heat_capacity(filename, N, runs, T):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([10**i for i in range(1, runs + 1)])
    results = []

    for length in simulations:
        cur = main[: length - 1]
        E_squared = np.sum(np.power(cur, 2)) / length
        E = (np.sum(cur) / length) ** 2
        results.append((E_squared - E) / N / T / T)

    plt.plot(simulations, results)
    plt.xscale("log")
    plt.show()


def plot_susceptibility(filename, N, runs, T):
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([10**i for i in range(1, runs + 1)])
    results = []

    for length in simulations:
        cur = main[: length - 1]
        M_squared = np.sum(np.power(cur, 2)) / length
        M = (np.sum(np.abs(cur)) / length) ** 2
        results.append((M_squared - M) / N / T / T)

    plt.plot(simulations, results)
    plt.xscale("log")
    plt.show()


if __name__ == "__main__":
    # plot_energy_per_spin("Energy_states_l_2_T_1_100000_random.txt", 4, 6)
    # plot_magnetisation_per_spin("Magnetism_states_l_2_T_1_100000_random.txt", 4, 6)
    # plot_specific_heat_capacity("Energy_states_l_2_T_1_100000_random.txt", 4, 6, 1)
    plot_susceptibility("Magnetism_states_l_2_T_1_100000_random.txt", 4, 6, 1)
