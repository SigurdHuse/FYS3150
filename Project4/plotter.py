import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["figure.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams["legend.fontsize"] = "medium"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 16


def plot_energy_per_spin(filename, l, runs, random, T):
    """Script to plot energy per spin from file"""

    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []
    N = l * l

    for length in simulations:
        cur = main[:length]
        results.append(np.sum(cur / length / N))

    title = f"Approximation of expected energy per spin in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Expected energy per spin [J]")
    labelname = "Expected energy per spin"
    labelname += (
        " with random initial state" * random
        + " with ordered inital state" * (not random)
    )

    color_rand = "midnightblue" * random + (not random) * "red"

    plt.plot(simulations, results, label=labelname, color=color_rand)
    plt.legend()
    plt.xscale("log")


def plot_magnetisation_per_spin(filename, l, runs, random, T):
    """Script to plot magnetisation per spin from file"""

    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[:length]
        results.append(np.sum(np.abs(cur)) / length / N)

    title = f"Approximation of expected magnetisation per spin in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Expected magnetisation per spin")

    labelname = "Expected magnetisation per spin"
    labelname += (
        " with random initial state" * random
        + " with ordered inital state" * (not random)
    )

    color_rand = "midnightblue" * random + (not random) * "red"

    plt.plot(
        simulations,
        results,
        label=labelname,
        color=color_rand,
    )
    plt.legend()
    plt.xscale("log")


def plot_specific_heat_capacity(filename, l, runs, T, random):
    """Script to plot specific heat capacity per spin from file"""

    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[:length]
        E_squared = np.sum(np.power(cur, 2)) / length
        E = (np.sum(cur) / length) ** 2
        results.append((E_squared - E) / N / T / T)

    title = f"Approximation of specific heat capacity in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    labelname = "Specific heat capacity"
    labelname += (
        " with random initial state" * random
        + " with ordered inital state" * (not random)
    )

    color_rand = "midnightblue" * random + (not random) * "red"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Specific heat capacity")
    plt.plot(simulations, results, label=labelname, color=color_rand)
    plt.legend()
    plt.xscale("log")


def plot_susceptibility(filename, l, runs, T, random):
    """Script to plot susceptibility per spin from file"""

    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []

    N = l * l

    for length in simulations:
        cur = main[:length]
        M_squared = np.sum(np.power(cur, 2)) / length
        M = (np.sum(np.abs(cur)) / length) ** 2
        results.append((M_squared - M) / N / T / T)

    title = f"Approximation of susceptibility in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    labelname = "Susceptibility"
    labelname += (
        " with random initial state" * random
        + " with ordered inital state" * (not random)
    )

    color_rand = "midnightblue" * random + (not random) * "red"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles")
    plt.ylabel("Susceptibility")
    plt.plot(simulations, results, label=labelname, color=color_rand)
    plt.legend()
    plt.xscale("log")


def generate_historgram(filename, l, T, bins, log):
    """Script to plot distributions of energy per spin from file"""

    N = l * l
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    main /= N
    plt.hist(main, density=True, bins=bins)
    if log == 1:
        plt.yscale("log")
    plt.grid()
    plt.title(f"Estimated probability distribution of $\epsilon$ when T = {T}, L = {l}")
    plt.ylabel("Probability")
    plt.xlabel("Energy per spin")


def plot_by_temp():
    """Script to plot quantities as a function of temperature from files"""

    sizes = [40, 60, 80, 100]
    temps = np.arange(2.1, 2.41, 0.02)
    runs = 10000000
    file_runs = 1000000
    energy = []
    magnet = []
    specific = []
    susc = []

    for L in sizes:
        tmp_energy = []
        tmp_magnet = []
        tmp_specific = []
        tmp_susc = []
        N = L * L
        for T in temps:
            T_squared = T * T
            filename1 = (
                f"data{L}/Energy_states_L_{L}_T_{T :.2f}0000_{file_runs}_random.txt"
            )
            filename2 = (
                f"data{L}/Magnetism_states_L_{L}_T_{T :.2f}0000_{file_runs}_random.txt"
            )
            main_E = np.loadtxt(filename1, delimiter=",", usecols=range(1), skiprows=1)
            main_M = np.loadtxt(filename2, delimiter=",", usecols=range(1), skiprows=1)

            m = np.sum(np.abs(main_M)) / N / runs
            epsilon = np.sum(main_E / N / runs)

            E_squared = np.sum(np.power(main_E, 2)) / runs
            E = (epsilon * L * L) ** 2
            spec = (E_squared - E) / N / T_squared

            M_squared = np.sum(np.power(main_M, 2)) / runs
            M = (m * L * L) ** 2
            suspect = (M_squared - M) / N / T_squared

            tmp_energy.append(epsilon)
            tmp_magnet.append(m)
            tmp_specific.append(spec)
            tmp_susc.append(suspect)

        energy.append(tmp_energy)
        magnet.append(tmp_magnet)
        specific.append(tmp_specific)
        susc.append(tmp_susc)

    for i in range(len(energy)):
        plt.plot(temps, energy[i], label=f"L = {sizes[i]}")
    plt.legend()
    plt.grid()
    plt.title("Energy")
    plt.show()
    plt.clf()

    for i in range(len(magnet)):
        plt.plot(temps, magnet[i], label=f"L = {sizes[i]}")
    plt.legend()
    plt.grid()
    plt.title("Magnet")
    plt.show()
    plt.clf()

    for i in range(len(specific)):
        plt.plot(temps, specific[i], label=f"L = {sizes[i]}")
    plt.legend()
    plt.show()
    plt.title("Specific")
    plt.grid()
    plt.clf()

    for i in range(len(susc)):
        plt.plot(temps, susc[i], label=f"L = {sizes[i]}")
    plt.legend()
    plt.show()
    plt.title("Susceptibility")
    plt.grid()
    plt.clf()


if __name__ == "__main__":
    plt.rc("pgf", texsystem="pdflatex")
    # Problem 4
    # plot_energy_per_spin(
    #     "data_5_6/Energy_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # # plt.show()
    # plt.savefig("plots/Energy_L_2_T_1.pgf")
    # plt.clf()

    # plot_magnetisation_per_spin(
    #     "data_5_6/Magnetism_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # plt.savefig("plots/Magnet_L_2_T_1.pgf")
    # plt.clf()

    # plot_specific_heat_capacity(
    #     "data_5_6/Energy_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # plt.savefig("plots/HC_L_2_T_1.pgf")
    # plt.clf()

    # plot_susceptibility(
    #     "data_5_6/Magnetism_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    # )
    # plt.savefig("plots/susceptibility_L_2_T_1.pgf")
    # plt.clf()

    # Problem 5
    plot_energy_per_spin(
        "data_5_6/Energy_states_L_20_T_1.000000_1000000_random.txt", 20, 7, 1, 1
    )

    plot_energy_per_spin(
        "data_5_6/Energy_states_L_20_T_1.000000_1000000_positiv.txt", 20, 7, 0, 1
    )
    plt.grid()
    plt.savefig("plots/Energy_L_20_T_1.pgf")
    plt.clf()

    plot_energy_per_spin(
        "data_5_6/Energy_states_L_20_T_2.400000_1000000_random.txt", 20, 7, 1, 2.4
    )

    plot_energy_per_spin(
        "data_5_6/Energy_states_L_20_T_2.400000_1000000_positiv.txt", 20, 7, 0, 2.4
    )
    plt.grid()
    plt.savefig("plots/Energy_L_20_T_2.4.pgf")
    plt.clf()

    plot_magnetisation_per_spin(
        "data_5_6/Magnetism_states_L_20_T_1.000000_1000000_random.txt", 20, 7, 1, 1
    )
    plot_magnetisation_per_spin(
        "data_5_6/Magnetism_states_L_20_T_1.000000_1000000_positiv.txt", 20, 7, 0, 1
    )
    plt.grid()
    plt.savefig("plots/Magnet_L_20_T_1.pgf")
    plt.clf()

    plot_magnetisation_per_spin(
        "data_5_6/Magnetism_states_L_20_T_2.400000_1000000_random.txt", 20, 7, 1, 2.4
    )
    plot_magnetisation_per_spin(
        "data_5_6/Magnetism_states_L_20_T_2.400000_1000000_positiv.txt", 20, 7, 0, 2.4
    )
    plt.grid()
    plt.savefig("plots/Magnet_L_20_T_2.4.pgf")
    plt.clf()

    # Problem 6
    # generate_historgram(
    #     "data_5_6/Energy_states_L_20_T_1.000000_1000000_random.txt",
    #     20,
    #     1,
    #     np.arange(-2, -1.9, 0.005),
    #     1,
    # )
    # plt.savefig("plots/Probs_L_20_T_1.pgf")

    # generate_historgram(
    #     "data_5_6/Energy_states_L_20_T_2.400000_1000000_random.txt",
    #     20,
    #     2.4,
    #     np.arange(-1.9, -0.6, 0.005),
    #     0,
    # )
    # plt.savefig("plots/Probs_L_20_T_2.4.pgf")

    # Problem 8
    # plot_by_temp()
