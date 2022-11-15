import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
import os

mpl.rcParams["figure.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams["legend.fontsize"] = "medium"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12


def plot_energy_per_spin(filename, l, runs, random, T):
    """Script to plot energy per spin from file"""

    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    simulations = np.array([j * 10**i for i in range(runs) for j in range(2, 11)])
    results = []
    N = l * l

    for length in simulations:
        cur = main[:length]
        results.append(np.sum(cur / length / N))

    title = f"Approximation of <$\epsilon$> in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles [1]")
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

    title = f"Approximation of <$|m|$> in a {l} x {l} grid"
    # title += " with random initial state" * random + " with ordered inital state" * (
    #    not random
    # )
    title += f" with T = {T}"

    plt.grid()
    plt.title(title)
    plt.xlabel("Number of Monte Carlo cycles [1]")
    plt.ylabel("Expected absolute value of magnetisation per spin [1]")

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

    title = f"Approximation of $C_V$ in a {l} x {l} grid"
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
    plt.xlabel("Number of Monte Carlo cycles [1]")
    plt.ylabel(r"Specific heat capacity [k_B]")
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

    title = f"Approximation of $\chi$ in a {l} x {l} grid"
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
    plt.xlabel("Number of Monte Carlo cycles [1]")
    plt.ylabel("Susceptibility [1/J]")
    plt.plot(simulations, results, label=labelname, color=color_rand)
    plt.legend()
    plt.xscale("log")


def generate_historgram(filename, l, T, bins, log):
    """Script to plot distributions of energy per spin from file"""

    N = l * l
    main = np.loadtxt(filename, delimiter=",", usecols=range(1), skiprows=1)
    main /= N

    print(f"Estimated variance for a {l} x {l} grid with T = {T}:", np.var(main))

    bin_height, bin_boundary = np.histogram(main, bins=bins)
    width = bin_boundary[1] - bin_boundary[0]
    bin_height = bin_height / np.sum(bin_height)
    plt.bar(bin_boundary[:-1], bin_height, width=width)

    if log == 1:
        plt.yscale("log")

    plt.grid()
    plt.title(
        f"Estimated probability distribution of <$\epsilon$> when T = {T}, L = {l}"
    )
    plt.ylabel("Probability [1]")
    plt.xlabel("Expected energy per spin [J]")


def compute_energy(sizes, temps, file_runs, runs):
    for L in sizes:
        values = []
        N = L * L
        filename_save = f"data/Computed_Energy_states_L_{L}.txt"
        for T in temps:
            filename = f"data/Energy_states_L_{L}_T_{T :.3f}000_{file_runs}_random.txt"
            cur = np.loadtxt(filename, usecols=range(1), skiprows=1)
            values.append((T, np.sum(cur) / N / runs))
        np.savetxt(filename_save, values)


def compute_magnetism(sizes, temps, file_runs, runs):
    for L in sizes:
        values = []
        N = L * L
        filename_save = f"data/Computed_Magnetism_states_L_{L}.txt"
        for T in temps:
            filename = (
                f"data/Magnetism_states_L_{L}_T_{T :.3f}000_{file_runs}_random.txt"
            )
            cur = np.loadtxt(filename, usecols=range(1), skiprows=1)
            values.append((T, np.sum(np.abs(cur)) / N / runs))
        np.savetxt(filename_save, values)


def compute_HC(sizes, temps, file_runs, runs):
    for L in sizes:
        values = []
        N = L * L
        filename_save = f"data/Computed_HC_states_L_{L}.txt"
        for T in temps:
            filename = f"data/Energy_states_L_{L}_T_{T :.3f}000_{file_runs}_random.txt"
            cur = np.loadtxt(filename, usecols=range(1), skiprows=1)
            E_squared = (np.sum(cur) / runs) ** 2
            E_exp_squared = np.sum(np.power(cur, 2)) / runs
            values.append((T, (E_exp_squared - E_squared) / N / T / T))
        np.savetxt(filename_save, values)


def compute_susc(sizes, temps, file_runs, runs):
    for L in sizes:
        values = []
        N = L * L
        filename_save = f"data/Computed_susc_states_L_{L}.txt"
        for T in temps:
            filename = (
                f"data/Magnetism_states_L_{L}_T_{T :.3f}000_{file_runs}_random.txt"
            )
            cur = np.loadtxt(filename, usecols=range(1), skiprows=1)
            M_squared = (np.sum(np.abs(cur)) / runs) ** 2
            M_exp_squared = np.sum(np.power(cur, 2)) / runs
            values.append((T, (M_exp_squared - M_squared) / N / T / T))
        np.savetxt(filename_save, values)


def plot_by_temp():
    """Script to plot quantities as a function of temperature from files"""

    sizes = [40, 60, 80, 100]
    temps = np.arange(2.1, 2.401, 0.006)
    runs = 2000000
    file_runs = 200000

    compute_energy(sizes, temps, file_runs, runs)
    compute_magnetism(sizes, temps, file_runs, runs)
    compute_HC(sizes, temps, file_runs, runs)
    compute_susc(sizes, temps, file_runs, runs)

    energy = []
    magnet = []
    specific = []
    susc = []

    for L in sizes:
        cur_E = np.loadtxt(f"data/Computed_Energy_states_L_{L}.txt")
        cur_M = np.loadtxt(f"data/Computed_Magnetism_states_L_{L}.txt")
        cur_specific = np.loadtxt(f"data/Computed_HC_states_L_{L}.txt")
        cur_susc = np.loadtxt(f"data/Computed_susc_states_L_{L}.txt")

        energy.append(cur_E)
        magnet.append(cur_M)
        specific.append(cur_specific)
        susc.append(cur_susc)

    for i in range(len(energy)):
        plt.plot(energy[i][:, 0], energy[i][:, 1], label=f"L = {sizes[i]}")

    plt.legend()
    plt.grid()
    plt.title("Approximation of <$\epsilon$> for different grid lengths L")
    plt.ylabel("Expected energy per spin [J]")
    plt.xlabel("Temperature T [J / $k_b$]")
    plt.savefig("plots/Energy_by_tmp.pgf")
    plt.clf()

    for i in range(len(magnet)):
        plt.plot(magnet[i][:, 0], magnet[i][:, 1], label=f"L = {sizes[i]}")
    plt.legend()
    plt.grid()
    plt.title("Approximation of <$|m|$> for different grid lengths L")
    plt.ylabel("Expected absolute value of magnetisation per spin [1]")
    plt.xlabel("Temperature T [J / $k_b$]")
    plt.savefig("plots/Magnet_by_tmp.pgf")
    plt.clf()

    for i in range(len(specific)):
        plt.plot(specific[i][:, 0], specific[i][:, 1], label=f"L = {sizes[i]}")
    plt.legend()
    plt.title("Approximation of $C_V$ for different grid lengths L")
    plt.ylabel(r"Specific heat capacity [$k_B$]")
    plt.xlabel("Temperature T [J / $k_b$]")
    plt.grid()
    plt.savefig("plots/HC_by_tmp.pgf")
    plt.clf()

    for i in range(len(susc)):
        plt.plot(susc[i][:, 0], susc[i][:, 1], label=f"L = {sizes[i]}")
    plt.legend()
    plt.title("Approximation of $\chi$ for different grid lengths L")
    plt.ylabel("Susceptibility [1/J]")
    plt.xlabel("Temperature T [J / $k_b$]")
    plt.grid()
    plt.savefig("plots/Susc_by_tmp.pgf")
    plt.clf()


def estimate_T_inf():
    """Script to compute T_c and estimate T_c for L = infinity using Specific Heat capacity"""

    temps = np.arange(2.1, 2.401, 0.006)
    sizes = np.array([40, 60, 80, 100])
    values = 1 / sizes
    critical_temps = np.zeros(len(sizes))

    for i, L in enumerate(sizes):
        cur = np.loadtxt(f"data/Computed_HC_states_L_{L}.txt")
        critical_temps[i] = cur[np.argmax(cur[:, 1])][0]
    estimate = sc.stats.linregress(values, critical_temps)

    for L, T_c in zip(sizes, critical_temps):
        print(f"The critical temperature for L = {L} is {T_c : .3f}")

    x = np.linspace(0, np.max(values), 10000)
    print("Estimate of T_c for infinite grid:", estimate.intercept)
    plt.plot(
        x,
        estimate.intercept + estimate.slope * x,
        label=f"Linear regression using $T_c$",
        color="r",
    )

    plt.plot(values, critical_temps, "ro", label=f"Observed $T_c$")
    plt.plot(
        0,
        estimate.intercept,
        "bo",
        label=f"Approximation of $T_c(\infty) = 2.2689$",
        markersize=10,
        color="midnightblue",
    )

    plt.yticks(np.arange(2.269, 2.29, 0.004))
    plt.grid()
    plt.legend()
    plt.xlabel(r"$L^{-1}$ [1]")
    plt.ylabel(r"$T_c [J / k_B]$")
    plt.title(
        f"Estimate of critical temperature $T_c$ for $\infty$x$\infty$-grid using $C_V$"
    )


if __name__ == "__main__":
    newpath = r"plots"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    plt.rc("pgf", texsystem="pdflatex")

    # Problem 4
    plot_energy_per_spin(
        "data/Energy_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    )
    # plt.show()
    plt.savefig("plots/Energy_L_2_T_1.pgf")
    plt.clf()

    plot_magnetisation_per_spin(
        "data/Magnetism_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    )
    plt.savefig("plots/Magnet_L_2_T_1.pgf")
    plt.clf()

    plot_specific_heat_capacity(
        "data/Energy_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    )
    plt.savefig("plots/HC_L_2_T_1.pgf")
    plt.clf()

    plot_susceptibility(
        "data/Magnetism_states_L_2_T_1.000000_1000000_random.txt", 2, 7, 1, 1
    )
    plt.savefig("plots/susceptibility_L_2_T_1.pgf")
    plt.clf()

    # Problem 5
    plot_energy_per_spin(
        "data/Energy_states_L_20_T_1.000000_1000000_random.txt", 20, 7, 1, 1
    )

    plot_energy_per_spin(
        "data/Energy_states_L_20_T_1.000000_1000000_positiv.txt", 20, 7, 0, 1
    )
    plt.grid()
    plt.savefig("plots/Energy_L_20_T_1.pgf")
    plt.clf()

    plot_energy_per_spin(
        "data/Energy_states_L_20_T_2.400000_1000000_random.txt", 20, 7, 1, 2.4
    )

    plot_energy_per_spin(
        "data/Energy_states_L_20_T_2.400000_1000000_positiv.txt", 20, 7, 0, 2.4
    )
    plt.grid()
    plt.savefig("plots/Energy_L_20_T_2.4.pgf")
    plt.clf()

    plot_magnetisation_per_spin(
        "data/Magnetism_states_L_20_T_1.000000_1000000_random.txt", 20, 7, 1, 1
    )
    plot_magnetisation_per_spin(
        "data/Magnetism_states_L_20_T_1.000000_1000000_positiv.txt", 20, 7, 0, 1
    )
    plt.grid()
    plt.savefig("plots/Magnet_L_20_T_1.pgf")
    plt.clf()

    plot_magnetisation_per_spin(
        "data/Magnetism_states_L_20_T_2.400000_1000000_random.txt", 20, 7, 1, 2.4
    )
    plot_magnetisation_per_spin(
        "data/Magnetism_states_L_20_T_2.400000_1000000_positiv.txt", 20, 7, 0, 2.4
    )
    plt.grid()
    plt.savefig("plots/Magnet_L_20_T_2.4.pgf")
    plt.clf()

    # Problem 6
    generate_historgram(
        "data/Energy_states_L_20_T_1.000000_1000000_random.txt",
        20,
        1,
        np.arange(-2, -1.9, 0.005),
        1,
    )
    plt.savefig("plots/Probs_L_20_T_1.pgf")
    plt.clf()

    generate_historgram(
        "data/Energy_states_L_20_T_2.400000_1000000_random.txt",
        20,
        2.4,
        np.arange(-1.9, -0.6, 0.005),
        0,
    )
    plt.savefig("plots/Probs_L_20_T_2.4.pgf")
    plt.clf()

    # Problem 8
    plot_by_temp()
    plt.clf()
    estimate_T_inf()
    plt.savefig("plots/Critical_temperature_L_inf.pgf")
