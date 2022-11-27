import pyarma as pa
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
import numpy as np
import os

mpl.rcParams["figure.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams["legend.fontsize"] = "medium"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12


def plot_probability(name, M, steps, T, slits):
    """Plots sum of probabilities deviation from 1"""

    A = pa.cx_cube()

    values = np.loadtxt(f"probs_{slits}.txt")
    # values = np.zeros(steps + 1)
    times = np.arange(0, steps + 1)

    """
    A.load(f"data/{name}_M_{M}_dt_{steps}_slits_{slits}.bin")

    for idx in times:
        matrix = np.zeros((M, M), dtype="complex")
        for i in range(M):
            for j in range(M):
                matrix[i, j] = A[i, j, idx]

        matrix = np.power(np.abs(matrix), 2)
        values[idx] = np.sum(matrix)

    np.savetxt(f"probs_{slits}.txt", values)
    """

    plt.plot(
        times,
        np.abs(1.0 - values),
        color="midnightblue",
        label="Sum of probabilities deviation from 1",
    )
    plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel("Time step [1]")
    plt.ylabel("Deviation from 1 [1]")
    plt.legend()
    # plt.yticks([10 ** (-16), 10 ** (-15), 10 ** (-14), 10 ** (-13)])
    plt.grid()


def plot_color_map(delta_T, M, T, name, slits):
    """Plots color map at different time steps"""

    A = pa.cx_cube()
    A.load(f"data/{name}_M_{M}_dt_{int(T/delta_T)}_slits_{slits}.bin")

    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=1.0)
    grids = [0] * 4
    times = np.array([0, 0.001, 0.0015, 0.002])

    for idx in range(4):
        matrix = np.zeros((M, M), dtype="complex")
        for i in range(M):
            for j in range(M):
                matrix[i, j] = A[i, j, int(times[idx] / delta_T)]

        grids[idx] = matrix

    for i in range(2):
        for j in range(2):
            axs[i, j].set_xticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_xticklabels([i / 10 for i in range(0, 12, 2)])
            axs[i, j].set_yticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_yticklabels([i / 10 for i in range(10, -1, -2)])
            im = axs[i, j].imshow(np.power(np.abs(grids[2 * i + j]), 2))
            axs[i, j].set_title(f"T = {times[2*i + j]}")
            plt.colorbar(im, ax=axs[i, j])

    plt.savefig("plots/Time_evolution_norm.pdf")

    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=1.0)

    for i in range(2):
        for j in range(2):
            axs[i, j].set_xticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_xticklabels([i / 10 for i in range(0, 12, 2)])
            axs[i, j].set_yticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_yticklabels([i / 10 for i in range(10, -1, -2)])
            im = axs[i, j].imshow(np.real(grids[2 * i + j]))
            axs[i, j].set_title(f"T = {times[2*i + j]}")
            plt.colorbar(im, ax=axs[i, j])

    plt.savefig("plots/Time_evolution_real.pdf")

    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=1.0)

    for i in range(2):
        for j in range(2):
            axs[i, j].set_xticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_xticklabels([i / 10 for i in range(0, 12, 2)])
            axs[i, j].set_yticks(np.linspace(0, M - 1, 6))
            axs[i, j].set_yticklabels([i / 10 for i in range(10, -1, -2)])
            im = axs[i, j].imshow(np.imag(grids[2 * i + j]))
            axs[i, j].set_title(f"T = {times[2*i + j]}")
            plt.colorbar(im, ax=axs[i, j])

    plt.savefig("plots/Time_evolution_imag.pdf")


def plot_screen(delta_T, M, T, name, slits):
    """Plots one line of the grid at x = 0.8 for different time steps with probability normalized to one"""

    A = pa.cx_cube()

    # values = np.loadtxt(f"probs_{slits}.txt")
    # values = np.zeros(steps + 1)
    # times = np.arange(0, steps + 1)

    desired_x_idx = int(0.8 * M)
    time = int(T / delta_T)

    A.load(f"data/{name}_M_{M}_dt_{time}_slits_{slits}.bin")

    values = np.zeros((M, time), dtype="complex")

    for idx in range(time):
        for i in range(M):
            values[i, idx] = A[i, desired_x_idx, idx]

        values[:, idx] /= np.sum(np.abs(values[:, idx]))

    plt.imshow(np.abs(values))
    plt.yticks(np.linspace(0, M - 1, 6), [i / 10 for i in range(0, 12, 2)])
    plt.xticks([0, 40, 80])
    plt.ylabel("Position")
    plt.xlabel("Time step[1]")
    # plt.xticklabels()
    plt.colorbar()
    plt.savefig(f"plots/Time_evolution_slits{slits}.pdf")


if __name__ == "__main__":
    newpath = r"plots"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    plt.rc("pgf", texsystem="pdflatex")

    plot_probability("No_slit_sigma_y_005", 200, 320, 0.008, 0)
    plt.savefig("plots/Probs_slits_0.pgf")
    plt.clf()
    plot_probability("Two_slits_sigma_y_0.1", 200, 320, 0.008, 2)
    plt.savefig("plots/Probs_slits_2.pgf")
    plt.clf()

    plot_color_map(2.5 * 10 ** (-5), 200, 0.002, "Two_slits_sigma_y_0.2", 2)
    plt.clf()
    plot_screen(2.5 * 10 ** (-5), 200, 0.002, "Two_slits_sigma_y_0.2", 2)
    plt.clf()
    plot_screen(2.5 * 10 ** (-5), 200, 0.002, "One_slit_sigma_y_0.2", 1)
    plt.clf()
    plot_screen(2.5 * 10 ** (-5), 200, 0.002, "Three_slits_sigma_y_0.2", 3)
