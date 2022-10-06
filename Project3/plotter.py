import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


def plot_one_particle():
    main = np.loadtxt("one_particle_n_10000.txt", delimiter=",", usecols=range(1))
    time = np.loadtxt(
        "time_one_particle_n_10000.txt",
        delimiter=",",
    )
    x, y, z = main[::3], main[1::3], main[2::3]
    plt.plot(x, y, label="One particle")
    plt.xlabel(f"time ($\mu$s)")
    plt.ylabel("z-position ($\mu$m)")


def plot_two_particles(interaction):
    if interaction:
        filename = "two_particles_with_interaction.txt"
        title = "Plot of two particles with interaction"
    else:
        filename = "two_particles_without_interaction.txt"
        title = "Plot of two particles without interaction"

    main = np.loadtxt(filename, delimiter=",", usecols=range(2))
    x, y = main[::3], main[1::3]
    p1_x, p1_y = x[:, 0], y[:, 0]
    p2_x, p2_y = x[:, 1], y[:, 1]

    plt.subplot(2, 1, 1)
    plt.plot(p1_x, p1_y, label="first particle", color="r")
    plt.plot(p1_x[0], p1_y[0], marker="o", color="midnightblue")
    plt.xlim([-25, 25])
    plt.ylim([-25, 25])
    plt.title(title)
    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(p2_x, p2_y, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], p2_y[0], marker="o", color="r")
    plt.xlim([-25, 25])
    plt.ylim([-25, 25])
    plt.legend()
    plt.grid()


def plot_phase_two_particles(interaction):
    if interaction:
        filename1 = "two_particles_with_interaction.txt"
        filename2 = "two_particles_with_interaction_vel.txt"
    else:
        filename1 = "two_particles_without_interaction.txt"
        filename2 = "two_particles_without_interaction_vel.txt"

    pos = np.loadtxt(filename1, delimiter=",", usecols=range(2))
    vel = np.loadtxt(filename2, delimiter=",", usecols=range(2))

    filename = "phase_space"
    if interaction:
        filename += "_interaction"
    filename += "_x.pgf"

    x, y, z = pos[:-3:3], pos[1:-2:3], pos[2:-1:3]
    p1_x, p1_y, p1_z = x[:, 0], y[:, 0], z[:, 0]
    p2_x, p2_y, p2_z = x[:, 1], y[:, 1], z[:, 1]

    vx, vy, vz = vel[:-3:3], vel[1:-2:3], vel[2:-1:3]
    v1_x, v1_y, v1_z = vx[:, 0], vy[:, 0], vz[:, 0]
    v2_x, v2_y, v2_z = vx[:, 1], vy[:, 1], vz[:, 1]

    plt.subplot(2, 1, 1)
    plt.plot(p1_x, v1_x, label="first particle", color="r")
    plt.plot(p1_x[0], v1_x[0], marker="o", color="midnightblue")
    plt.title("Phase space plot in x-plane" + " with interaction" * (interaction == 1))
    plt.xlim([-16, 16])
    plt.ylim([-16, 16])

    plt.subplot(2, 1, 2)
    plt.plot(p2_x, v2_x, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], v2_x[0], marker="o", color="r")
    plt.xlim([-16, 16])
    plt.ylim([-16, 16])

    plt.savefig(filename)
    plt.clf()

    plt.subplot(2, 1, 1)
    plt.plot(p1_y, v1_y, label="first particle", color="r")
    plt.plot(p1_y[0], v1_y[0], marker="o", color="midnightblue")
    plt.title("Phase space plot in y-plane" + " with interaction" * (interaction == 1))
    plt.xlim([-16, 16])
    plt.ylim([-16, 16])

    plt.subplot(2, 1, 2)
    plt.plot(p2_y, v2_y, label="second particle", color="midnightblue")
    plt.plot(p2_y[0], v2_y[0], marker="o", color="r")

    filename = filename[:-5] + "y" + filename[-4:]

    plt.savefig(filename)
    plt.clf()

    plt.subplot(2, 1, 1)
    plt.plot(p1_z, v1_z, label="first particle", color="r")
    plt.plot(p1_z[0], v1_z[0], marker="o", color="midnightblue")
    plt.title("Phase space plot in z-plane" + " with interaction" * (interaction == 1))
    plt.xlim([-16, 16])
    plt.ylim([-16, 16])

    plt.subplot(2, 1, 2)
    plt.plot(p2_z, v2_z, label="second particle", color="midnightblue")
    plt.plot(p2_z[0], v2_z[0], marker="o", color="r")

    filename = filename[:-5] + "z" + filename[-4:]

    plt.savefig(filename)
    plt.clf()


def generate_exact_solution(x0, v0, m, q, d, B0, V0, T, n):
    omega0 = q * B0 / m
    omegaz = 2 * q * V0 / m / d / d
    omega_minus = (omega0 - np.sqrt(omega0**2 - 2 * omegaz)) / 2
    omega_plus = (omega0 + np.sqrt(omega0**2 - 2 * omegaz)) / 2
    A_plus = (v0 + omega_minus * x0) / (omega_minus - omega_plus)
    A_minus = -(v0 + omega_plus * x0) / (omega_minus - omega_plus)
    phi_plus = 0
    phi_minus = 0
    t = np.linspace(0, T, n + 1)
    values = A_plus * np.exp(-1j * (omega_plus * t + phi_plus)) + A_minus * np.exp(
        -1j * (omega_minus * t + phi_minus)
    )
    return values, t, omegaz


def plot_exact_solution(x0, z0, v0, m, q, d, B0, V0, T, plot_z):
    values, t, omegaz = generate_exact_solution(x0, v0, m, q, d, B0, V0, T, T * 10000)
    z = z0 * np.cos(np.sqrt(omegaz) * t)
    x, y = np.real(values), np.imag(values)
    if plot_z:
        plt.plot(t, z, label="Exact solution for one particle")
    else:
        plt.plot(x, y, label="Exact solution for one particle")


def save_plot(title, xlim=[], ylim=[], filename=""):
    plt.grid()
    plt.title(title)
    plt.legend()

    plt.xlim(xlim)
    plt.ylim(ylim)
    # plt.rc("pgf", texsystem="pdflatex")
    plt.show()
    # plt.savefig(filename)
    plt.clf()


def plot_3D_two_particles(interaction):
    if interaction:
        filename = "two_particles_with_interaction.txt"
    else:
        filename = "two_particles_without_interaction.txt"

    main = np.loadtxt(filename, delimiter=",", usecols=range(2))
    x, y, z = main[::3], main[1::3], main[2::3]
    p1_x, p1_y, p1_z = x[:, 0], y[:, 0], z[:, 0]
    p2_x, p2_y, p2_z = x[:, 1], y[:, 1], z[:, 1]

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.plot3D(p1_x, p1_y, p1_z, label="Particle one", color="midnightblue")
    ax.plot3D(p1_x[0], p1_y[0], p1_z[0], color="r", marker="o")

    ax.plot3D(p2_x, p2_y, p2_z, label="Particle two", color="r")
    ax.plot3D(p2_x[0], p2_y[0], p2_z[0], color="midnightblue", marker="o")


def plot_relative_error(RK4_method):
    nvals = [10**i for i in range(1, 6)]
    fig, axs = plt.subplots(3, 2)
    fig.tight_layout(pad=3.0)
    x_plot, y_plot = 0, 0
    for n in nvals:
        RK4 = np.loadtxt(
            "one_particle_n_" + str(n) + "_method_RK4.txt",
            delimiter=",",
            usecols=range(1),
        )
        euler = np.loadtxt(
            "one_particle_n_" + str(n) + "_method_euler.txt",
            delimiter=",",
            usecols=range(1),
        )
        exact, t, omegaz = generate_exact_solution(
            10, 5, 40.0775, 1, 10**4, 96.5, 9.65 * 10**8, 10, n
        )

        x_RK4, y_RK4, z_RK4 = RK4[3::3], RK4[4::3], RK4[5::3]
        x_euler, y_euler, z_euler = euler[3::3], euler[4::3], euler[5::3]
        x, y, z = (
            np.real(exact[1:]),
            np.imag(exact[1:]),
            10 * np.cos(np.sqrt(omegaz) * t[1:]),
        )

        error_RK4 = (
            np.abs((x_RK4 - x) / x) + np.abs((y_RK4 - y) / y) + np.abs((z_RK4 - z) / z)
        )
        error_euler = (
            np.abs((x_euler - x) / x)
            + np.abs((y_euler - y) / y)
            + np.abs((z_euler - z) / z)
        )
        t = t[1:]
        if RK4_method:
            axs[y_plot, x_plot].plot(t, error_RK4, label="RK4")

        else:
            axs[y_plot, x_plot].plot(t, error_euler, label="Euler")

        axs[y_plot, x_plot].set_title(f"Relative error when h = {10/n}")
        axs[y_plot, x_plot].grid()
        axs[y_plot, x_plot].set_xlim([0, 10])
        axs[y_plot, x_plot].set_yscale("log")
        x_plot += 1
        if x_plot == 2:
            x_plot = 0
            y_plot += 1

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc="lower right")


def compute_error_convergence_rate(RK4):
    nvals = [10**i for i in range(1, 6)]
    delta_values = []
    h_values = []
    for n in nvals:
        if RK4:
            data = np.loadtxt(
                "one_particle_n_" + str(n) + "_method_RK4.txt",
                delimiter=",",
                usecols=range(1),
            )
        else:
            data = np.loadtxt(
                "one_particle_n_" + str(n) + "_method_euler.txt",
                delimiter=",",
                usecols=range(1),
            )
        exact, t, omegaz = generate_exact_solution(
            10, 5, 40.0775, 1, 10**4, 96.5, 9.65 * 10**8, 10, n
        )
        x_data, y_data, z_data = data[3::3], data[4::3], data[5::3]
        x, y, z = (
            np.real(exact[1:]),
            np.imag(exact[1:]),
            10 * np.cos(np.sqrt(omegaz) * t[1:]),
        )
        r_max = np.max(np.abs(x_data - x) + np.abs(y_data - y) + np.abs(z_data - z))
        delta_values.append(r_max)
        h_values.append(10 / n)
    r_err = 0.0
    for k in range(1, 5):
        r_err += np.log(delta_values[k] / delta_values[k - 1]) / np.log(
            h_values[k] / h_values[k - 1]
        )
    r_err /= 4
    if RK4:
        print("Error convergence rate for RK4 is ", r_err)
    else:
        print("Error convergence rate for forward Euler is ", r_err)


def plot_number_of_particles_trapped():
    f_vals = [0.1, 0.4, 0.7]
    number = 1
    for f in f_vals:
        main = np.loadtxt(f"hundred_particles_f_{f}.txt", usecols=range(1))
        omega, particles = main[::2], main[1::2]
        plt.subplot(3, 1, number)
        plt.plot(omega, particles, label=f"f = {f}")
        plt.legend()
        plt.grid()
        if number == 2:
            plt.ylabel(f"Number of particles left after 500 $\mu$s")
        plt.xlabel(f"$\omega_V$")
        number += 1
    plt.show()


if __name__ == "__main__":
    plt.rc("pgf", texsystem="pdflatex")
    """
    plot_one_particle()
    save_plot("Numerical approximation of one particle")
    plot_3D_two_particles(False)
    save_plot("3D plot")
    """
    # plot_one_particle()
    # plt.show()
    # save_plot(
    #    f"Numerical approximation of one particle for 100 $\mu$s",
    #    [0, 100],
    #    [-11, 11],
    #    "one_particle_time_100.pgf",
    # )

    # plot_two_particles(1)
    # plt.savefig("two_particles_with_interaction.pgf")
    # plt.clf()
    # plot_two_particles(0)
    # plt.savefig("two_particles_without_interaction.pgf")
    # plt.clf()

    # plot_phase_two_particles(1)
    # plot_phase_two_particles(0)

    # plot_exact_solution(10, 10, 5, 40.0775, 1, 10**4, 96.5, 9.65 * 10**8, 100, True)
    # save_plot("Exact solution for one particle")
    #
    # plot_relative_error(True)
    # plt.show()
    # compute_error_convergence_rate(True)
    # compute_error_convergence_rate(False)
    plot_number_of_particles_trapped()
