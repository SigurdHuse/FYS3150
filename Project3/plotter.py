import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["figure.titlesize"] = 16
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["legend.fontsize"] = "large"
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 16


def plot_one_particle():
    main = np.loadtxt("one_particle_n_10000.txt", delimiter=",", usecols=range(1))
    time = np.loadtxt(
        "time_one_particle_n_10000.txt",
        delimiter=",",
    )
    x, y, z = main[::3], main[1::3], main[2::3]
    plt.xlim([0, 50])
    plt.xticks([i * 5 for i in range(11)])
    plt.ylim([-21, 21])
    plt.grid()
    plt.plot(time, z, label="One particle simulated with RK4", color="midnightblue")
    plt.title("Plot of one particle in a Penning trap simulated with RK4")
    plt.legend(loc="upper center")
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

    plt.plot(p1_x, p1_y, label="first particle", color="r")
    plt.plot(p1_x[0], p1_y[0], marker="o", color="midnightblue")
    plt.plot(p2_x, p2_y, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], p2_y[0], marker="o", color="r")
    plt.axis("equal")
    plt.grid()
    plt.legend()
    plt.xlabel("x-position ($\mu$m)")
    plt.ylabel("y-position ($\mu$m)")
    plt.title(title)
    """
    fig, ax = plt.subplots(2, 1)
    fig.tight_layout(pad=3.0)

    ax[0].plot(p1_x, p1_y, label="first particle", color="r")
    ax[0].plot(p1_x[0], p1_y[0], marker="o", color="midnightblue")
    ax[0].set_aspect("equal", adjustable="box", anchor="C")
    ax[0].set_xlabel("x-position ($\mu$m)")
    ax[0].set_ylabel("y-position ($\mu$m)")
    ax[0].set_title(title)
    ax[0].legend()
    ax[0].grid()

    ax[1].plot(p2_x, p2_y, label="second particle", color="midnightblue")
    ax[1].plot(p2_x[0], p2_y[0], marker="o", color="r")
    ax[1].set_xlabel("x-position ($\mu$m)")
    ax[1].set_ylabel("y-position ($\mu$m)")
    ax[1].set_aspect("equal", adjustable="box", anchor="C")
    ax[1].legend()
    ax[1].grid()
    """


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

    # Plot x-plane
    plt.subplot(2, 1, 1)
    plt.plot(p1_x, v1_x, label="first particle", color="r")
    plt.plot(p1_x[0], v1_x[0], marker="o", color="midnightblue")
    plt.title(
        "Phase space plot in x-plane"
        + " with"
        + "out" * (interaction == 0)
        + " interaction"
    )
    plt.ylabel(f"v-x ($\mu m / \mu s)$")
    plt.xlim([-70, 50])
    plt.ylim([-50, 50])
    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(p2_x, v2_x, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], v2_x[0], marker="o", color="r")
    plt.ylabel(f"v-x ($\mu m / \mu s$)")
    plt.xlabel(f"x-position ($\mu m$)")
    plt.xlim([-70, 80])
    plt.ylim([-60, 60])
    plt.legend()
    plt.grid()

    plt.savefig(filename)
    plt.clf()

    # Plot y-plane
    plt.subplot(2, 1, 1)
    plt.plot(p1_y, v1_y, label="first particle", color="r")
    plt.plot(p1_y[0], v1_y[0], marker="o", color="midnightblue")
    plt.title(
        "Phase space plot in y-plane"
        + " with"
        + "out" * (interaction == 0)
        + " interaction"
    )
    plt.xlim([-60, 15])
    plt.ylim([-40, 40])
    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(p2_y, v2_y, label="second particle", color="midnightblue")
    plt.plot(p2_y[0], v2_y[0], marker="o", color="r")
    plt.xlim([-80, 60])
    plt.ylim([-60, 60])
    plt.legend()
    plt.grid()

    filename = filename[:-5] + "y" + filename[-4:]

    plt.savefig(filename)
    plt.clf()

    # Plot z - plane
    plt.subplot(2, 1, 1)
    plt.plot(p1_z, v1_z, label="first particle", color="r")
    plt.plot(p1_z[0], v1_z[0], marker="o", color="midnightblue")
    plt.title(
        "Phase space plot in z-plane"
        + " with"
        + "out" * (interaction == 0)
        + " interaction"
    )
    plt.xlim([-50, 50])
    plt.ylim([-40, 40])
    if interaction:
        plt.xlim([-30, 30])
        plt.ylim([-30, 30])
    plt.ylabel(f"v-z ($\mu m / \mu s$)")
    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(p2_z, v2_z, label="second particle", color="midnightblue")
    plt.plot(p2_z[0], v2_z[0], marker="o", color="r")
    plt.xlim([-11, 11])
    plt.ylim([-11, 11])
    if interaction:
        plt.xlim([-25, 25])
        plt.ylim([-20, 20])
    plt.ylabel(f"v-z ($\mu m / \mu s$)")
    plt.xlabel(f"z-position ($\mu m$)")
    plt.legend()
    plt.grid()

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
        plt.plot(t, z, label="Analytic solution for one particle", color="midnightblue")
    else:
        plt.plot(x, y, label="Analytic solution for one particle", color="midnightblue")
    plt.grid()
    plt.legend(loc="upper center")
    plt.title("Plot of analytic solution for one particle in Penning trap")
    plt.xlabel(f"time ($\mu$s)")
    plt.ylabel("z-position ($\mu$m)")


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
    ax.set_title(
        "3d plot of two particles with" + "out" * (interaction == 0) + " interaction"
    )
    ax.set_zlabel("z-position ($\mu$m)")
    ax.set_ylabel("y-position ($\mu$m)")
    ax.set_xlabel("x-position ($\mu$m)")
    ax.legend()


def plot_relative_error(RK4_method):
    nvals = [4000, 8000, 16000, 32000]
    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=4.0)
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
            20, 25, 40.0775, 1, 500, 9.65 * 1e1, 2.41 * 1e6, 50, n
        )

        x_RK4, y_RK4, z_RK4 = RK4[3::3], RK4[4::3], RK4[5::3]
        x_euler, y_euler, z_euler = euler[3::3], euler[4::3], euler[5::3]
        x, y, z = (
            np.real(exact[1:]),
            np.imag(exact[1:]),
            20 * np.cos(np.sqrt(omegaz) * t[1:]),
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
            if n == 4000:
                axs[y_plot, x_plot].plot(t, error_RK4, label="RK4")
            else:
                axs[y_plot, x_plot].plot(t, error_RK4)
        else:
            if n == 4000:
                axs[y_plot, x_plot].plot(t, error_euler, label="Euler")
            else:
                axs[y_plot, x_plot].plot(t, error_euler)

        axs[y_plot, x_plot].set_title(f"Relative error when n = {n}")
        axs[y_plot, x_plot].grid()
        axs[y_plot, x_plot].set_xlim([0, 50])
        axs[y_plot, x_plot].set_yscale("log")
        axs[y_plot, x_plot].set_xlabel(f"time ($\mu s$)")
        axs[y_plot, x_plot].set_ylabel(f"Relative error")
        x_plot += 1
        if x_plot == 2:
            x_plot = 0
            y_plot += 1

    # fig.suptitle("Plot of relative error for different values of h")
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc="lower center")


def compute_error_convergence_rate(RK4):
    nvals = [4000, 8000, 16000, 32000]
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
            20, 25, 40.0775, 1, 500, 9.65 * 1e1, 2.41 * 1e6, 50, n
        )
        x_data, y_data, z_data = data[3::3], data[4::3], data[5::3]
        x, y, z = (
            np.real(exact[1:]),
            np.imag(exact[1:]),
            20 * np.cos(np.sqrt(omegaz) * t[1:]),
        )

        r_max = np.max(np.abs(x - x_data) + np.abs(y - y_data) + np.abs(z - z_data))
        delta_values.append(r_max)
        h_values.append(50 / n)
    r_err = 0.0
    for k in range(1, 4):
        r_err += np.log(delta_values[k] / delta_values[k - 1]) / np.log(
            h_values[k] / h_values[k - 1]
        )
    r_err /= 3
    if RK4:
        print("Error convergence rate for RK4 is ", r_err)
    else:
        print("Error convergence rate for forward Euler is ", r_err)


def plot_number_of_particles_trapped():
    f_vals = [0.1, 0.4, 0.7]
    number = 1
    nr_of_particles = 100
    for f in f_vals:
        main = np.loadtxt(f"hundred_particles_f_{f}.txt", usecols=range(1))
        omega, particles = main[::2], main[1::2]
        plt.subplot(3, 1, number)
        plt.plot(omega, particles / nr_of_particles, "bo", label=f"f = {f}")
        plt.legend()
        plt.grid()
        if number == 2:
            plt.ylabel(f"Fraction of trapped particles after 500 $\mu$s")
        plt.xlabel(f"$\omega_V$ (MHz)")
        number += 1
        plt.xticks(np.arange(0.2, 2.5, 0.2))


def plot_number_of_particles_trapped_fine():
    without = np.loadtxt("hundred_particles_f_0.4_fine.txt", usecols=range(1))
    interaction = np.loadtxt(
        "hundred_particles_f_0.4_fine_interaction.txt", usecols=range(1)
    )
    nr_of_particles = 100
    omega1, omega2 = without[::2], interaction[::2]
    particle1, particle2 = without[1::2], interaction[1::2]

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(omega1, particle1 / nr_of_particles, "bo", label=f"without interaction")
    ax[0].legend()
    ax[0].set_xlabel(f"$\omega_V$ (MHz)")

    ax[1].plot(omega2, particle2 / nr_of_particles, "bo", label=f"with interaction")
    ax[1].legend()
    ax[1].set_xlabel(f"$\omega_V$ (MHz)")
    fig.supylabel(f"Fraction of trapped particles after 500 $\mu$s")


if __name__ == "__main__":
    plt.rc("pgf", texsystem="pdflatex")
    plot_one_particle()
    plt.savefig("one_particle_time_50.pgf")
    plt.clf()

    plot_3D_two_particles(True)
    plt.savefig("3d_plot_with_interaction.pgf")
    plt.clf()
    plot_3D_two_particles(False)
    plt.savefig("3d_plot_without_interaction.pgf")
    plt.clf()

    plot_two_particles(1)
    plt.savefig("two_particles_with_interaction.pgf")
    plt.clf()
    plot_two_particles(0)
    plt.savefig("two_particles_without_interaction.pgf")
    plt.clf()

    plot_phase_two_particles(1)
    plot_phase_two_particles(0)

    plot_exact_solution(20, 20, 25, 40.0775, 1, 500, 96.5, 2.41 * 1e6, 50, True)
    plt.savefig("Exact solution for one particle.pgf")
    plt.clf()

    plot_relative_error(True)
    plt.savefig("Relative_error_RK4.pgf")
    plt.clf()
    plot_relative_error(False)
    plt.savefig("Relative_error_Euler.pgf")
    plt.clf()

    compute_error_convergence_rate(True)
    compute_error_convergence_rate(False)
    plot_number_of_particles_trapped()
    plt.savefig("Numbers_of_particles_trapped_0.2_2.5_.pgf")
    # plt.clf()
    plot_number_of_particles_trapped_fine()
    plt.savefig("Numbers_of_particles_trapped_2.0_2.2_.pgf")
    # plt.show()
