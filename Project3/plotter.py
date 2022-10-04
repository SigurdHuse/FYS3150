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
    plt.plot(time, z)


def plot_two_particles(interaction):
    if interaction:
        filename = "two_particles_with_interaction.txt"
    else:
        filename = "two_particles_without_interaction.txt"

    main = np.loadtxt(filename, delimiter=",", usecols=range(2))
    x, y = main[::3], main[1::3]
    p1_x, p1_y = x[:, 0], y[:, 0]
    p2_x, p2_y = x[:, 1], y[:, 1]

    plt.subplot(2, 1, 1)
    plt.plot(p1_x, p1_y, label="first particle", color="r")
    plt.plot(p1_x[0], p1_y[0], marker="o", color="r")

    plt.subplot(2, 1, 2)
    plt.plot(p2_x, p2_y, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], p2_y[0], marker="o", color="midnightblue")


def plot_phase_two_particles(interaction):
    if interaction:
        filename1 = "two_particles_with_interaction.txt"
        filename2 = "two_particles_with_interaction_vel.txt"
    else:
        filename1 = "two_particles_without_interaction.txt"
        filename2 = "two_particles_without_interaction_vel.txt"

    pos = np.loadtxt(filename1, delimiter=",", usecols=range(2))
    vel = np.loadtxt(filename2, delimiter=",", usecols=range(2))

    x, y = pos[::3], pos[1::3]
    p1_x, p1_y = x[:, 0], y[:, 0]
    p2_x, p2_y = x[:, 1], y[:, 1]

    vx, vy = vel[::3], vel[1::3]
    v1_x, v1_y = vx[:, 0], vy[:, 0]
    v2_x, v2_y = vx[:, 1], vy[:, 1]

    plt.subplot(2, 1, 1)
    plt.plot(p1_x, v1_x, label="first particle", color="r")
    plt.plot(p1_x[0], v1_x[0], marker="o", color="midnightblue")

    plt.subplot(2, 1, 2)
    plt.plot(p2_x, v2_x, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], v2_x[0], marker="o", color="r")


def generate_exact_solution(x0, v0, m, q, d, B0, V0, T):
    omega0 = q * B0 / m
    omegaz = 2 * q * V0 / m / d / d
    omega_minus = (omega0 - np.sqrt(omega0**2 - 2 * omegaz)) / 2
    omega_plus = (omega0 + np.sqrt(omega0**2 - 2 * omegaz)) / 2
    A_plus = (v0 + omega_minus * x0) / (omega_minus - omega_plus)
    A_minus = -(v0 + omega_plus * x0) / (omega_minus - omega_plus)
    phi_plus = 0
    phi_minus = 0
    t = np.linspace(0, T, T * 1000)
    values = A_plus * np.exp(-1j * (omega_plus * t + phi_plus)) + A_minus * np.exp(
        -1j * (omega_minus * t + phi_minus)
    )
    return values, t, omegaz


def plot_exact_solution(x0, z0, v0, m, q, d, B0, V0, T, plot_z):
    values, t, omegaz = generate_exact_solution(x0, v0, m, q, d, B0, V0, T)
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
    plt.show()
    # plt.xlim(xlim)
    # plt.ylim(ylim)
    # plt.rc("pgf", texsystem="pdflatex")
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


if __name__ == "__main__":
    """
    plot_two_particles(1)
    save_plot("with interaction")
    plot_two_particles(0)
    save_plot("without interaction")
    plot_phase_two_particles(1)
    save_plot("with interaction")
    plot_phase_two_particles(0)
    save_plot("without interaction")
    plot_one_particle()
    save_plot("Numerical approximation of one particle")
    plot_3D_two_particles(False)
    save_plot("3D plot")
    """
    plot_exact_solution(10, 10, 5, 40.0775, 1, 10**4, 96.5, 9.65 * 10**8, 100, True)
    save_plot("Exact solution for one particle")
    # plot_one_particle()
    # save_plot("Numerical approximation of one particle")
