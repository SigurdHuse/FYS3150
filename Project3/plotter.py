import numpy as np
import matplotlib.pyplot as plt


def plot_one_particle():
    main = np.loadtxt("one_particle_n_10000.txt", delimiter=",", usecols=range(1))
    time = np.loadtxt(
        "time_one_particle_n_10000.txt",
        delimiter=",",
    )
    x, y, z = main[::3], main[1::3], main[2::3]
    plt.plot(time, z)
    plt.show()


def plot_two_particles(interaction, title):
    if interaction:
        filename = "two_particles_with_interaction.txt"
    else:
        filename = "two_particles_without_interaction.txt"

    main = np.loadtxt(filename, delimiter=",", usecols=range(2))
    x, y = main[::3], main[1::3]
    p1_x, p1_y = x[:, 0], y[:, 0]
    p2_x, p2_y = x[:, 1], y[:, 1]

    plt.plot(p1_x, p1_y, label="first particle", color="r")
    plt.plot(p1_x[0], p1_y[0], marker="o", color="r")
    plt.plot(p2_x, p2_y, label="second particle", color="midnightblue")
    plt.plot(p2_x[0], p2_y[0], marker="o", color="midnightblue")
    plt.title(title)
    plt.show()


def plot_phase_two_particles(interaction, title):
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

    plt.plot(p1_x, v1_x, label="first particle", color="r")
    plt.plot(p2_x, v2_x, label="second particle", color="midnightblue")
    plt.show()


if __name__ == "__main__":
    # plot_two_particles(1, "with interaction")
    plt.clf()
    # plot_two_particles(0, "without interaction")
    plt.clf()
    plot_phase_two_particles(1, "with interaction")
    plt.clf()
    plot_phase_two_particles(0, "without interaction")
    plt.clf()
    # plot_one_particle()
