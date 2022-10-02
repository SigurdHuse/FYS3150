import numpy as np
import matplotlib.pyplot as plt


def plot_one_particle():
    main = np.loadtxt("one_particle_n_10000.txt", delimiter=",", usecols=range(1))
    t, z = main[::4], main[3::4]
    plt.plot(t, z)
    plt.show()


def plot_two_particles(interaction, title):
    if interaction:
        filename = "two_particles_with_interaction.txt"
    else:
        filename = "two_particles_without_interaction.txt"

    main = np.loadtxt(filename, delimiter=",", usecols=range(2))
    x, y = main[::4], main[1::4]
    p1_x, p1_y = x[0, :], y[0, :]
    p2_x, p2_y = x[1, :], y[1, :]

    plt.plot(p1_x, p1_y, label="first particle", color="r")
    plt.plot(p2_x, p2_y, label="second particle", color="midnightblue")
    plt.title(title)
    plt.show()


if __name__ == "__main__":
    plot_two_particles(1, "with interaction")
    plt.clf()
    plot_two_particles(0, "without interaction")
