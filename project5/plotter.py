import pyarma as pa
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import numpy as np


def plot_probability(M, steps, slits):
    A = pa.cx_cube()

    values = np.loadtxt("probs.txt")
    times = list(range(steps + 1))

    """
    A.load(f"Data_M_{M}_dt_{steps}_slits_{slits}.bin")

    for idx in times:
        matrix = np.zeros((M, M), dtype="complex")
        for i in range(M):
            for j in range(M):
                matrix[i, j] = A[i, j, idx]

        matrix = np.power(np.abs(matrix), 2)
        values[idx] = np.sum(matrix)

    np.savetxt("probs.txt", values)
    """
    # print(1 - values)
    plt.plot(times, np.abs(1.0 - values))
    plt.yscale("log")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    plot_probability(200, 320, 0)
