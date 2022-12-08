import pyarma as pa
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import numpy as np
import os

plt.rcParams["animation.ffmpeg_path"] = "/usr/bin/ffmpeg"


def create_animation(Name, M, time_steps, slits, filename, delta_t):
    """Creates animation of data from simulation"""

    A = pa.cx_cube()

    fig, ax = plt.subplots()

    A.load("data/" + Name + f"_M_{M}_dt_{time_steps}_slits_{slits}.bin")

    length = time_steps + 1
    size = M
    value = np.zeros((size, size), dtype="complex")

    fig = plt.figure()
    ax = plt.gca()

    for i in range(size):
        for j in range(size):
            value[i, j] = A[i, j, 0]

    value = np.power(np.abs(value), 2)
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(value))

    img = ax.imshow(value, cmap=plt.get_cmap("viridis"), norm=norm)
    ax.text(5, 5, "t = 0.000", bbox={"facecolor": "white", "pad": 10})
    cbar = fig.colorbar(img, ax=ax)

    def animation(idx):
        value = np.zeros((size, size), dtype="complex")

        for i in range(size):
            for j in range(size):
                value[i, j] = A[i, j, idx]

        value = np.power(np.abs(value), 2)
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(value))
        ax.text(5, 5, f"t = {delta_t*idx:e}", bbox={"facecolor": "white", "pad": 10})
        ax.set_xticks(np.linspace(0, M - 1, 6))
        ax.set_xticklabels([i / 10 for i in range(0, 12, 2)])
        ax.set_yticks(np.linspace(0, M - 1, 6))
        ax.set_yticklabels([i / 10 for i in range(10, -1, -2)])

        img.set_norm(norm)

        img.set_data(value)

        return img

    anim = FuncAnimation(
        fig,
        animation,
        interval=1000,
        frames=np.arange(0, length, 1),
        repeat=False,
        blit=0,
    )

    anim.save(filename, writer="ffmpeg", bitrate=-1, fps=30)


if __name__ == "__main__":
    newpath = r"animations"
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    # create_animation(
    #     "Two_slits_sigma_y_0.1", 200, 320, 2, newpath + "/two_slits_sigma_y_01.mp4"
    # )
    # create_animation(
    #     "No_slit_sigma_y_005", 200, 320, 0, newpath + "/No_slit_sigma_y_005.mp4"
    # )
    create_animation(
        "One_slit_sigma_y_0.2",
        201,
        80,
        1,
        newpath + "/one_slits_sigma_y_02.mp4",
        2.5 * 10 ** (-5),
    )
    create_animation(
        "Two_slits_sigma_y_0.2",
        201,
        80,
        2,
        newpath + "/two_slits_sigma_y_02.mp4",
        2.5 * 10 ** (-5),
    )
    create_animation(
        "Three_slits_sigma_y_0.2",
        201,
        80,
        3,
        newpath + "/three_slits_sigma_y_02.mp4",
        2.5 * 10 ** (-5),
    )
