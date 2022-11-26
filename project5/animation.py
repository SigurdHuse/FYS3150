import pyarma as pa
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import numpy as np

plt.rcParams["animation.ffmpeg_path"] = "/usr/bin/ffmpeg"


def create_animation(Name, M, time_steps, slits, filename):
    A = pa.cx_cube()

    fig, ax = plt.subplots()

    A.load(Name + f"_M_{M}_dt_{time_steps}_slits_{slits}.bin")

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

    img = ax.imshow(
        value,
        cmap=plt.get_cmap("viridis"),
        norm=norm,
    )

    cbar = fig.colorbar(img, ax=ax)

    def animation(idx):
        # Normalize the colour scale to the current frame?
        value = np.zeros((size, size), dtype="complex")
        for i in range(size):
            for j in range(size):
                value[i, j] = A[i, j, idx]

        value = np.power(np.abs(value), 2)
        print(idx)
        # print(np.sum(value))
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(value))

        img.set_norm(norm)

        # Update z data
        img.set_data(value)

        return img

    anim = FuncAnimation(
        fig,
        animation,
        interval=1000,
        frames=np.arange(0, length - 1, 1),
        repeat=False,
        blit=0,
    )
    # plt.show()
    anim.save(filename, writer="ffmpeg", bitrate=-1, fps=60)


if __name__ == "__main__":
    create_animation("T008Bigsigmay", 200, 320, 2, "twoslits.mp4")
