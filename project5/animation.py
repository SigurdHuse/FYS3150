import pyarma as pa
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import numpy as np


def create_animation(M, time_steps, filename):
    A = pa.cx_cube()

    fig, ax = plt.subplots()

    A.load(f"Data_M_{M}_dt_{time_steps}.bin")
    length = time_steps
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

    def animation(idx):
        # Normalize the colour scale to the current frame?
        value = np.zeros((size, size), dtype="complex")
        for i in range(size):
            for j in range(size):
                value[i, j] = A[i, j, idx]

        value = np.power(np.abs(value), 2)
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
        frames=np.arange(0, length, 1),
        repeat=False,
        blit=0,
    )

    anim.save(filename, writer="ffmpeg", bitrate=-1, fps=60)


if __name__ == "__main__":
    create_animation(200, 320, "animation1.mp4")
