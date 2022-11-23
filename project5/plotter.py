import pyarma as pa
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

A = pa.cx_cube()
B = pa.cx_mat()

fig, ax = plt.subplots()

A.load("Data_M_200_dt_320.bin")

value = np.zeros((200, 200), dtype="complex")

for idx in np.arange(0, 320, 20):
    for i in range(200):
        for j in range(200):
            value[i, j] = A[i, j, idx]

    new_value = np.power(np.abs(value), 2)
    print(np.sum(new_value))
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(new_value))

    im = plt.imshow(new_value, norm=norm, cmap=plt.get_cmap("viridis"))
    plt.colorbar()
    plt.show()
