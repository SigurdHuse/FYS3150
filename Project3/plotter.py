import numpy as np
import matplotlib.pyplot as plt

main = np.loadtxt("one_particle_n_10000.txt", delimiter=",", usecols=range(1))
t, x, y, z = main[::4], main[1::4], main[2::4], main[3::4]

plt.plot(x, y)
plt.show()
