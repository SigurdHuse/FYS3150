import numpy as np
import matplotlib.pyplot as plt

main = np.loadtxt("numeric_solution_n_1000.txt")
plt.plot(main[:,0], main[:,1])
plt.show()