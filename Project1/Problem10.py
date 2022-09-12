import numpy as np
import matplotlib.pyplot as plt


main = np.loadtxt(f"Run_times.txt", delimiter = ", ", usecols = range(2))

n_values = [10**i for i in range(1, 7)]

plt.plot(n_values, main[:,0], 'bo', label = "Not specialized")
plt.plot(n_values, main[:,1], 'ro', label = "Specialized")
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.legend()
plt.title("Plot of run time for specialized and general algorithm")
plt.xlabel("n")
plt.ylabel("Log of run time")
plt.savefig("Problem10.pdf")