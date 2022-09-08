import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('exact_solution.txt')

plt.plot(data[:,0], data[:,1], '-', label="Exact Solution")
plt.xlabel('x_i')
plt.ylabel('y')
plt.grid()
plt.title(f'Exact solution for the Poisson equation in 1.D for n = {len(data[:,0])}')
plt.legend(loc='upper right')

plt.savefig("figure_problem2.pdf")




























