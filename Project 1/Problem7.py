import numpy as np
import matplotlib.pyplot as plt


for n in [10, 100, 1000, 10000]:
    main = np.loadtxt(f"numeric_solution_n_{n}.txt", delimiter = ", ", usecols = range(2))
    plt.plot(main[:,0], main[:,1], label = f"N = {n}")  

main = np.loadtxt(f"exact_solution.txt", usecols = range(2))
plt.plot(main[:,0], main[:,1], label = f"Exact solution")
plt.xlabel("steps[x_i]")
plt.ylabel("y")
plt.title("analytical solution compared with numeric solution")
plt.grid()
plt.legend()
plt.show()


#print(main[:,0])
#print('................................................................')
#print(main[:,1])





