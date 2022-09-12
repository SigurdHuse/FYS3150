import numpy as np
import matplotlib.pyplot as plt



for n in [10, 100, 1000, 10000]:
    exact = np.loadtxt(f"exact_solution_n_{n}.txt", delimiter = ", ", usecols = range(2)) 
    numeric = np.loadtxt(f"numeric_solution_n_{n}.txt", delimiter = ", ", usecols = range(2))

    diff = np.log10(np.absolute(exact[:,1] - numeric[:,1]))
    plt.plot(exact[1:-1,0], diff[1:-1], label = f"N = {n}")  

plt.legend()
plt.grid()
plt.title('Absolute error')
plt.xlabel("x-values")
plt.ylabel("log 10")
plt.show()
#plt.savefig("problem_8a.pdf")
plt.clf()



for n in [10, 100,1000,10000]:
    exact = np.loadtxt(f"exact_solution_n_{n}.txt", delimiter = ", ", usecols = range(2), dtype = "float128") 
    numeric = np.loadtxt(f"numeric_solution_n_{n}.txt", delimiter = ", ", usecols = range(2), dtype = "float128")
    diff = np.log10(np.absolute((exact[1:-1,1] - numeric[1:-1,1])/exact[1:-1,1]))
    #max_rel.append(np.max(10**diff))
    plt.plot(exact[1:-1,0], diff, label = f"N = {n}")


plt.grid()
plt.legend()
plt.title("relative error")
plt.xlim([0,1])
plt.xlabel("x-values")
plt.ylabel("log 10")   
#plt.show() 


"""
max_rel = []    
#plt.savefig("problem_8b.pdf")

n_values = [10**i for i in range(1,8)]

for n in n_values:
    exact = np.loadtxt(f"exact_solution_n_{n}.txt", delimiter = ", ", usecols = range(2)) 
    numeric = np.loadtxt(f"numeric_solution_n_{n}.txt", delimiter = ", ", usecols = range(2))
    diff = np.abs(((exact[1:-1,1] - numeric[1:-1,1])/exact[1:-1,1]))
    max_rel.append(np.max(diff))


for n, val in zip(n_values, max_rel):
    print(f"n = {n :8}  max relative error = {val}")
"""






