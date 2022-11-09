# Simulating the The Ising model using the Metropolis algorithm

## How to build

### With parallelization
`c++ main.cpp -o main.exe -larmadillo -fopenmp -O2`

### Without parallelization
`c++ main.cpp -o main.exe -larmadillo -O2`

## How to run main.exe

### One temperature
`./main.exe <Length of grid> <Temperature of system> <Number of cycles> <Number of Markoc Chains> <Bool indicating random or ordered initial state>`

### Range of temperatures