# Simulating the The Ising model using the Metropolis algorithm

## How to build
### With parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -fopenmp -O2`

### Without parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -O2`

## How to run main.exe

### One temperature
Run command:  
`./main.exe <Length of grid> <Temperature of system> <Number of cycles> <Number of Markoc Chains> <Bool indicating random or ordered initial state>`

### Range of temperatures
Run command:  
`./main.exe <Length of grid> <Start temperature> <End temperature> <Increase in temperature for each simulation> <Number of cycles> <Number of Markoc Chains> <Bool indicating random or ordered initial state>`

## Run tests
Run command:   
`c++ test.cpp src/system.cpp -I include/ -larmadillo -O2 -o test.exe && ./test.exe`