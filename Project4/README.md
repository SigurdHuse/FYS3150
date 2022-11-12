# Simulating the The Ising model using the Metropolis algorithm

## How to build
### With parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -fopenmp -O2`

### Without parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -O2`

## How to run main.exe
## Generate data from report
Run command 
`./main.exe -a`

### One temperature
Run command:  
`./main.exe <Length of grid> <Temperature of system> <Number of cycles> <Number of Markov Chains> <Bool indicating random or ordered initial state>`

### Range of temperatures
Run command:  
`./main.exe <Length of grid> <Start temperature> <End temperature> <Increase in temperature for each simulation> <Number of cycles> <Number of Markov Chains> <Bool indicating random or ordered initial state>`

## Run tests
Run command:   
`c++ test.cpp src/system.cpp -I include/ -larmadillo -O2 -o test.exe && ./test.exe`

## Generate plots
Run command  `python3 plotter.py`  

## Dependecies
The code was run on Ubuntu 20.04.5 LTS, with a Intel® Core™ i7-9750H CPU @ 2.60GHz × 12. The dependecies and how to install are:  
* Armadillo version: 9.800.4 (Horizon Scraper), `sudo apt-get install libarmadillo-dev`
* Python3 3.8.10, `sudo apt install python3.8`
* c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0, `sudo apt install build-essential`