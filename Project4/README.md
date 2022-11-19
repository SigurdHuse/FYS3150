# Simulating the The Ising model using the Metropolis algorithm

The Ising model consists of a static square grid filled with "particles", which can only interact with the particles to the: left, right, up and down of it. The model is used to describe ferromagnetic properties, and uses quantized spin. It consits of a LxL-grid filled with particles which have either +1 or -1 spin. To simulate such a system the program in this repository uses the Metropolis algorithm, which simulates a Markov chain describing the system. Running enough cycles of simulation the Chain reaches equilibrium, and we get a good approximations of the distribution of the states of the system.

## How to build
### With parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -fopenmp -O2`

### Without parallelization
Run command:  
`c++ main.cpp src/system.cpp -I include/ -o main.exe -larmadillo -O2`

## How to run main.exe
If running with parallelization, make sure to specifiy number of cores to be used by running command `export OMP_NUM_THREADS= <number of cores>`
## Generate data from report
Run command 
`./main.exe`

### One temperature
Run command:  
`./main.exe <Length of grid> <Temperature of system> <Number of cycles> <Number of Markov Chains> <Bool indicating random or ordered initial state>`

### Range of temperatures
Run command:  
`./main.exe <Length of grid> <Start temperature> <End temperature> <Increase in temperature for each simulation> <Number of cycles> <Number of Markov Chains> <Bool indicating random or ordered initial state>`

### Time running time
Run command  
`time ./main.exe <Desired commands>`

## Run tests
Run command:   
`c++ test.cpp src/system.cpp -I include/ -larmadillo -O2 -o test.exe && ./test.exe`

## Generate plots and compute estimated values from report
Run command  `python3 plotter.py`  

## Dependecies
The code was run on Ubuntu 20.04.5 LTS, with a Intel® Core™ i7-9750H CPU @ 2.60GHz × 12. The dependecies and how to install are:  
* Armadillo version: 9.800.4 (Horizon Scraper), `sudo apt-get install libarmadillo-dev`
* Python3 3.8.10, `sudo apt install python3.8`
* c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0, `sudo apt install build-essential`
* Numpy version 1.23.3, `pip install numpy`
* Matplotlib version 3.4.3, `pip install matplotlib`
* Scipy version 1.91, `pip install scipy`