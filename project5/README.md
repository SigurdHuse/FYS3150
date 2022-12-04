# Simulating a Wave packet with zero or more slits
The schrodinger equation can be used to describe the movement of a wave packet. To able to simulate such a wave packet it is placed inside a square box with side length 1 which dirichlet boundary conditions. Using the Crank-Nicolson method it is then possible to numerically solve the schrodinger equation, and simulate the movement of a wave packet. Additionally, it is possible to add a wall with slits inside the box, and simulate how this affects the wave packet's movement. 

## How to build
Run command:  
`make`

## How to run main.exe
### Generate data from report
Run command  `./main.exe`.

### Run your own simulation
Run command:  
`./main.exe <M> <T> <steps> <sigma_x> <x_c> <p_x> <sigma_y> <y_c> <p_y> <v0> <Name of simulation>` 

Where:  
* M = Side length of system
* T = Total time
* steps = Number of time steps
* sigma_x = Initial width of wave packet in x-direction
* x_c = Initial x-position of wave packet
* p_x = Wave packet momenta in x-direction
* sigma_y = Initial width of wave packet in y-direction
* y_c = Initial y-position of wave packet
* p_y = Wave packet momenta in y-direction
* v0 = Potential value the barrier has


### How to change barrier configuration
Can be done by changing `config.txt`, explanations of parameters are in file.

## Run tests
Run command:   
`c++ test.cpp src/system.cpp -I include/ -larmadillo -O2 -o test.exe && ./test.exe`

## Generate plots and compute estimated values from report
Run command  `python3 plotter.py`.  

## Generate animations of data from report
Run command `python3 animation.py`.  

## Dependecies
The code was run on Ubuntu 20.04.5 LTS, with a Intel® Core™ i7-9750H CPU @ 2.60GHz × 12. The dependecies and how to install are:  
* Armadillo version: 9.800.4 (Horizon Scraper), `sudo apt-get install libarmadillo-dev`
* Python3 3.8.10, `sudo apt install python3.8`
* c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0, `sudo apt install build-essential`
* Numpy version 1.23.3, `pip install numpy`
* Matplotlib version 3.4.3, `pip install matplotlib`