# Simulating a Wave packet with zero or more slits

## How to build
Run command:  
`make`

## How to run main.exe
### Generate data from report
Run command 
`./main.exe`

### Run your own simulation
Run command:  
`./main.exe <side length of system> <Total time> <Number of time steps> <Initial width of wave packet in x-direction> <Initial x-position of wave packet> <Wave packet momenta in x-direction> <Initial width of wave packet in y-direction> <Initial y-position of wave packet> <Wave packet momenta in y-direction> <Value of barrier> <Name of simulation>` 

### How to change barrier configuration
Can be done by changing `config.txt`, explanations of parameters are in file.

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