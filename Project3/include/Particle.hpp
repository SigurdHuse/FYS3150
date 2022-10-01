#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
public:
    double q_, m_;
    arma::vec r_ = arma::vec(3), v_ = arma::vec(3);
    Particle(double q, double m, arma::vec r, arma::vec v);
};

#endif