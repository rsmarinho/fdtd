#ifndef _FDTD_1D_HPP_
#define _FDTD_1D_HPP_

#include <cstdlib>
#include <iostream>

// #define H5_USE_18_API
// #include <hdf5.h>

// #define ARMA_USE_HDF5
#include <armadillo>

void fdtd1d(int ke, int nsteps, int verbose, double frequency_signal, double delta_t);

#endif