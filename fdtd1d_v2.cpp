#include "fdtd1d.hpp"

using namespace std;
using namespace arma;

void fdtd1d(int ke, int nsteps, int verbose, double frequency_signal, double delta_t) {
    //initiate variables
    vec ex(ke,fill::zeros);
    vec dx(ke,fill::zeros);
    vec ix(ke,fill::zeros);
    vec hy(ke,fill::zeros);

    vec boundary_low(2,fill::zeros);
    vec boundary_high(2,fill::zeros);

    double ddx = 0.01; //Cell size
    double dt = ddx / delta_t; // Time step size

    double pulse = 0.0;

    //create dielectric profile
    float epsz = 8.854e-12;
    float epsr = 4.0;
    float sigma = 0.04;

    vec gax(ke,fill::ones);
    vec gbx(ke,fill::zeros);
    
    int k_start = int(ke * 0.5);    
   
    gax.rows(k_start, ex.n_elem - 1) = gax.rows(k_start, ex.n_elem - 1) / ( epsr + ( sigma * dt / epsz ) );
    gbx.rows(k_start, ex.n_elem - 1) = gbx.rows(k_start, ex.n_elem - 1) + ( sigma * dt / epsz ); 

    for(int i = 1; i <= nsteps; i++){
        // calculate Dx
        for(int k = 1; k < ke; k++){
            dx[k] = dx[k] + 0.5 * (hy[k - 1] - hy[k]);
        }
        
        // put a sinusoidal at the low end
        pulse = sin(2.0 * datum::pi * frequency_signal * dt * i);
        dx[1] = pulse + dx[1];
        
        // calculate Ex field from Dx
        for(int k = 1; k < ke; k++){
            ex[k] = gax[k] * ( dx[k] - ix[k] );
            ix[k] = ix[k] + ( gbx[k] * ex[k] );
        }

        // absorbing boundary conditions
        ex[0] = boundary_low[0];
        boundary_low = shift(boundary_low, -1);
        boundary_low[boundary_low.n_elem - 1] = ex[1];

        ex[ex.n_elem - 1] = boundary_high[0];
        boundary_high = shift(boundary_high, - 1);
        boundary_high[boundary_high.n_elem - 1] = ex[ke - 2];

        // calculate Hy field
        for(int k = 0; k < ke - 1; k++){
            hy[k] = hy[k] + 0.5 * (ex[k] - ex[k+1]);
        }
    }

    vec line = linspace(0, 1, ke);
    // Save matrices E and H:

    line.save(hdf5_name("fdtd_ans.hdf5", "line"), hdf5_binary);
    ex.save(hdf5_name("fdtd_ans.hdf5", "ex", hdf5_opts::append) );
    hy.save(hdf5_name("fdtd_ans.hdf5", "hy", hdf5_opts::append) );
    //hy.save("Hy_mat.hdf5", hdf5_binary);
}
