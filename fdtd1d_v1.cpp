#include "fdtd1d.hpp"

using namespace std;
using namespace arma;

void fdtd1d(int ke, int nsteps, int verbose, double frequency_signal, double delta_t) {
    //initiate variables
    vec ex(ke,fill::zeros); 
    vec hy(ke,fill::zeros);

    vec boundary_low(2,fill::zeros);
    vec boundary_high(2,fill::zeros);

    double ddx = 0.01; //Cell size
    double dt = ddx / delta_t; // Time step size

    double pulse = 0.0;

    //create dielectric profile
    float epsz = 8.854e-12;
    float epsilon = 4.0;
    float sigma = 0.002;//0.04;

    vec ca(ke,fill::ones);
    vec cb(ke,fill::ones);
    cb = 0.5 * cb;
    int cb_start = int(ke * 0.5);    
   
    double eaf = dt * sigma / ( 2 * epsz * epsilon );
    ca.rows(cb_start, ex.n_elem - 1) = ca.rows(cb_start, ex.n_elem - 1) * ( 1 - eaf ) / ( 1 + eaf );
    cb.rows(cb_start, ex.n_elem - 1) =  cb.rows(cb_start, ex.n_elem - 1) / ( epsilon * (1 + eaf) ); 

    if(verbose == 0)
        cout << "profile:\n" << cb.t() << endl;
    for(int i = 1; i <= nsteps; i++){
        for(int k = 1; k < ke; k++){
            ex[k] = ca[k] * ex[k] + cb[k] * (hy[k - 1] - hy[k]);
            if(verbose == 0)
                cout << "Itr: " << i << " - d_profile :" << cb[k] << " - Ex[" << k << "]:" << ex[k] << " = hy[" << k <<  "-1]:" << hy[k-1]<< " - hy[" << k <<"]:" << hy[k]<< endl;
        }

        pulse = sin(2.0 * datum::pi * frequency_signal * dt * i);
        ex[1] = pulse + ex[1];

        ex[0] = boundary_low[0];
        boundary_low = shift(boundary_low, -1);
        boundary_low[boundary_low.n_elem - 1] = ex[1];

        ex[ex.n_elem - 1] = boundary_high[0];
        boundary_high = shift(boundary_high, - 1);
        boundary_high[boundary_high.n_elem - 1] = ex[ke - 2];

        if(verbose == 0)
            cout << "------------" << endl; 

        for(int k = 0; k < ke - 1; k++){
            hy[k] = hy[k] + 0.5 * (ex[k] - ex[k+1]);
            if(verbose == 0)
                cout << "Itr: " << i << " - d_profile :" << cb[k] << " - Hy["<< k << "]:" << hy[k] << " = ex[" << k << "]:" << ex[k]<< " - ex[" << k << "+1]:" << ex[k+1] << endl;
        }
        if(verbose == 0)
            cout << "============" << endl;
    }

    vec line = linspace(0, 1, ke);
    // Save matrices E and H:

    line.save(hdf5_name("fdtd_ans.hdf5", "line"), hdf5_binary);
    ex.save(hdf5_name("fdtd_ans.hdf5", "ex", hdf5_opts::append) );
    hy.save(hdf5_name("fdtd_ans.hdf5", "hy", hdf5_opts::append) );
    //hy.save("Hy_mat.hdf5", hdf5_binary);
}
