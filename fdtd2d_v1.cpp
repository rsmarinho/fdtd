#include "fdtd1d.hpp"

using namespace std;
using namespace arma;

void fdtd1d(int ke, int nsteps, int verbose, double frequency_signal, double delta_t) {
    // initiate variables
    int ie = ke;
    int je = ke;
    int ic = int(ie/2);
    int jc = int(je/2);
    
    mat ez(ie, je, fill::zeros);
    mat dz(ie, je, fill::zeros);
    mat hx(ie, je, fill::zeros);
    mat hy(ie, je, fill::zeros);
    mat ihx(ie, je, fill::zeros);
    mat ihy(ie, je, fill::zeros);

    double ddx = 0.01; //Cell size
    double dt = ddx / delta_t; // Time step size

    double pulse = 0.0;

    float t0 = 50;
    float spread = 12;
    
    // create dielectric profile
    mat gaz(ie, je, fill::ones);

    // calculate the PML parameters
    vec gi2(ie, fill::ones);
    vec gi3(ie, fill::ones);
    vec fi1(ie, fill::zeros);
    vec fi2(ie, fill::ones);
    vec fi3(ie, fill::ones);

    vec gj2(ie, fill::ones);
    vec gj3(ie, fill::ones);
    vec fj1(ie, fill::zeros);
    vec fj2(ie, fill::ones);
    vec fj3(ie, fill::ones);

    // create the PML as described
    int npml = 8;
    int xnum, xd, xxn, xn;
    for(int n = 0; n < npml; n++){
        xnum = npml - n;
        xd = npml;
        xxn = xnum / xd;
        xn = 0.33 * pow(xxn, 3);

        gi2[n] = 1 / (1 + xn);
        gi2[ie - 1 - n] = 1 / (1 + xn);
        gi3[n] = (1 - xn) / (1 + xn);
        gi3[ie - 1 - n] = (1 - xn) / (1 + xn);

        gj2[n] = 1 / (1 + xn);
        gj2[je - 1 - n] = 1 / (1 + xn);
        gj3[n] = (1 - xn) / (1 + xn);
        gj3[je - 1 - n] = (1 - xn) / (1 + xn);

        xxn = (xnum - 0.5) / xd;
        xn = 0.33 * pow(xxn, 3);

        fi1[n] = xn;
        fi1[ie - 2 - n] = xn;
        fi2[n] = 1 / (1 + xn);
        fi2[ie - 2 - n] = 1 / (1 + xn);
        fi3[n] = (1 - xn) / (1 + xn);
        fi3[ie - 2 - n] = (1 - xn) / (1 + xn);

        fj1[n] = xn;
        fj1[je - 2 - n] = xn;
        fj2[n] = 1 / (1 + xn);
        fj2[je - 2 - n] = 1 / (1 + xn);
        fj3[n] = (1 - xn) / (1 + xn);
        fj3[je - 2 - n] = (1 - xn) / (1 + xn);
    }

    for(int i = 1; i <= nsteps; i++){
        
        // calculate Dz
        for(int k = 1; k < je - 1; k++){
            for(int j = 1; j < ie - 1; j++){
                dz(j,k) = gi3(j) * gj3(k) * dz(j,k) + 
                          gi2(j) * gj2(k) * 0.5 * 
                          (hy(j,k) - hy(j-1,k) - 
                           hx(j,k) + hx(j,k-1));
            }
        }
        
        // put a sinusoidal at the low end
//        pulse = exp(-0.5 * pow(((t0 - i) / spread), 2) );
        pulse = sin(2.0 * datum::pi * frequency_signal * dt * i);
        dz(ic-10,jc-10) = pulse;

//        ez = gaz * dz; // calculate the Ez field from Dz
    
        // calculate Ez field from Dz
        for(int k = 1; k < je - 1; k++){
           for(int j = 1; j < ie - 1; j++){
                ez(j,k) = gaz(j,k) *  dz(j,k);
            }
        }
    
        // calculate Hx field
        for(int k = 0; k < je - 1; k++){
            for(int j = 0; j < ie - 1; j++){
                double acurl_e = ez(j,k) - ez(j,k + 1);
                ihx(j,k) = ihx(j,k) + acurl_e;
                hx(j,k) = fj3(k) * hx(j,k) + fj2(k) * 
                          (0.5 * acurl_e + fi1(j) * ihx(j, k));
            }
        }
    
        // calculate Hy field
        for(int k = 0; k < je - 1; k++){
            for(int j = 0; j < ie - 1; j++){
                double acurl_e = ez(j,k) - ez(j + 1,k);
                ihy(j,k) = ihy(j,k) + acurl_e;
                hy(j,k) = fi3(j) * hy(j,k) - fi2(j) * 
                          (0.5 * acurl_e + fj1(k) * ihy(j, k));
            }
        }
    }

//    mat plane = linspace(0, 1, ke);
    // Save matrices E and H:

    ez.save(hdf5_name("fdtd_ans.hdf5", "ez"), hdf5_binary);
    hx.save(hdf5_name("fdtd_ans.hdf5", "hx", hdf5_opts::append) );
    hy.save(hdf5_name("fdtd_ans.hdf5", "hy", hdf5_opts::append) );
}
