#include <cstdlib>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, const char **argv) {
    //initiate variables
    int ke = 200;

    vec ex(ke,fill::zeros); 
    vec hy(ke,fill::zeros);

    int kc = int(ke/2);
    double t0 = 40;
    double spread =12;

    vec boundary_low(2,fill::zeros);
    vec boundary_high(2,fill::zeros);

    int nsteps = atoi(argv[1]);
    double pulse = 0.0;

    for(int i = 0; i < nsteps; i++){
        for(int k = 1; k <= ke; k++){
            ex[k] = ex[k] + 0.5 * (hy[k-1] - hy[k]);
        }

        pulse = exp(-0.5 * pow( (t0 - i) / spread, 2) );
        ex[kc] = pulse;

        ex[0] = boundary_low[0];
        boundary_low = shift(boundary_low, -1);
        boundary_low[boundary_low.n_elem - 1] = ex[1];

        ex[ex.n_elem - 1] = boundary_high[0];
        boundary_high = shift(boundary_high, -1);
        boundary_high[boundary_high.n_elem - 1] = ex[ke - 2];

        for(int k = 0; k < ke; k++){
            hy[k] = hy[k] + 0.5 * (ex[k] - ex[k+1]);
        }
    }
    cout << "Ex:\n" << ex.t() << "\n";
    cout << "Hy:\n" << hy.t() << "\n";

    // Save matrices E and H:
    ex.save("Ex_mat.hdf5", hdf5_binary);
    hy.save("Hy_mat.hdf5", hdf5_binary);
    return 0;
}
