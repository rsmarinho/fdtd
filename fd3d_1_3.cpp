#include <cstdlib>
#include <iostream>
#include <armadillo>
#include <getopt.h>

using namespace std;
using namespace arma;


void print_usage() {
    printf("Usage: fdtd -k num -f num -t num -n num\n");
}

int main(int argc, char *argv[]) {
    int options = 0;
    int ke = -1, nsteps = -1, verbose = -1;
    double frequency_signal = -1.0, delta_t = -1.0;

    //Specifying the expected options
    //The two options l and b expect numbers as argument
    static struct option long_options[] = {
        {"length division",  required_argument, 0,  'k' },
        {"frequency",        required_argument, 0,  'f' },
        {"time division",    required_argument, 0,  't' },
        {"time steps",       required_argument, 0,  'n' },
        {"verbose",          no_argument      , 0,  'v' },
        {0,                  0,                 0,  0   }
    };

    int long_index =0;
    while ((options = getopt_long(argc, argv,"vk:f:t:n:",
                   long_options, &long_index )) != -1) {
        switch (options) {
             case 'v' : verbose = 0;
                 break;
             case 'k' : ke = atoi(optarg);
                 break;
             case 'f' : frequency_signal = atof(optarg);
                 break;
             case 't' : delta_t = atof(optarg);
                 break;
             case 'n' : nsteps = atoi(optarg);
                 break;
             default: print_usage();
                 exit(EXIT_FAILURE);
        }
    }
    if (ke == -1 || nsteps == -1 || frequency_signal == -1.0 || delta_t == -1.0) {
        print_usage();
        exit(EXIT_FAILURE);
    }

    //initiate variables
    vec ex(ke,fill::zeros); 
    vec hy(ke,fill::zeros);

    //int kc = 0;
    //double t0 = 40;
    //double spread =12;

    vec boundary_low(2,fill::zeros);
    vec boundary_high(2,fill::zeros);

    //create dielectric profile
    vec cb(ke,fill::ones);
    cb = 0.5 * cb;

    float epsilon = 4.0;
    int cb_start = int(ke * 0.5);
    cb.rows(cb_start, ex.n_elem - 1) =  cb.rows(cb_start, ex.n_elem - 1) / epsilon; 

    double pulse = 0.0;

    double ddx = 0.01; //Cell size
    double dt = ddx / delta_t; // Time step size

    if(verbose == 0)
        cout << "profile:\n" << cb.t() << endl;
    for(int i = 1; i <= nsteps; i++){
        for(int k = 1; k < ke; k++){
            ex[k] = ex[k] + cb[k] * (hy[k-1] - hy[k]);
            if(verbose == 0)
                cout << "Itr: " << i << " - d_profile :" << cb[k] << " - Ex[" << k << "]:" << ex[k] << " = hy[" << k <<  "-1]:" << hy[k-1]<< " - hy[" << k <<"]:" << hy[k]<< endl;
        }

        //pulse = exp(-0.5 * pow( (t0 - i) / spread, 2) );
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
 //           cout << "Itr: " << i << " Ex:" << ex.t();
            hy[k] = hy[k] + 0.5 * (ex[k] - ex[k+1]);
            if(verbose == 0)
                cout << "Itr: " << i << " - d_profile :" << cb[k] << " - Hy["<< k << "]:" << hy[k] << " = ex[" << k << "]:" << ex[k]<< " - ex[" << k << "+1]:" << ex[k+1] << endl;
        }
        if(verbose == 0)
            cout << "============" << endl;
    }
    //cout << "Ex:\n" << ex.t() << "\n";
    //cout << "Hy:\n" << hy.t() << "\n"

    vec line = linspace(0, 1, ke);
    // Save matrices E and H:

    line.save(hdf5_name("fdtd_ans.hdf5", "line"), hdf5_binary);
    ex.save(hdf5_name("fdtd_ans.hdf5", "ex", hdf5_opts::append) );
    hy.save(hdf5_name("fdtd_ans.hdf5", "hy", hdf5_opts::append) );
    //hy.save("Hy_mat.hdf5", hdf5_binary);
    return 0;
}
