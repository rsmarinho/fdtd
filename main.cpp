#include <iostream>
#include <getopt.h>

#include "fdtd1d.hpp"

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

    fdtd1d(ke, nsteps, verbose, frequency_signal, delta_t);

    return 0;
}
