# fdtd

## compile with
make FILE="fd3d_1_5.cpp"

## run with
clear; ./fdtd -k 200 -f 700e6 -t 6e8 -n 500

## or multiple runs
for i in {100..2000..50}; do clear; ./fdtd -k 200 -f 700e6 -t 6e8 -n 2000; sleep 0.5; done
