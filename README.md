# Sparsehub

Sparsehub solves the linear advection equation for multiple materials at once.

## Build Instructions

    git clone git@gitlab.lanl.gov:charest/sparsehub
    cd sparsehub
    mkdir build
    cd build
    cmake ../.
    make

## Running Problems

To run a sampe problem:

    ./sparsehub ../examples/input.inp

You can use different sparse data represenations us=ing the '-l' command line argument.  For example, 

    ./sparsehub ../examples/input.inp -l 2

will run the problem using compressed row storage.  Currently, dense (0), fixed-bandwith (1), and compressed row (2) storage is available.

## Command Line Arguments

    ./sparsehub -h

    help: sparsehub [-h] [-l LAYOUT] [FILENAME]

    Sparsehub solves the advection equation using a sparse matrix representation of the data.

    Options:
      -h              Display this help message.
      -l LAYOUT       Specify the sparse layout type [0=dense (default), 1=fixed bandwith, 2=compressed row].

    Arguments:
      FILENAME        The input file name.
