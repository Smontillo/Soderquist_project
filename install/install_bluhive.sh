module load armadillo/10.7.1
module load gcc/11.2.0/b1
module load cmake/3.17.0
module load hdf5/1.8.17/b1
module load json11/1.0.0
module load lapack/3.7.1/b1
module load openblas/0.2.20/b1

export CC=gcc
export CXX=g++

rm -rf build
mkdir build
cd build
cmake ..
make rhot
cd ..
