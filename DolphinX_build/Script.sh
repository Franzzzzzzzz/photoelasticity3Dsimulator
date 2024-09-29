#!/bin/bash

#Basix
if [ 1 -eq 0 ]
then 
git clone https://github.com/FEniCS/basix.git
cd basix/cpp
cmake -DCMAKE_BUILD_TYPE=Release -B build-dir -S .
cmake --build build-dir
sudo cmake --install build-dir
cd ../..

# HDF5
wget 'https://objects.githubusercontent.com/github-production-release-asset-2e65be/258591100/1a8a2f09-544a-465f-a7bf-ed04ad7eb8f0?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=releaseassetproduction%2F20240929%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20240929T120204Z&X-Amz-Expires=300&X-Amz-Signature=c52a3bbd935c5cfb15882f616b1b60432f697ebd03028e23c1d8a5134a14e53c&X-Amz-SignedHeaders=host&response-content-disposition=attachment%3B%20filename%3Dhdf5-1.14.4-2.tar.gz&response-content-type=application%2Foctet-stream' -nc -O hdf5-1-14-4-2.tar.gz
tar -xf hdf5-1-14-4-2.tar.gz
cd hdf5-1.14.4-2
CC=/usr/bin/mpicc ./configure --enable-parallel --prefix=$HOME/local
make
make install 
cd ..


# GKlib
git clone https://github.com/KarypisLab/GKlib.git
cd GKlib
make config 
make
cd ..

#METIS
git clone https://github.com/KarypisLab/METIS.git
cd METIS
make config shared=1 cc=mpicc prefix=$HOME/local
make install
cd ..

#ffcx 
git clone https://github.com/FEniCS/ffcx.git
cd ffcx
cmake -B build-dir -S cmake/
cmake --build build-dir
sudo cmake --install build-dir
cd ..
fi
# dolphinx
git clone https://github.com/FEniCS/dolfinx.git
cd dolfinx/cpp
mkdir build 
cd build 
xterm -e 'echo "Variable to setup for cmake:
DOLFINX_UFCX_PYTHON = OFF
HDF5_LIBRARIES = /home/franz/local/lib/libhdf5.so
HDF5_INCLUDE_DIRS = /home/franz/local/include/
PARMETIS_INCLUDE_DIRS = /home/franz/local/include
PARMETIS_LIBRARY = /home/franz/local/lib"; bash' &
ccmake ..
