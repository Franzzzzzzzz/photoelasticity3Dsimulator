#!/bin/bash

if [ 1 -eq 0 ]
then 
#Basix
git clone https://github.com/FEniCS/basix.git
cd basix/cpp
cmake -DCMAKE_BUILD_TYPE=Release -B build-dir -S .
cmake --build build-dir
sudo cmake --install build-dir
cd ../..

# HDF5
cp ~/Téléchargements/hdf5-1.14.4-2.tar.gz .
tar -xf hdf5-1.14.4-2.tar.gz
cd hdf5-1.14.4-2
CC=/usr/bin/mpicc ./configure --enable-parallel --prefix=$HOME/local
make -j 4
make install 
cd ..

# GKlib
git clone https://github.com/KarypisLab/GKlib.git
cd GKlib
make config 
make -j4
make install
cd ..

#METIS
git clone https://github.com/KarypisLab/METIS.git
cd METIS
make config shared=1 cc=mpicc prefix=$HOME/local
make -j4
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

# PETSc & SLEPc (recommended
sudo apt install petsc-dev slepc-dev

# dolphinx
git clone https://github.com/FEniCS/dolfinx.git
cd dolfinx/cpp
mkdir build 
cd build 
xterm -e 'echo "Variable to setup for cmake:
DOLFINX_UFCX_PYTHON = OFF
HDF5_LIBRARIES = /home/franz/local/lib/libhdf5.so
HDF5_INCLUDE_DIRS = /home/franz/local/include
_REQUIRE_PARMETIS = ON
PARMETIS_INCLUDE_DIRS = /home/franz/local/include
PARMETIS_LIBRARY = /home/franz/local/lib"; bash' &
ccmake ..
make -j4
make install 
cd ..

### Setting up the python environment 

# We are gonna need the python part for ffcx ...
conda create -n fenicsx-source python=3.10 #create a new environment
conda activate fenicsx-source
export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/ #for pip to find gfortran in blas ........
cd basicx/python 
pip install .
cd ../..

git clone https://github.com/FEniCS/ufl.git
cd ufl 
pip install .
cd ..

cd ffcx
pip install .
cd ..

pip install pyvista petsc4py

cd dolphinx
source /home/franz/local/lib/dolfinx/dolfinx.conf
pip install -r build-requirements.txt
pip install --check-build-dependencies --no-build-isolation .


if [ 1 -eq 0 ]
then 
fi
