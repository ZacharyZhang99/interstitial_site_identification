# interstitial_site_identification

Copyright (2025) Zuoyong (Zachary) Zhang. This repo shares the python scripts for interstitial site identification in alloy systems.

Requirements:
1. lammps with python package
2. python 3.xx or higher version

Input files:
1. M.txt, M represents the metal element of solvent, such as Al, Ni, Mg, and Cu, etc. This file includes four columns which are atom id, x, y, and z, respectively.
2. potential file
3. in.min, lammps script for site screening.
4. M.data, lammps data file, please use the format in examples.

Compile lammps with python:
1. create a build folder in your lammps-xxxxxx folder

2. compile lammps using cmake command:

cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -DCMAKE_C_COMPILER=/usr/bin/mpicc -D BUILD_SHARED_LIBS=yes -D LAMMPS_MACHINE=mpi -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DPKG_PYTHON=yes -D DOWNLOAD_KIM=yes -D PKG_KIM=on -D LAMMPS_EXCEPTIONS=yes -D PKG_MANYBODY=yes ../cmake
make
make install-python
