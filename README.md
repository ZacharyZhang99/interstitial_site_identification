# interstitial_site_identification

Copyright © 2025 Zuoyong (Zachary) Zhang. This repository provides Python scripts designed for identifying interstitial sites in alloy systems.

Requirements:
LAMMPS with the Python package installed
Python 3.x or later

Input Files:
M.txt: Represents the metal element of the solvent (e.g., Al, Ni, Mg, Cu, etc.). This file contains four columns: atom ID, x-coordinate, y-coordinate, and z-coordinate.
Potential file: Required for defining interactions in the system.

in.min: A LAMMPS script used for screening interstitial sites.

M.data: A LAMMPS data file. Please follow the format provided in the examples.


Compiling LAMMPS with Python:
Create a build folder within your lammps-xxxxxx directory.
Compile LAMMPS using the following CMake command:


cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -DCMAKE_C_COMPILER=/usr/bin/mpicc -D BUILD_SHARED_LIBS=yes -D LAMMPS_MACHINE=mpi -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DPKG_PYTHON=yes -D DOWNLOAD_KIM=yes -D PKG_KIM=on -D LAMMPS_EXCEPTIONS=yes -D PKG_MANYBODY=yes ../cmake
make
make install-python
