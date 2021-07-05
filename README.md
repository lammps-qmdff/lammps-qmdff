**paper-2021-data** directory contains examples of input files for the use of QMDFF force field in LAMMPS.

**lammps-modified** directory contains a modified version of LAMMPS supporting QMDFF. This is not an active repository, it is just a copy of the source code of the LAMMPS master branch from 30.05.2017, with modifications. To install it, one can use make utility and follow standard instructions of LAMMPS installation.

Quick installation guide that activates modifications with all required dependencies:

Enter lammps-qmdff directory, change to **"src"** subdirectory, enable packages:

<code>
make yes-USER-EXCITED yes-USER-MISC yes-MOLECULE yes-KSPACE yes-MISC yes-RIGID
</code>

Optionally, include any other packages you need for the simulations.

Build LAMMPS (for example, on 12 CPUs with default MPI options):

<code>
make -j12 mpi
</code>

