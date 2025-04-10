# DFT-CES2
### Density functional theory in classical explicit solvents 2
doi.org/10.1021/jacsau.5c00176

## How to use patch files
### To edit the original source codes from Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) and quantum-espresso (QE)
For LAMMPS 11 Aug 2017,

&gt;&gt; pwd <br /> $USER/lammps-11Aug17

&gt;&gt; ls <br /> src etc

&gt;&gt; patch -p0 < lammps.patch

For LAMMPS 2 Aug 2023,

&gt;&gt; pwd <br /> $USER/lammps-2Aug2023

&gt;&gt; ls <br /> src etc

&gt;&gt; patch -p0 < lammps-2Aug2023.patch

<br />
For QE,

&gt;&gt; pwd <br /> $USER/qe-6.3

&gt;&gt; ls <br /> PP Modules PW etc

&gt;&gt; patch -p0 < mod.patch

&gt;&gt; patch -p0 < pp.patch

&gt;&gt; patch -p0 < pw.patch
