This repository contains my own public branches.

To get the official plumed versions, please refer to the repository at http://github.com/plumed/plumed2

Below you can see the content of the local branches. These branches might be merged to official plumed repository
as soon as they are ready for production.

* _gio-master_: This README file.

* _v2.0-hrex_: An implementation of Hamiltonian replica exchange for GROMACS.
There is an alternate GROMACS patch and a script to edit topologies. Both are preliminary and might not work with
all GROMACS settings. If you with to use it, I suggest you to get in touch with me (bussi@sissa.it)

* _v2.0-plumed1_: An extra module (src/plumed1) which, in combination with the development version of PLUMED 1.3,
can be used to emulate collective variables and biases in PLUMED 2. Notice that not all the features
are supported (e.g. replica exchange is not supported), and that efficiency is poor.
Mostly to be used to test when porting variables.

These branches are not maintained since they were already merged to official plumed repository:

* _v2.0-gpu-balance_: Modified gromacs patch which allows for a better load balance with GPUs
