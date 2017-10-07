Welcome to the Truspace developer documenation!

This site contains all the information you need to modify or extend Truspace.

For compilation and user-guide, please refer to http://bitbucket.org/truspacedevelopers/truspace


### Code design and structure

TruSpace is written in C++ using and Object Oriented design pattern.

The available classes are

-   \link parameters \endlink : class that handles the input reading, the calculation
    of the structure constants and the analysis of symmetry subspaces.
-   \link worker \endlink : class that takes care of running the different parts of the
    code after all parameters are set
-   \link states \endlink: class that calculates and holds all the informations
    regarding the computation basis
-   chirME : class that calculates and holds all the information about
    matrix element of the fields on the computational basis
-   \link filenames \endlink : helper class to manage output and swap files names
-   Hamiltonian : class that handles the computation and the matrix form
    of the TSA Hamiltonian of the perturbed model
-   \link diag \endlink: helper class that contains Numerical Recipies diagonalization
    routines for long double (required to calculate the conformal
    computational basis)
-   tensor2D, tensor3D, tensor4D, tensor5D, tensor6D : helper template class to handle 2D, 3D, 4D, 5D and
    6D tensors. It ensures memory contiguity and automatic deallocation
    when out of scope, without forcing compilation-time sizes.
-   NRG : model-agnostic full diagonalization, `NRG` and `Sweep` routine

Moreover, the are a series of helper functions

-   \link config.h \endlink Configuration file holding version number, some hard-coded limits that effect array dimensions, constants and macros 
-   \link block.h \endlink Contains a routine to identify symmetry subspace by the analysis of the interaction matrix. 
-   \link characters.h \endlink Set of helper functions to calculate integers partitions and manage GiNaC computation of character expansion 
-   \link cpt.h \endlink Set of helper function to calculate structure constant of Minimal CFTs 
-   \link gamma.h \endlink Computation of gamma function 

The model-agnostic nature of the NRG class should allow easy extension of the current code to other Quantum Hamiltonians.

