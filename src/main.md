Welcome to the Truspace developer documenation!

This site contains all the information you need to modify or extend Truspace.

For compilation and user-guide, please refer to http://bitbucket.org/truspacedevelopers/truspace



## Code design and structure

TruSpace is written in C++ using and Object Oriented design pattern.

Table of content 

  * [Classes and utlities](#classes-and-utlities)
  * [Code architecture](#code-architecture)
  * [How to extened Truspace](#how-to-extened-truspace)


### Classes and utlities


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


### Code architecture

The following diagram explains how the various objects are created and called.

\image html ../img/diagram.png

The flow is the following:

- The main, called truspace.cpp
  - creates an instance of \link parameters \endlink from the input file, then calls parameters::analyse_subspaces that calculates the size of Hilbert space, the structure constans and tries to find symmetry subspaces in the Hilbert space
  - creates and instace of \link worker \endlink and runs it
- \link worker \endlink 
  - iterates through momentum subspaces and symmetry subspaces. For each of them 
  
    - creates an instance of \link states \endlink, which calculates (or reads from file) and hold the computational conformal basis. The basis is identified by removing null-vectors and then performing a orthgonalization
    - creates an instance if \link filenames \endlink, which holds the strings identifying the various output files
    - creates an instance of \link chirME \endlink, which calculates (or reads from file)  and hold the matrix elements of the pertubation(s) by calling \link chirME::calculate_or_read_chirME \endlink  Matrix element are stored on the chiral states to save memory, which means that the complete matrix element will be the product of the contribution for chiral and anti-chiral components. They are stored in \link chirME::me_chi_ortho \endlink. If the expectation value of an operator is asked in the input file, it also calculates chiral matrix elements for the specified operator. The matrix elements are stored in \link chirME::me_chi_ortho_op_l \endlink and \link chirME::me_chi_ortho_op_r \endlink for the chiral and anti-chiral component of the operator. 
    - iterates through values of `R`, for each of them
      - creates and instance of \link Hamiltonian \endlink, which hold the hamiltonian matrix and a set of routines to calculate the parts of the hamiltonian relevant for NRG/Sweep
      - creates and \link NRG \endlink object and runs it 
        - \link NRG \endlink performs NRG iteration and Sweeps, calculates eigenvalues for each NRG/Sweep iteration, and, if requested in the input file, eigenvectors and expectation values of an operator. For eigenvectors and expectation values, the printed value are relative to the last NRG iteration. 


### How to extened Truspace

If you want to extend truspace, you will need to create your own version of \link states \endlink, \link chirME \endlink and \link Hamiltonian \endlink. You don't really need to respect the public interfaces for \link states \endlink, \link chirME \endlink, as these are only used by  \link Hamiltonian \endlink, that you are rewirting anyway. Instead, the interfaces of  \link Hamiltonian \endlink are important to  \link NRG \endlink, in particular 




- \link Hamiltonian::init \endlink - Initialization of the Hamiltonian, meaning the filling of the matrix element for the first NRG iteration
- \link Hamiltonian::NRG_block_A \endlink - Fills top-left part of the Hamiltonian matrix, usually diagonal which contains the eigenvaues of the previous iteration
- \link Hamiltonian::NRG_block_BC \endlink - Fills the bottom-left and top-right parts of the Hamiltonian matrix, which contains the expectaion value of the perturbing operator(s) between previous step eigenstates and the newly added states of the computational basis
- \link Hamiltonian::NRG_block_D  \endlink - Fills the bottom-right part of the Hamiltonian matrix, which contains the expectaion value of the perturbing operator(s) on the newly added states of the computational basis 
- \link Hamiltonian::sweep_block_BC \endlink - Fills the bottom-left and top-right parts of the Hamiltonian matrix, which contains the expectaion value of the perturbing operator(s) between previous step eigenstates and eigenstates of the previous step same iteration 
- \link Hamiltonian::Msize \endlink - Size of the Hamiltonian matrix at each NRG iteration
- \link Hamiltonian::print_eigenValues \endlink - Prints eigenvalues to file
- \link Hamiltonian::print_eigenVectors \endlink - Prints eigenvectors to file
- \link chirME::calculate_ExpVal \endlink - Calculates expecation value on the final eigenstates of a given operator. Only interface of \link chirME \endlink to be strictly respected. 
