Welcome to the Truspace developer documenation!

This site contains all the information you need to modify or extend Truspace.

For compilation and user-guide, please refer to https://github.com/gbrandino/truspace 



## Code design and structure

TruSpace is written in C++ using and Object Oriented design pattern.

Table of content 

  * [Classes and utlities](#classes-and-utlities)
  * [Code architecture](#code-architecture)
  * [How to extened Truspace](#how-to-extened-truspace)


### Classes and utlities

The code is currently strctured to produce 3 executables and a library

- libNRG - source code in src/nrg: the core of truspace, contains the NRG routines as well as some base class to derive for specialization, namely:
  - \link NRG \endlink : The core nrg engine, model agnostic, supporting also 'Sweep'
  - \link Hamiltonian \endlink : Abstract class to handle the computation of the TSA Hamiltonian
  - \link states \endlink : Abstract class that calculates and holds all the informations
    regarding the computation basis
  - \link matrix_elements \endlink : Abstract class that calculates and holds all the information about
    matrix element of the fields on the computational basis
  - \link parameters \endlink : Class that holds the nrg parameters
  - \link filenames \endlink :  Helper class to manage output and swap files names
  - tensor2D, tensor3D, tensor4D, tensor5D, tensor6D : helper template class to handle 2D, 3D, 4D, 5D and
    6D tensors. It ensures memory contiguity and automatic deallocation
    when out of scope, without forcing compilation-time sizes.
  - sparse_tensor3D and sparse_tensor6D :helper template class to handle 3D, and
    6D sparse tensors. 

- truspace.x and truspace_MPI.x - the main executables, serial (with some multithreading support) and MPI respectively. The available classes and files are 
  - \link minimal_model_parameters \endlink : class that handles the minimal model input reading, the reading of fields and structure constants 
    and the analysis of symmetry subspaces.
  - \link worker \endlink : class that takes care of running the different parts of the
    code after all parameters are set
  - \link minimal_model_states \endlink: class that calculates and holds all the informations
    regarding on the conformal computation basis. Derived from \link states \endlink
  - chirME : class that calculates and holds all the information about
    matrix element of the fields on the conformal computational basis. Derived from \link matrix_elements \endlink
  - \link minimal_model_filenames \endlink : helper class to manage output and files names specific for the minimal model. 
    Derived from \link filenames \endlink
  - minimal_model_Hamiltonian : class that handles the computation and the matrix form
    of the TSA Hamiltonian of the perturbed minimal model. Derived from \link Hamiltonian \endlink
  - \link diag \endlink: helper class that contains Numerical Recipies diagonalization
    routines for long double (required to calculate the conformal
    computational basis)
  - \link minimal_model_config.h \endlink Configuration file holding version number, some hard-coded limits that effect array dimensions, constants and macros 
  - \link block.h \endlink Contains a routine to identify symmetry subspace by the analysis of the interaction matrix. 
  - \link characters.h \endlink Set of helper functions to calculate integers partitions and manage GiNaC computation of character expansion 
  - \link cpt.h \endlink Set of helper function to calculate structure constant of Minimal CFTs 
  - \link gamma.h \endlink Computation of gamma function 

- fields_and_structure_constants.x - executable that calculates the fields list and structure constants for arbitrary unitary minimal models. It also produces
  fields file and structure constants file that can be used in truspace.x and truspace_MPI.x


The model-agnostic nature of the NRG class should allow easy extension of the current code to other Quantum Hamiltonians.


### Code architecture

The following diagram explains how the various objects are created and called.

\image html ../img/diagram.png

The flow is the following:

  - The main, called truspace.cpp
    - creates an instance of \link minimal_model_parameters \endlink from the input file, then calls parameters::analyse_subspaces that calculates the size of Hilbert space, the structure constans and tries to find symmetry subspaces in the Hilbert space
    - creates and instace of \link worker \endlink and runs it
    - \link worker \endlink iterates through momentum subspaces and symmetry subspaces. For each of them 
  
      - creates an instance of \link minimal_model_states \endlink, which calculates (or reads from file) and hold the computational conformal basis. The basis is identified by removing null-vectors and then performing a orthgonalization
      - creates an instance if \link minimal_model_filenames \endlink, which holds the strings identifying the various output files
      - creates an instance of \link chirME \endlink, which calculates (or reads from file)  and hold the matrix elements of the pertubation(s) by calling \link chirME::calculate_or_read_chirME \endlink  Matrix element are stored on the chiral states to save memory, which means that the complete matrix element will be the product of the contribution for chiral and anti-chiral components. They are stored in \link chirME::me_chi_ortho \endlink. If the expectation value of an operator is asked in the input file, it also calculates chiral matrix elements for the specified operator. The matrix elements are stored in \link chirME::me_chi_ortho_op \endlink. 
      - iterates through values of `R`, for each of them
        - creates and instance of \link minimal_model_Hamiltonian \endlink, which holds the hamiltonian matrix and a set of routines to calculate the parts of the hamiltonian relevant for NRG/Sweep
        - creates an instance of \link NRG \endlink and runs it 
          - \link NRG \endlink performs NRG iteration and Sweeps, calculates eigenvalues for each NRG/Sweep iteration, and, if requested in the input file, eigenvectors and expectation values of an operator. For eigenvectors and expectation values, the printed value are relative to the last NRG iteration. 


### How to extened Truspace

If you want to extend truspace, you will need to create your own derived version of \link states \endlink, \link matrix_elements \endlink, \link Hamiltonian \endlink, \link parameters \endlink and \link filenames \endlink. You don't really need to respect the public interfaces for \link states \endlink, \link chirME \endlink, as these are only used by  \link Hamiltonian \endlink, that you are rewirting anyway. Instead, the interfaces of  \link Hamiltonian \endlink are important to  \link NRG \endlink, in particular 

- \link Hamiltonian::init \endlink - Initialization of the Hamiltonian, that corresponds to filling of the matrix element for the first NRG iteration
- \link Hamiltonian::NRG_block_A \endlink - Fills top-left part of the Hamiltonian matrix, usually diagonal, which contains the eigenvaues of the previous iteration
- \link Hamiltonian::NRG_block_BC \endlink - Fills the bottom-left and top-right parts of the Hamiltonian matrix, which contains the expectaion value of the perturbing operator(s) between previous step eigenstates and the newly added states of the computational basis
- \link Hamiltonian::NRG_block_D  \endlink - Fills the bottom-right part of the Hamiltonian matrix, which contains the expectaion value of the perturbing operator(s) on the newly added states of the computational basis 
- \link Hamiltonian::sweep_block_BC \endlink - Fills the bottom-left and top-right parts of the Hamiltonian matrix, which contains the expectaion value of the hamiltonian between previous step eigenstates and eigenstates of the previous sweep at the same NRG iteration 
- \link Hamiltonian::Msize \endlink - Size of the Hamiltonian matrix at each NRG iteration
- \link Hamiltonian::print_eigenValues \endlink - Prints eigenvalues to file
- \link Hamiltonian::print_eigenVectors \endlink - Prints eigenvectors to file
- \link matrix_elements::calculate_ExpVal \endlink - Calculates expextation value on the final eigenstates of a given operator. The only interface of \link matrix_elements \endlink to be strictly respected.

An example of such an application to the Quantum Harmonic oscillator in 1D is available in *src/truspace_extension_example*. In this example, we use the plain waves as computational basis and the quadratic term as perturbation. 

The driver program of the Harmonic is \link truspace_harmonic.cpp \endlink. This executable reads the two input files, the first for model-specific parameters (in this case, highest and lowest pertubed states, coupling constant and range of system sizes) and the nrg input files. Samples of thoses files for the Harmonic Oscillator are available in the folder *extension_example*.

The code then creates an instance of \link harmonic_filenames \endlink. This helper class holds the filenames used by the code. The parent class contains mandatory filenames that are used by the nrg library. In particular 

- std::string EigenName;           File name of eigenvalues file 
- std::string EigenNameFinal;      File name of the eigenvalues for last NRG iteration 
- std::string EigenVectName;       File name of eigenvectors file 
- std::string MeOpName;            File name of Expectation value on the final eigenstates

The constructor \link harmonic_filenames::harmonic_filenames \endlink should define the values of those strings. 


Using the input files, the code builds a \link custom_parameters \endlink object. This class inherits from the \link parameters \endlink class, such that the parameters from the nrg input files are stored in the inherited member variables, while the model specific are stored in the custom member variables. 
In the constructor \link custom_parameters::custom_parameters \endlink it shuold be implemented the parser of the model input file. The \link custom_parameters::print \endlink function should be implemented only if used in the main for input parameter visualization.

The driver program the continues building the computational basis, creating an \link harmonic_states \endlink object. This class inherits from \link states \endlink class. No interface should be respected, as the unpertubed states are used only by the Hamiltonian class and the matrix_elements class, that have to be custom written. We usually define a calculate method like \link harmonic_states::calculate \endlink, but those calculation could be done in the constructor as well. 

After the computational basis is created, the driver program proceeds to calculate the perturbation matrix elements. The sample class \link harmonic_matrix_elements \endlink has once again a conventional \link harmonic_matrix_elements::calculate \endlink method (once again, those calculation could be done in the constructor). As stated above, the only interface to be respected is the \link harmonic_matrix_elements::calculate_ExpVal \endlink, since it is used by the NRG in case the calculation of the expectation value of some operators are required. If you are not planning to calculate expectation values, just implement this function such that it does not do anything. \link harmonic_matrix_elements::calculate_for_ExpVal \endlink is just an helper function to store the operator matrix elements to be used later for expectation value ( in particular, in the case in which the  operator whose expectation values needs to be calculated are different than the perturbation).

The code then creates and instance of \link harmonic_Hamiltonian \endlink using the parameters, the states and the matrix element. As states above, there are a few function in this class that needs to be implemented with the correct signature. In the present example those are 

- \link harmonic_Hamiltonian::init \endlink
- \link harmonic_Hamiltonian::NRG_block_A \endlink
- \link harmonic_Hamiltonian::NRG_block_BC \endlink 
- \link harmonic_Hamiltonian::NRG_block_D  \endlink 
- \link harmonic_Hamiltonian::sweep_block_BC \endlink 
- \link harmonic_Hamiltonian::Msize \endlink 
- \link harmonic_Hamiltonian::print_eigenValues \endlink 
- \link harmonic_Hamiltonian::print_eigenVectors \endlink 

Finally, the code creates an instance of the \link NRG \endlink class and invokes the \link NRG::run \endlink method.
