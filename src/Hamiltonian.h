/** @file Hamiltonian.h
 * \brief Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "minimal_model_parameters.h"
#include "tensors.h"
#include "chirME.h"
#include "states.h"
#include <string>

/** Hamiltonian Class.  
 * Class that calculates and holds the Hamiltonian matrix and the function to calculates the elements needed for NRG and Sweep.
 * For NRG/Sweep the structure of the matrix is 
 * 
 * |   |   |
 * |---|---| 
 * | A | B |
 * | C | D |
 * 
 * Where A is a matrix of size parameters::Nsize**2, B and C are of size  parameters::Nsize* parameters::Nstep and C is of size parameters::Nstep**2
 */
class Hamiltonian
{
    public: 
        tensor2D< double > H; /**< Tensor2D holding the Hamiltonian matrix elements */
        int Msize; /**< Size of the Hamiltonian matrix */
        double R;  /**< Current value of the system size R */
        int Momentum; /**< Current momentum subspace. */
        minimal_model_parameters& p; /**< Reference to parameters object */
        states& st; /**< Reference to state obejct */
        chirME& me; /**< Reference to chirMe object */
        int subspace; /**< Current symmetry subspace */
        /** Constructor. 
         * Allocates the correct size of the tensor H, which hold the Hamiltonian matrix elements.
         * @param[in] p Reference to a parameter object.
         * @param[in] st Reference to a state object.
         * @param[in] me Reference to a chirME object.
         * @param[in] subspace Integer descrbing the subspace currently under investigation.
         */

        Hamiltonian(minimal_model_parameters& p, states& st, chirME& me, int subspace);
        /** Destructor
         */
        ~Hamiltonian();
        /** Initializes the first iteration Hamiltonian.
         *  The first iteration hamiltonian has conformal energy on diagonal and expectation value of the perturbing operator(s) on the conformal states and stores it in the tensor H.
         *  Internally, it uses the state object member and the chirME object member
         */
        void init();

        /** Routine to print eigenvalues to file.
         * @param[in] iteration NRG iteration number.
         * @param[in] sweeep Sweep iteration.
         * @param[in] eigenvalues Array containing the eigenvalues to print.
         * @param[in] fname std::string containing the file name to print. Usually, a member variable of filenames.
         * @param[in] maxState Integer describing the maximum state considered in the current NRG iteration. Used to define effective energy of the truncation.
 
         */ 
        void print_eigenValues(int iteration, int sweep, double* eigenvalues, std::string fname, int maxState);

        /** Routine to print eigenvectors to file.
         * @param[in] eigenvectors Matrix containing the eigenvectors to print, one per row
         * @param[in] fname std::string containing the file name to print. Usually, a member variable of filenames.
         * @param[in] size Number of component of each eigenvectors (row size of the eigenvectors array(.
         */
        void print_eigenVectors(double* eigenvectors, std::string fname, int size);

        /** Routine to fill the top left corner of the NRG matrix (A block). 
         * The top left corner of the NRG/Sweep matrix usually contains the eigenvaues of the previous iteration on the diagonal, 0 otherwise
         *
         * @param[in]  eigenvalues Array containing the eigenvalues of the previous NRG iteration.
         */
        void NRG_block_A(double* eigenvalues);
        
        /** Routine to fill the bottom left and top right corners of the NRG matrix (B and C blocks).
         * The the bottom left and top right corners contain the expectaion value of the perturbing operator(s) between previous step eigenstates and the newly added states of the computational basis.
         * @param[in] iteration Index of the NRG iteration
         * @param[in] coeff Matrix containing the eigenvectors of the previous NRG iteration, one per row
         */
        void NRG_block_BC(int iteration, double* coeff);

        /** Routine to fill the bottom right corner of the NRG matrix (D block).
         *  The bottom right corner contains the expectaion value of the perturbing operator(s) on the newly added states of the computational basis
         *  @param[in] iteration Index of the NRG iteration
         */
        void NRG_block_D(int iteration);
  
        /** Routine to fill the bottom left and top right corners of the Sweep  matrix (B and C blocks).
         *  The the bottom left and top right corners contain the expectaion value of the Hamiltonian between previous step eigenstates and eigenstates of the previous step same iteration
         *  @param[in] iteration Index of the NRG iteration
         *  @param[in] coeff Matrix containing the eigenvectors of the previous NRG iteration, one per row
         *  @param[in] coeffSweep Matrix containing the eigenvectors of the previous sweep and same NRG iteration, one per row
         *  @param[in] NRGspacesize Current size of the NRG hilbert space, i.e. the number of components of eigenvectors in coeff
         *  @param[in] NRGspacesize Size of the NRG hilbert space at the previous sweep same NRG iteration, i.e. the number of components of eigenvectors in coeffSweep
         */
        void sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize);

    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */ 
        Hamiltonian();

};

#endif 
