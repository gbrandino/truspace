/** @file minimal_model_states.h
 * \brief Class that contains and calculates chiral states by removing null vectors from each verma module
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * \date November 2017
 */

#ifndef MINIMAL_MODEL_STATES_H
#define MINIMAL_MODEL_STATES_H
#include "states.h"
#include "minimal_model_parameters.h"
#include "filenames.h"

/** Class that contains and calculates chiral states by removing null vectors from each verma module
 */
class minimal_model_states: public states
{
    public:
        /** Constructor 
         * @param[in] p Reference to a parameters object
         */
        minimal_model_states(minimal_model_parameters& p);
        /** Destructor 
         */
        ~minimal_model_states();
   
        /** Routine to calculate chiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] fname Reference to a filenames object
         */ 
        void calculate_chiral(minimal_model_parameters& p, int subspace, filenames& fname );        
        /** Routine to calculate nonchiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] mom Momentum subspace to work on
         */
        void calculate_nonchiral(minimal_model_parameters& p, int subspace, int mom);
     
        /** Routine to read chiral states from file, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] fname Reference to a filenames object
         */
        void read_chiral(minimal_model_parameters& p, int subspace, filenames& fname);
        

        int MaxChiralStates; /**< The maximum number of chiral states, calculated from character. Used for array allocation */
        int MaxChiralVerma; /**< The maximum number of chiral states in a Verma module, calculated from character. Used for array allocation */
        int MaxStates;    /**< The maximum number of nonchiral states, calculated from character. Used for array allocation */
        double* state_en; /**< Array contaning the conformal energy of every nonchiral state */
        double* state_eff_en; /**< Array contaning the effective energy of every nonchiral state */
        int* state_left;  /**< Array that, given the global index of nonchiral states, returns the index in the correct Verma module, for the bra state */
        int* state_right; /**< Array that, given the global index of nonchiral states, returns the index in the correct Verma module, for the ket state */
        //int* state_vmod;  /**< Array that, given the global index of nonchiral states, returns the index of the correct Verma module */
        int* state_vmod_left;  /**< Array that, given the global index of nonchiral states, returns the index of the correct Verma module for the chiral part*/
        int* state_vmod_right;  /**< Array that, given the global index of nonchiral states, returns the index of the correct Verma module for the antichiral part*/
        long double* h; /**< array that contains the conformal dimension of the primary field labelling the verma module for each verma module in the subspace */
        tensor3D< long double > cstate_ortho_modecoeff; /**< 3D Tensor containing the coefficient of orthogonal conformal states on the non-orthogonal basis */ 
        tensor3D< long double > me_chi_0; /**< 3D Tensor containing the overlaps of conformal states, used to remove null vectors */ 
        tensor2D< long double > norm; /**< 2D Tensor containing norms of */ 
        tensor3D< long double > cstate_mode; /**< 3D Tensor that contains the strings of lowering operators that define a non-orthogonal chiral state. The first index runs on Verma modules in the subspace, the second is the index of the chiral state. The third on the string of indeces */
        tensor2D< long double > cstate_norm; /**< Norm of conformal chiral states */
        int* num_chir;  /**< Array to store the number of chiral non-null states at each level */
        tensor2D< int > cstate_lev; /**< Array that returns the level of a chiral non-orthogonal states, given its index.  The first index runs on Verma modules in the subspace, the second is the index of the chiral state */ 
        tensor2D< int > cstate_nmode; /**< Array that returns the number of modes for a chiral non-orthogonal states , given its index.  The first index runs on Verma modules in the subspace, the second is the index of the chiral state */
        tensor2D< int > cstate_ortho_lev; /**< Array that returns the level of a chiral orthogonal states, given its index.  The first index runs on Verma modules in the subspace, the second is the index of the chiral state */ 
        tensor2D< int > cstate_ortho_nmode; /**< Array that returns the number of modes for a chiral orthogonal states, given its index.  The first index runs on Verma modules in the subspace, the second is the index of the chiral state */
        tensor3D< int > cstate_ortho_mode; /**< 3D Tensor that contains the strings of lowering operators that define an orthogonal chiral state. The first index runs on Verma modules in the subspace, the second is the index of the chiral state. The third on the string of indeces */
        int MaxKac; /**< Stores the size of the Kac determinant at each level */ 
        tensor2D< int > nstates_lev; /**< Number of chiral state up to a given level. The first index runs on Verma modules in the subspace, the second on levels  */
        tensor2D< int > nstates_lev_each; /**< Number of chiral states in each level. The first index runs on Verma modules in the subspace, the second on levels  */
        tensor3D< int > state_right_lev; /**< Used only with sweep. Tensor that given the verma module, the level and the index inside a level return the global index of the states of the chiral level in the verma module, for the anti-chiral part */
        tensor3D< int > state_left_lev; /**< Used only with sweep. Tensor that given the verma module, the level and the index inside a level return the global index of the states of the chiral level in the verma module, for the chiral part */
        tensor4D< int > state_nonchiral_lev; /**< Used only with sweep. Tensor that given the verma module, the level of a non chiral states, and the index of its chiral and antichiral components gives the index on the non chiral state in the global indexing */
    //    char overlapname[100]; 
    //    char statename[100];

    private:
        minimal_model_states();

       
};
#endif 
