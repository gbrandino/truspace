/** Class that contains the matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: October 2019 
 */

#include "matrix_elements.h"
#include "parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include "config.h"

matrix_elements::matrix_elements()
{

}

matrix_elements::~matrix_elements()
{

}

matrix_elements::matrix_elements(parameters& p, states& st)
{
    tensor3D< double > me_unperturbed_tmp(p.numPert, st.nstate, st.nstate );
    me_unperturbed = me_unperturbed_tmp;

    // allocating tensors used for matrix element calculation
    if ( p.calcExpVal )
        { 
        tensor2D< double >  me_unperturbed_for_expval_tmp(st.nstate, st.nstate );
        me_unperturbed_for_expval = me_unperturbed_for_expval_tmp;
        }
}

/*
 * This routine calculates chiral matrix elements for all the perturbing operators as well and the operator for the calculation of expectation values. 
 */

void matrix_elements::calculate(parameters& p, filenames& fnames, states& st)
{

    // here we are using just one pertubationm but the example is written such that it should be easily extended to multiple ones
    for (int i = 0; i < st.nstate; i++)
    {
        for (int j = 0; j < st.nstate; j++)
        {
            int q_i = st.state_quantum_number[i];
            int q_j = st.state_quantum_number[j];
            if (q_i == q_j){
                me_unperturbed(0,i,j) = p.m * p.pConst[0]*p.pConst[0]/24;              //the factor L^2 is added later on
            }
            else{
                int fact = abs(q_i - q_j)%2 == 0? 1 : -1;
                me_unperturbed(0,i,j) = p.m * p.pConst[0]*p.pConst[0]/((2*PI)*(2*PI))* fact/((q_i-q_j)*(q_i-q_j));   //the factor L^2 is added later on
            }
        }
    }

}

void matrix_elements::calculate_for_ExpVal(parameters& p, states& st)
{
    for (int i = 0; i < st.nstate; i++)
    {
        for (int j = 0; j < st.nstate; j++)
        {
            int q_i = st.state_quantum_number[i];
            int q_j = st.state_quantum_number[j];
            if (q_i == q_j){
                me_unperturbed_for_expval(i,j) = p.m * p.pConst[0]*p.pConst[0]/24;              //the factor L^2 is added later on
            }
            else{
                int fact = abs(q_i - q_j)%2 == 0? 1 : -1;
                me_unperturbed_for_expval(i,j) = p.m * p.pConst[0]*p.pConst[0]/((2*PI)*(2*PI))* fact/((q_i-q_j)*(q_i-q_j));   //the factor L^2 is added later on
            }
        }
    }


}


void matrix_elements::calculate_ExpVal(parameters& p, states& st, double* coeff, int size, std::string fname, double R)
{
    //calculating Expectation values on perturbed eigenstates
    std::ofstream me_op_output;
    me_op_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

    if ( p.onlyDiagExpV)
        {
        for (int i=0; i < p.maxEigen; ++i)
            {
            double temp=0;
            for (int ii=0; ii < size; ++ii)
                {
                for (int jj=0;jj < size; ++jj)
                    {
                    temp += R*R * coeff[i*size + ii] * me_unperturbed_for_expval(ii,jj) * coeff[i*size + jj];
                    }
                }
                me_op_output << std::setiosflags(std::ios::fixed) << std::setprecision(14) << temp <<" ";
            }
        me_op_output<< std::endl;
        }
    else
        {
        for (int i=0; i < p.maxEigen; ++i)
            {
            for (int j=0; j < p.maxEigen; ++j)
                {
                double temp=0;
                for (int ii=0; ii < size; ++ii)
                    {
                    for (int jj=0;jj < size; ++jj)
                        {
                        temp += R*R * coeff[i*size + ii] * me_unperturbed_for_expval(ii,jj) * coeff[j*size + jj];
                        }
                    }
                //me_op_output << std::setiosflags(std::ios::fixed) << std::setprecision(14)<<temp*pow(2.*pi/R,Del_op_tot)<<" ";
                me_op_output << std::setiosflags(std::ios::fixed) << std::setprecision(14) << temp <<" ";
                }
            me_op_output<< std::endl;
            }
        }
        me_op_output.close();
}

