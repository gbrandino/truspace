/** Class that contains the matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: October 2019 
 */

#include "matrix_elements.h"
#include "custom_parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include "config.h"

double Int_XsquareCosAxCosBx ( double a, double b, double x)
{

    return 1/2. * (   ((pow(a,2) * pow(x,2)  - 2*a*b*pow(x,2) + pow(b,2) * pow(x,2) - 2 ) * sin ( x * (a - b) ) ) / pow (a - b, 3) 
                    + ((pow(a,2) * pow(x,2)  + 2*a*b*pow(x,2) + pow(b,2) * pow(x,2) - 2 ) * sin ( x * (a + b) ) ) / pow (a + b, 3)   
                    + (2 * x * cos(x*(a - b))) / pow(a - b,2) +
                    + (2 * x * cos(x*(a + b))) / pow(a + b,2) );

}

double Int_XsquareSinAxSinBx ( double a, double b, double x)
{

    return 1/2. * (   ((pow(a,2) * pow(x,2)  - 2*a*b*pow(x,2) + pow(b,2) * pow(x,2) - 2 ) * sin ( x * (a - b) ) ) / pow (a - b, 3)  
                    - ((pow(a,2) * pow(x,2)  + 2*a*b*pow(x,2) + pow(b,2) * pow(x,2) - 2 ) * sin ( x * (a + b) ) ) / pow (a + b, 3) 
                    + (2 * x * cos(x*(a - b))) / pow(a - b,2)  
                    - (2 * x * cos(x*(a + b))) / pow(a + b,2) );

}


double Int_XsquareCosAxCosAx ( double a, double x)
{   

    return ( 4 * pow(a,3) * pow(x,3) + (6*pow(a,2)*pow(x,2) - 3) * sin(2*a*x) + 6*a*x*cos(2*a*x)) / (24 * pow(a,3));

}


double Int_XsquareSinAxSinAx ( double a, double x)
{

    return ( 4 * pow(a,3) * pow(x,3) - (6*pow(a,2)*pow(x,2) - 3) * sin(2*a*x) - 6*a*x*cos(2*a*x)) / (24 * pow(a,3));

}

matrix_elements::matrix_elements()
{

}

matrix_elements::~matrix_elements()
{

}

matrix_elements::matrix_elements(custom_parameters& p, states& st)
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

void matrix_elements::calculate(custom_parameters& p, filenames& fnames, states& st)
{

    // here we are using just one pertubationm but the example is written such that it should be easily extended to multiple ones
    for (int i = 0; i < st.nstate; i++)
    {
        bool i_even = false;
        if (i%2 == 0)
            i_even = true;
        for (int j = 0; j < st.nstate; j++)
        {
            bool j_even = false;
            if (j%2 == 0)
                j_even = true;
            double q_i = st.state_quantum_number[i];
            double q_j = st.state_quantum_number[j];
            if (( i_even and !j_even ) or (!i_even and j_even)){
                me_unperturbed(0,i,j) = 0.;
            }
            if ( i_even and j_even)  {
                if (fabs(q_i - q_j) < 0.000001){
                    me_unperturbed(0,i,j) = 2 *p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosAx(q_i, 1/2.0) - Int_XsquareCosAxCosAx(q_i, -1/2.0));  //the factor L^2 is added later on
                    }
                else if (fabs(q_i + q_j) < 0.000001){
                    me_unperturbed(0,i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosAx(q_i, 1/2.0) - Int_XsquareCosAxCosAx(q_i, -1/2.0));
                    }
                else{
                    me_unperturbed(0,i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosBx(q_i, q_j, 1/2.0) - Int_XsquareCosAxCosBx(q_i, q_j, -1/2.0));
                    }
            }
            if ( !i_even and !j_even)  {
                if (fabs(q_i - q_j) < 0.000001){
                    me_unperturbed(0,i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinAx(q_i, 1/2.0) - Int_XsquareSinAxSinAx(q_i, -1/2.0));  //the factor L^2 is added later on
                    }
                else if (fabs(q_i + q_j) < 0.000001){
                    me_unperturbed(0,i,j) = - 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinAx(q_i, 1/2.0) - Int_XsquareSinAxSinAx(q_i, -1/2.0));
                    }
                else{
                    me_unperturbed(0,i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinBx(q_i, q_j, 1/2.0) - Int_XsquareSinAxSinBx(q_i, q_j, -1/2.0));
                    }
            }
        }
    }
}

void matrix_elements::calculate_for_ExpVal(custom_parameters& p, states& st)
{
    for (int i = 0; i < st.nstate; i++)
    {
        bool i_even = false;
        if (i%2 == 0)
            i_even = true;
        for (int j = 0; j < st.nstate; j++)
        {
            bool j_even = false;
            if (j%2 == 0)
                j_even = true;
            double q_i = st.state_quantum_number[i];
            double q_j = st.state_quantum_number[j];
            if (( i_even and !j_even ) or (!i_even and j_even)){
                me_unperturbed(0,i,j) = 0.;
            }
            if ( i_even and j_even)  {
                if (fabs(q_i - q_j) < 0.000001){
                    me_unperturbed_for_expval(i,j) = 2 *p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosAx(q_i, 1/2.0) - Int_XsquareCosAxCosAx(q_i, -1/2.0));  //the factor L^2 is added later on
                    }
                else if (fabs(q_i + q_j) < 0.000001){ 
                    me_unperturbed_for_expval(i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosAx(q_i, 1/2.0) - Int_XsquareCosAxCosAx(q_i, -1/2.0));
                    }
                else{
                    me_unperturbed_for_expval(i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareCosAxCosBx(q_i, q_j, 1/2.0) - Int_XsquareCosAxCosBx(q_i, q_j, -1/2.0));
                    }
            }
            if ( !i_even and !j_even)  {
                if (fabs(q_i - q_j) < 0.000001){
                    me_unperturbed_for_expval(i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinAx(q_i, 1/2.0) - Int_XsquareSinAxSinAx(q_i, -1/2.0));  //the factor L^2 is added later on
                    }
                else if (fabs(q_i + q_j) < 0.000001){ 
                    me_unperturbed_for_expval(i,j) = - 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinAx(q_i, 1/2.0) - Int_XsquareSinAxSinAx(q_i, -1/2.0));
                    }
                else{
                    me_unperturbed_for_expval(i,j) = 2 * p.m * pow(p.pConst[0],2)/2 *  (Int_XsquareSinAxSinBx(q_i, q_j, 1/2.0) - Int_XsquareSinAxSinBx(q_i, q_j, -1/2.0));
                    }
            }
        }
    }

}


void matrix_elements::calculate_ExpVal(parameters* in_p, states& st, double* coeff, int size, std::string fname, double R)
{
    custom_parameters* p = (custom_parameters*)in_p;
    //calculating Expectation values on perturbed eigenstates
    std::ofstream me_op_output;
    me_op_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

    if ( p->onlyDiagExpV)
        {
        for (int i=0; i < p->maxEigen; ++i)
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
        for (int i=0; i < p->maxEigen; ++i)
            {
            for (int j=0; j < p->maxEigen; ++j)
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

