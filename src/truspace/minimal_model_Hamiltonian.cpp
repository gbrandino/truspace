/** minimal_model_Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "minimal_model_Hamiltonian.h"
#include "minimal_model_parameters.h"
#include "tensors.h"
#include "chirME.h"
#include "minimal_model_states.h"
#include "minimal_model_config.h"
#include <cblas.h>
#include <iostream>
#include <fstream>
#include <iomanip>

minimal_model_Hamiltonian::~minimal_model_Hamiltonian()
{


}


minimal_model_Hamiltonian::minimal_model_Hamiltonian(minimal_model_parameters& p_in, minimal_model_states& st_in, chirME& me_in)
:
Hamiltonian(p_in, st_in, me_in)
{
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    if (!p.NRG)
        {
        Msize = st_der.nstate;
        tensor2D< double> H_tmp(Msize, Msize);
        H = H_tmp;
        }
    else
        {
        Msize = p.Nsize + p.Nstep;
        tensor2D< double> H_tmp(Msize, Msize);
        H = H_tmp;
        }

    // checking printing size 
    if ( p.maxEigen == 0 )
        {
        if ( p.NRG )
            {
            p.maxEigen = p.Nsize;
            }
        else
            {
            p.maxEigen = Msize;
            }
        }

    if ( p.NRG )
        {
        if (p.maxEigen > p.Nsize)
            {
            p.maxEigen = p.Nsize;
            }
        }
    else
        {
        if (p.maxEigen > Msize)
            {
            p.maxEigen = Msize;
            }
        }


}

void minimal_model_Hamiltonian::init()
{

    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);
    for (int i=0;i<Msize*Msize;i++)
         {
         H.data[i]=0.;
         }
 
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        //effective strength of each perturbation
        double c = p_der.pConst[pert]*2.*pi*pow(R/2./pi,1. - \
                   (p_der.deltas[pert].first + \
                    p_der.deltas[pert].second)); 
        #pragma omp parallel for 
        for(int i=0; i < Msize; ++i) 
            {
            int st_vmod_left  = st_der.state_vmod_left[i];
            int st_vmod_right = st_der.state_vmod_right[i];
            int st_left = st_der.state_left[i];
            int st_right  = st_der.state_right[i];
            // Kac_vmod_{right,left}_{i,j} contains the Kac indexes of the chiral and antichiral components, as oppposed to st_vmod_{left,right} which contain the LINEAR index in the SUBSPACE
            // of the chiral and antichiral component

            std::pair <int,int> vmod_left_i  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
            std::pair <int,int> vmod_right_i = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];
            //off-diagonal elements
            #pragma omp parallel for
            for(int j=0; j < Msize; ++j) 
                {
                if( ( st_der.cstate_lev(st_der.state_vmod_left[i],st_der.state_left[i]) -      \
                      st_der.cstate_lev(st_der.state_vmod_right[i],st_der.state_right[i]) )==  \
                    ( st_der.cstate_lev(st_der.state_vmod_left[j],st_der.state_left[j]) -      \
                      st_der.cstate_lev(st_der.state_vmod_right[j],st_der.state_right[j]) ) )
                    {
                    int st_vmodj_left = st_der.state_vmod_left[j];
                    int st_vmodj_right = st_der.state_vmod_left[j];
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_right];
                    int field_i = -1;
                    int field_j = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if (p_der.fields[f] == std::make_pair(vmod_left_i,vmod_right_i) ){
                            field_i = f+1;
                        }
                        if (p_der.fields[f] == std::make_pair(vmod_left_j,vmod_right_j) ){
                            field_j = f+1;
                        }
                    }
                    if (field_i == -1){
                        std::cout << " Error, field not found! " << field_i << " (" << vmod_left_i.first << "," << vmod_left_i.second << ") (" << vmod_right_i.first << "," << vmod_right_i.second << ")"<< std::endl;
                    }
                    if (field_j == -1){
                        std::cout << " Error, field not found! " << field_j << " (" << vmod_left_j.first << "," << vmod_left_j.second << ") (" << vmod_right_j.first << "," << vmod_right_j.second <<  ")"<< std::endl;
                    }
                    long double SCtmp = p_der.cptsmall(field_i, p_der.r[pert], field_j);

                    double temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodj_left, st_left, st_der.state_left[j]).first;
                    temp *= me_der.me_chi_ortho(pert, st_vmod_right, st_vmodj_right, st_right, st_der.state_right[j]).second;
                    H(i,j) += (SCtmp * c * temp);
                    }
                }

            // diagonal elements (check 2pi factor)
            H(i,i) += 2.*pi/R*(st_der.state_en[i]-p_der.C/12.)/p_der.numPert;
    
            }
        }
    //check for hermiticity
    for(int i=0;i<Msize;++i) 
        {
        for(int j=0;j<Msize;++j) 
            {
            if (fabs(H(i,j)-H(j,i)) > .000001)
                std::cout << " Warning - hermiticity problem " << std::endl;
            }
        } 



}

void minimal_model_Hamiltonian::NRG_block_A(double* eigenvalues)
{
    // called first, clears the entire Hamiltonian, then stores the eigenvalues on the diagonal 
    // for the Nsize*Nsize submatrix
    for (int i=0;i<Msize*Msize;i++)
         {
         H.data[i]=0.;
         }
    for( int i=0; i < p.Nsize;++i)
        H(i,i) = eigenvalues[i];

}

void minimal_model_Hamiltonian::NRG_block_BC(int iteration, double* coeff)
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

    int shift = p.Nsize + p.Nstep*iteration;
    double* me_needed = new double[p.Nstep*shift];
    double* results = new double[p.Nstep*p.Nsize];
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        double c = p_der.pConst[pert]*2.*pi*pow(R/2./pi,1.- (p_der.deltas[pert].first + p_der.deltas[pert].second));
        
        for(int j=0;j < p.Nstep;++j) 
            {
            int jshift = shift + j;
            int  jshift1 = j*shift;
            for(int k=0; k < shift;++k) 
                {
                double temp;
                int st_vmod_left = st_der.state_vmod_left[jshift];
                int st_vmod_right = st_der.state_vmod_right[jshift];
                int st_left = st_der.state_left[jshift];
                int st_right = st_der.state_right[jshift];
                int st_vmodk_left = st_der.state_vmod_left[k];
                int st_vmodk_right = st_der.state_vmod_right[k];
                if( (st_der.cstate_lev(st_vmod_left,st_left) - st_der.cstate_lev(st_vmod_right, st_right)) == \
                    (st_der.cstate_lev(st_vmodk_left,st_der.state_left[k]) - st_der.cstate_lev(st_vmodk_right, st_der.state_right[k])))
                    {
                    // Kac_vmod_{right,left}_{i,j} contains the Kac indexes of the chiral and antichiral components, as oppposed to st_vmod_{left,right} which contain the LINEAR index in the SUBSPACE
                    // of the chiral and antichiral component
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];
                    std::pair <int,int> vmod_left_k  = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_left];
                    std::pair <int,int> vmod_right_k = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_right];

                    int field_j = -1;
                    int field_k = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if( p_der.fields[f] == std::make_pair(vmod_left_j, vmod_right_j) ){
                            field_j = f+1;
                        }
                        if( p_der.fields[f] == std::make_pair(vmod_left_k, vmod_right_k) ){
                            field_k = f+1;
                        }
                    }
                    if (field_j == -1 || field_k == -1){
                        std::cout << " Error, fields not found! " << field_j << " " << field_k << std::endl;
                    }

                    long double SCtmp = p_der.cptsmall(field_j, p_der.r[pert], field_k);
                    temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodk_left, st_left, st_der.state_left[k]).first;
                    temp *= SCtmp * me_der.me_chi_ortho(pert, st_vmod_right, st_vmodk_right, st_right, st_der.state_right[k]).second;
                    }
                else
                    {
                    temp=0.;
                    }
                me_needed[jshift1+k]= c*temp;
                }
            }

        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nstep,p.Nsize,shift,1.0,me_needed,shift,coeff,shift,0.0,results,p.Nsize);

        for(int i=0;i<p.Nstep;++i)
            {
            int ishift = i + p.Nsize;
            int ishift1 = i*p.Nsize;
            for(int j=0;j< p.Nsize;++j)
                {
                H(ishift, j) = results[ishift1+j];
                H(j,ishift) = results[ishift1+j];
                }
            }
        }

    delete [] results;    
    delete [] me_needed;
}

void minimal_model_Hamiltonian::NRG_block_D(int iteration)
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

   for (int pert=0; pert < p_der.numPert; pert++)
        {
        //effective strength of each perturbation
        double c = p_der.pConst[pert]*2.*pi*pow(R/2./pi,1.-(p_der.deltas[pert].first + p_der.deltas[pert].second));
        int shift = p.Nsize + p.Nstep*iteration;
        for(int i=0; i < p.Nstep; ++i)
            {
            int ishift = shift + i; 
            int st_vmod_left = st_der.state_vmod_left[ishift];
            int st_vmod_right = st_der.state_vmod_right[ishift];
            int st_left = st_der.state_left[ishift];
            int st_right  = st_der.state_right[ishift];
            std::pair <int,int> vmod_left_i  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
            std::pair <int,int> vmod_right_i = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];

            //off-diagonal elements
            for(int j=0; j < p.Nstep; ++j)
                {
                int jshift = shift + j;
                if( ( st_der.cstate_lev(st_vmod_left, st_left) - st_der.cstate_lev(st_vmod_right, st_right) )== \
                    ( st_der.cstate_lev(st_der.state_vmod_left[jshift],st_der.state_left[jshift]) - \
                    st_der.cstate_lev(st_der.state_vmod_right[jshift],st_der.state_right[jshift]) ) )
                    {
                    int st_vmodj_left = st_der.state_vmod_left[jshift];
                    int st_vmodj_right = st_der.state_vmod_right[jshift];
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_right];
                    int field_i = -1;
                    int field_j = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if (p_der.fields[f] == std::make_pair(vmod_left_i,vmod_right_i) ){
                            field_i = f+1;
                        }
                        if (p_der.fields[f] == std::make_pair(vmod_left_j,vmod_right_j) ){
                            field_j = f+1;
                        }
                    }
                    if (field_i == -1){
                        std::cout << " Error, field not found! " << field_i << " (" << vmod_left_i.first << "," << vmod_left_i.second << ") (" << vmod_right_i.first << "," << vmod_right_i.second << ")"<< std::endl;
                    }
                    if (field_j == -1){
                        std::cout << " Error, field not found! " << field_j << " (" << vmod_left_j.first << "," << vmod_left_j.second << ") (" << vmod_right_j.first << "," << vmod_right_j.second <<  ")"<< std::endl;
                    }
                    long double SCtmp = p_der.cptsmall(field_i, p_der.r[pert], field_j);

                    double temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodj_left, st_left, st_der.state_left[jshift]).first;
                    temp *= me_der.me_chi_ortho(pert, st_vmod_right, st_vmodj_right, st_right, st_der.state_right[jshift]).second;
                    H(i + p.Nsize, j + p.Nsize) += (SCtmp*c*temp);
                    }
                }

            // diagonal elements (check 2pi factor)
            H(i + p.Nsize, i + p.Nsize ) += 2.*pi/R*(st_der.state_en[shift + i]-p_der.C/12.)/p_der.numPert;

            }
        }


}
 

void minimal_model_Hamiltonian::sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize)
{

    /* 

    // The following part is correct only for zero momentum states and diagonal partition functions
    int subsizeMax = 0;
    for ( int i = 0; i<p.subdim[subspace]; ++i)
        {
        for (int j = 0; j<=p.maxLev; ++j)
            {
            if ( st.nstates_lev_each( i, j) > subsizeMax)
                subsizeMax = st.nstates_lev_each(i, j);
            }
        }

    tensor6D<double> m1(p.numPert, p.subdim[subspace], p.subdim[subspace], p.maxLev+1, p.maxLev+1, subsizeMax*subsizeMax);
    tensor6D<double> m2(p.numPert, p.subdim[subspace], p.subdim[subspace], p.maxLev+1, p.maxLev+1, subsizeMax*subsizeMax);

    //storing relevant chiral matrix elements
    for ( int pert=0; pert<p.numPert; ++pert)
        {
        for ( int vmod2=0; vmod2 < p.subdim[subspace]; ++vmod2)
            {
            for ( int lev2=p.minLev; lev2 <=p.maxLev; ++lev2)
                {
                 for ( int vmod=0; vmod < p.subdim[subspace]; ++vmod)
                    {
                    for ( int lev=p.minLev; lev <=p.maxLev; ++lev)
                        {
                        int subsize = st.nstates_lev_each(vmod,lev);
                        int subsize2 = st.nstates_lev_each(vmod2,lev2);
                        if (subsize!=0 and subsize2!=0)
                            {
                            for (int ii=0;ii<subsize2;ii++)
                                {
                                for (int jj=0;jj<subsize;jj++)
                                    {
                                    int ind1 = st.state_left_lev(vmod2,lev2,ii);
                                    int ind2 = st.state_left_lev(vmod,lev,jj);
                                    m1(pert, vmod2, vmod, lev2, lev, ii*subsize+jj)=me.me_chi_ortho(pert, vmod2, vmod, ind1, ind2).first;
                                        
                                    ind1 = st.state_right_lev(vmod2, lev2, ii);
                                    ind2 = st.state_right_lev(vmod, lev, jj);
                                    m2(pert, vmod2, vmod, lev2, lev, ii*subsize+jj) = me.me_chi_ortho(pert, vmod2, vmod, ind1, ind2).second;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
   

    double* coeffSweep_j;
    double* results;
    double* results2;
    double* H_part= new double[p.Nstep*NRGspacesize];

    for (int i=0; i<p.Nstep*NRGspacesize; ++i)
        H_part[i]=0.;        

    for (int state=0;state < p.Nstep;state++)
        {
        //double temp=0.; 
        for ( int vmod2=0; vmod2 < p.subdim[subspace]; ++vmod2)
            {
            for (int lev2 = p.minLev; lev2 <= p.maxLev; ++lev2)
                {
                int subsize2 = st.nstates_lev_each(vmod2,lev2);
                for ( int vmod=0; vmod < p.subdim[subspace]; ++vmod)
                    {
                    for(int lev = p.minLev ; lev<= p.maxLev; ++lev)
                        {
                        int subsize = st.nstates_lev_each(vmod,lev);
                        if (subsize!=0 and subsize2!=0)
                            {
                            coeffSweep_j=new double[subsize*subsize];
                            for (int ii=0;ii<subsize;ii++)
                                {
                                for (int jj=0;jj<subsize;jj++)
                                    {
                                    if (st.state_nonchiral_lev(vmod, lev, ii, jj) < NRGstoresize)
                                        {
                                        coeffSweep_j[ii*subsize+jj] = coeffSweep[ state*NRGstoresize + st.state_nonchiral_lev(vmod,lev,ii,jj) ];
                                        }
                                    else
                                        {
                                        coeffSweep_j[ii*subsize+jj]=0.;
                                        }
                                    }
                                }


                            for(int pert=0; pert < p.numPert; ++pert)
                                {
                                results=new double[subsize*subsize2];
                                results2=new double[subsize2*subsize2];

                                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,subsize2,subsize,subsize,1.0,&m1(pert,vmod2,vmod,lev2,lev,0),subsize,coeffSweep_j,subsize,0.0,results,subsize);
                                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,subsize2,subsize2,subsize,1.0,results,subsize,&m2(pert,vmod2,vmod,lev2,lev,0),subsize,0.0,results2,subsize2);

                                double c = p.pConst[pert]*2.*pi*pow(R/2./pi,1.-(p.deltas[pert].first + p.deltas[pert].second));
                                for (int ii=0;ii<subsize2;ii++)
                                    {
                                    for (int jj=0;jj<subsize2;jj++)
                                        {
                                        int ind1 = st.state_nonchiral_lev(vmod2, lev2, ii, jj);
                                        if (ind1 < NRGspacesize)
                                            {
                                            H_part[state*NRGspacesize+ind1] += c*results2[ii*subsize2+jj];
                                            }
                                        }
                                    }
                                delete [] results;
                                delete [] results2;
                                }//closing loop over perturbation
                            delete [] coeffSweep_j;
                            }//closing if
                        } //closing loop over level -->lev
                    }// closing loop over verma modules -->vmod
                } // closing loop over lev -->lev2
            } // closing loop over verma modules --> vmod2
        } //closing loop over state 
     

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGspacesize,1.0,coeff,NRGspacesize,H_part,NRGspacesize,0.0,H_part2,p.Nstep);    

    delete [] H_part;  */

    // basic version for non diagonal model / higher momentum states
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

    double* H_part = new double[p.Nsize * NRGstoresize ];   
    for (int i=0; i<p.Nsize*NRGstoresize; i++)
       H_part[i] = 0.0;

    for (int pert=0; pert < p_der.numPert; pert++)
        {
        double c = p_der.pConst[pert]*2.*pi*pow(R/2./pi,1.- (p_der.deltas[pert].first + p_der.deltas[pert].second));

        for(int i=0;i < p.Nsize;++i)
            {
            //for(int j=0;j < p.Nstep;++j)
            //    {
                for(int k=0; k < NRGspacesize;++k )
                    {
                    for (int  k1=0; k1 < NRGstoresize; ++k1)
                        {
                        double temp;
                        int st_vmod_left = st_der.state_vmod_left[k];
                        int st_vmod_right = st_der.state_vmod_right[k];
                        int st_left = st_der.state_left[k];
                        int st_right = st_der.state_right[k];
                        int st_vmodk_left = st_der.state_vmod_left[k1];
                        int st_vmodk_right = st_der.state_vmod_right[k1];
                        if( (st_der.cstate_lev(st_vmod_left,st_left) - st_der.cstate_lev(st_vmod_right, st_right)) == \
                            (st_der.cstate_lev(st_vmodk_left,st_der.state_left[k1]) - st_der.cstate_lev(st_vmodk_right, st_der.state_right[k1])))
                            {
                    // Kac_vmod_{right,left}_{i,j} contains the Kac indexes of the chiral and antichiral components, as oppposed to st_vmod_{left,right} which contain the LINEAR index in the SUBSPACE
                    // of the chiral and antichiral component
                            std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
                            std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];
                            std::pair <int,int> vmod_left_k  = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_left];
                            std::pair <int,int> vmod_right_k = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_right];

                            int field_j = -1;
                            int field_k = -1;
                            for (int f = 0; f < p_der.n_fields; f++)
                                {
                                if( p_der.fields[f] == std::make_pair(vmod_left_j, vmod_right_j) ){
                                    field_j = f+1;
                                }
                                if( p_der.fields[f] == std::make_pair(vmod_left_k, vmod_right_k) ){
                                    field_k = f+1;
                                }
                            }
                            if (field_j == -1 || field_k == -1){
                                std::cout << " Error, fields not found! " << field_j << " " << field_k << std::endl;
                            }

                            long double SCtmp = p_der.cptsmall(field_j, p_der.r[pert], field_k);
                            temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodk_left, st_left, st_der.state_left[k1]).first;
                            temp *= SCtmp * me_der.me_chi_ortho(pert, st_vmod_right, st_vmodk_right, st_right, st_der.state_right[k1]).second;
                            }
                        else
                            {
                            temp=0.;
                            }
                        H_part[i*NRGstoresize + k1] += c * coeff[i*NRGspacesize + k] *temp;
                        }
                    }
                }
            //}
        }

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,H_part,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part2,p.Nstep);

    delete [] H_part;
    // basic version till here
 
    //conformal energy contribution
    double* coefftemp= new double[p.Nsize*NRGstoresize];
    for (int i=0;i < p.Nsize;++i)
        {
        for (int k=0; k < NRGstoresize;k++)
            {
            coefftemp[i*NRGstoresize + k]=coeff[i*NRGspacesize+k]*(2.*pi/R*(st_der.state_en[k]-p_der.C/12.));
            }
        }

    double* H_part3=new double[p.Nsize*p.Nstep];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,coefftemp,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part3,p.Nstep);

    delete [] coefftemp;

    //sum up perturbation and conformal contribution
    for (int i=0; i < p.Nsize;i++)
        {
        for (int j=0; j < p.Nstep;j++)
            {
            H_part2[i*p.Nstep + j] += H_part3[i*p.Nstep+j];
            }
        }
    delete [] H_part3;

    for (int i=0; i < p.Nsize;i++)
        {
        for (int j=0; j < p.Nstep;j++)
            {
            H(i,j + p.Nsize) = H_part2[i*p.Nstep +j];
            H(j +p.Nsize,i) = H_part2[i*p.Nstep +j];
            }
        }

    delete [] H_part2;
}



void minimal_model_Hamiltonian::print_eigenValues(int iteration, int sweep, double* eigenValues, std::string fname, int maxState)
{
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    std::ofstream eigen_output;
    if ( Momentum == -1)
        {
        eigen_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
        }
    else
        {
        std::stringstream temp;
        temp << fname << "_Mom_" << Momentum;
        eigen_output.open(temp.str().c_str(), std::ofstream::out | std::ofstream::app);
        }

    eigen_output << iteration  << " " << sweep << " " << std::setiosflags(std::ios::fixed) << std::setprecision(14) << st_der.state_eff_en[maxState] << " " << R << " ";

    for (int i=0; i < p.maxEigen; ++i)
        {
        eigen_output  << eigenValues[i] << " ";
        }
    eigen_output << std::endl;
    

}

void minimal_model_Hamiltonian::print_eigenVectors(double* eigenVectors, std::string fname, int size)
{
    std::ofstream eigenVec_output;
    if ( Momentum == -1)
        {
        eigenVec_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
        }
    else
        {
        std::stringstream temp;
        temp << fname << "_Mom_" << Momentum; 
        eigenVec_output.open(temp.str().c_str(), std::ofstream::out | std::ofstream::app);
        }

    for (int i=0; i < p.maxEigen; ++i)
        {
        for (int j=0; j < size; ++j)
            {
            eigenVec_output  << eigenVectors[i*size + j]  << " ";
            }
        eigenVec_output << std::endl;
        }

}
