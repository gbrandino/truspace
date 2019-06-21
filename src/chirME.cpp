/** Class that contains the chiral matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */

#include "chirME.h"
#include "parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"
#include "element_computation.h"
#include "cpt.h"
#include <fstream>
#include <iomanip>
#include "config.h"

chirME::chirME()
{

}

chirME::~chirME()
{

}

chirME::chirME(parameters& p, states& st, int subspace)
{
    tensor5D< std::pair<long double, long double> > me_chi_ortho_tmp(p.numPert, p.IndexToPair[subspace].size(),p.IndexToPair[subspace].size(), st.MaxChiralStates, st.MaxChiralStates );
    me_chi_ortho = me_chi_ortho_tmp;

    // allocating tensors used for matrix element calculation
    if ( p.calcExpVal )
        { 
        tensor4D< std::pair<long double, long double> >  me_chi_ortho_op_tmp(p.IndexToPair[subspace].size(),p.IndexToPair[subspace].size(), st.MaxChiralStates, st.MaxChiralStates );
        me_chi_ortho_op = me_chi_ortho_op_tmp;
        }
}

/*
 * This routine calculates chiral matrix elements for all the perturbing operators as well and the operator for the calculation of expectation values. 
 * If instructed to read matrix elements from file, it first checks is the file exists, if not it calculates them anyway
 * 
 */

void chirME::calculate_or_read_chirME(parameters& p, filenames& fnames, states& st, int subspace)
{

    long double grade[MaxNumGrade];
    long double temp; 
    long double temp1;
  
    tensor2D<int> new1(p.maxLev+1,p.maxLev+1);
    int mop;
    FILE *state_file;

    char fname[100];
       
    tensor5D< std::pair<long double, long double> > me_chi(p.numPert, p.IndexToPair[subspace].size(), p.IndexToPair[subspace].size(), st.MaxChiralStates, st.MaxChiralStates );
    tensor5D<int> me_exists(p.numPert, p.IndexToPair[subspace].size(), p.IndexToPair[subspace].size(), p.maxLev+1, p.maxLev+1);
    tensor2D<int> me0_exists(p.IndexToPair[subspace].size(), p.maxLev+1);

    //int rflag=0;
    int eq_r;
    int eq_s;
    int status;
    if (p.READ_ME)
        {
        for (int pert = 0; pert < p.numPert; ++pert)
            {
            for(int i = 0; i < p.IndexToPair[subspace].size();++i) 
                {
                for(int j = 0; j < p.IndexToPair[subspace].size(); ++j)                   
                    {
                    for(int l=0;l <= p.maxLev;++l) 
                        {
                        for(int k=0;k <= p.maxLev;++k) 
                            {
                            snprintf(fname,100,"%s_v%d_%d_v%d_%d_l%d_l%d", fnames.MeChirName[pert].c_str(), p.IndexToPair[subspace][i].first, p.IndexToPair[subspace][i].second, 
                                    p.IndexToPair[subspace][j].first, p.IndexToPair[subspace][j].second, l, k);
                            state_file = fopen(fname,"r");
                            if ((state_file == NULL)) 
                                {
                                me_exists(pert,i,j,l,k) = 0;
                                }
                            else 
                                {
                                //rflag=1;
                                me_exists(pert,i,j,l,k) = 1;
                                int num_l = (l==0)?1:st.nstates_lev(i,l)-st.nstates_lev(i,l-1);
                                int num_lindex = (l==0)?0:st.nstates_lev(i,l-1);
                                int num_k = (k==0)?1:st.nstates_lev(j,k)-st.nstates_lev(j,k-1);
                                int num_kindex = (k==0)?0:st.nstates_lev(j,k-1);
                                for(int r=0;r<num_l;++r) 
                                    { 
                                    for(int s=0;s<num_k;++s) 
                                        {
                                        status = fscanf(state_file,"%d %d %Lf %Lf",&eq_r,&eq_s,&temp, &temp1);
                                        me_chi(pert,i,j,num_lindex+r,num_kindex+s).first = temp;
                                        me_chi(pert,i,j,num_lindex+r,num_kindex+s).second = temp1;
                                        }
                                    }
                                fclose(state_file);
                                }
                            }   
                        }
                    }
                }
            }
        }
    else
        {
        for (int pert=0; pert < p.numPert; ++pert)
            {
            for(int i = 0; i < p.IndexToPair[subspace].size(); ++i) 
                {
                for(int j = 0; j < p.IndexToPair[subspace].size(); ++j) 
                    {
                    for(int l=0;l <= p.maxLev;++l) 
                        {
                        for(int k=0;k <= p.maxLev;++k) 
                            {
                            me_exists(pert,i,j,l,k) = 0;
                            }
                        } 
                    }
                }
            }
        }


    /* // copy relevant strcture constant
    tensor3D<long double> SC(p.numPert, p.IndexToPair[subspace].size(), p.IndexToPair[subspace].size());
    for (int pert=0; pert < p.numPert; ++pert)
        {
        for(int i=0;i < p.IndexToPair[subspace].size();++i)
            {
            for(int j=0;j < p.IndexToPair[subspace].size();++j) 
                {
                int R1 = p.finalsubspace(subspace,i) / p.maxq;
                int S1 = p.finalsubspace(subspace,i) % p.maxq;
                int R2 = p.finalsubspace(subspace,j) / p.maxq;
                int S2 = p.finalsubspace(subspace,j) % p.maxq;
                SC(pert,i,j)= std::sqrt(p.cpt(R1, S1, p.r[pert], p.s[pert], R2, S2));
                }
            }
        } */


    // calculate chiral matrix element. If READ_ME is specified but some files are missing, the missing ones are calculated here
    for (int pert=0; pert < p.numPert; ++pert)
        {
        for(int i=0;i < p.IndexToPair[subspace].size();++i) 
            {
            for(int j=0;j < p.IndexToPair[subspace].size();++j) 
                {
                //if (fabs(SC(pert,i,j)) > .0000001) 
                //    {
                    for(int l=0;l < st.num_chir[i];++l)  /*modes to left of op*/
                        {
                        for(int k=0;k < st.num_chir[j];++k)  /*modes to right of op*/
                            {
                            /*check if matrix element hasn't been computed, compute it*/
                            if (me_exists(pert, i, j, st.cstate_lev(i,l), st.cstate_lev(j,k) ) == 0) 
                                {
                                int nm = 1 + st.cstate_nmode(i,l) + st.cstate_nmode(j,k);
                                mop = 1 + st.cstate_nmode(i,l);
                                for(int r=1;r <= st.cstate_nmode(i,l);++r)
                                    grade[r] = -1.*st.cstate_mode(i,l,r);
                                grade[mop] = 0.;
                                for(int r=1;r <= st.cstate_nmode(j,k);++r)
                                    grade[r+mop] = st.cstate_mode(j, k, st.cstate_nmode(j,k)-r+1);
                                temp = element_computation( p.C, grade, nm, mop, st.h[i], st.h[j], 1., p.deltas[pert].first)/sqrt(st.norm(i,l)*st.norm(j,k));
                                me_chi(pert, i, j, l, k).first = temp;
                                temp = element_computation( p.C, grade, nm, mop, st.h[i], st.h[j], 1., p.deltas[pert].second)/sqrt(st.norm(i,l)*st.norm(j,k));
                                me_chi(pert, i, j, l, k).second = temp;
                                }
                            }
                        }
                //    }
                //else 
                //    {
                //    for(int l=0;l< st.num_chir[i];++l) 
                //        {
                //        for(int k=0;k<st.num_chir[j];++k) 
                //            {
                //            me_chi(pert, i, j, l, k).first  = 0.;
                //            me_chi(pert, i, j, l, k).second = 0.;
                //            }
                //        }
                //    }
                }
            }


        /*compute matrix elements of orthonormal states*/

        for(int i=0;i< p.IndexToPair[subspace].size();++i) 
            {
            for(int j=0;j < p.IndexToPair[subspace].size();++j) 
                {
                for(int k=0;k < st.num_chir[i];++k) 
                    {
                    for(int l=0;l < st.num_chir[j];++l) 
                        {
                        temp1 = 0.0;
                        for(int r=0;r < st.cstate_ortho_nmode(i,k);++r) 
                            {
                            for(int s=0;s< st.cstate_ortho_nmode(j, l);++s) 
                                {
                                temp = me_chi(pert, i, j, st.cstate_ortho_mode(i, k, r), st.cstate_ortho_mode(j,l,s) ).first;
                                temp *= (st.cstate_ortho_modecoeff(i, k, r)*st.cstate_ortho_modecoeff(j, l, s));
                                temp1 += temp;
                                }
                            }
                        me_chi_ortho(pert,i, j, k, l).first = temp1;

                        temp1 = 0.0;
                        for(int r=0;r < st.cstate_ortho_nmode(i,k);++r)
                            {
                            for(int s=0;s< st.cstate_ortho_nmode(j, l);++s)
                                {
                                temp = me_chi(pert, i, j, st.cstate_ortho_mode(i, k, r), st.cstate_ortho_mode(j,l,s) ).second;
                                temp *= (st.cstate_ortho_modecoeff(i, k, r)*st.cstate_ortho_modecoeff(j, l, s));
                                temp1 += temp;
                                }
                            }
                        me_chi_ortho(pert,i, j, k, l).second = temp1;


                        }
                    }
                }
            } 

        /*check for Hermiticity*/

        for(int i=0;i < p.IndexToPair[subspace].size();++i) 
            {
            for(int j=0;j< p.IndexToPair[subspace].size();++j) 
                {
                for(int k=0;k < st.num_chir[i];++k) 
                    {
                    for(int l=0;l < st.num_chir[j];++l) 
                        {
                        temp = me_chi(pert, i, j, k, l).first;
                        temp1 = me_chi(pert, j, i, l, k).first;
                        if (fabs(temp-temp1) > .00000001) 
                            {
                            std::cout << " Warning: Hermiticity not satisfied on non-orthogonal chiral ME" << " Vmod " << p.finalsubspace(subspace,i) << " state " << k << " Vmod " << p.finalsubspace(subspace,j) << " state " << l << " values " << temp << " " << temp1 << std::endl;
                            std::cout << st.cstate_nmode(i,k) << " " << st.cstate_nmode(j,l) << std::endl;
                            for(int r=1;r <= st.cstate_nmode(i,k);++r)
                                std::cout << r << " " << st.cstate_mode(i,k,r) << std::endl ;
                            for(int r=1;r <= st.cstate_nmode(j,l);++r)
                                std::cout << r << " " << st.cstate_mode(j,l,r) << std::endl;
                            }
                        temp = me_chi_ortho(pert, i, j, k, l).first;
                        temp1 = me_chi_ortho(pert, j, i, l, k).first;
                        if (fabs(temp-temp1) > .00001) 
                            {
                            std::cout << " Warning: Hermiticity not satisfied on orthogonal chiral ME" << " Vmod " << p.finalsubspace(subspace,i) << " state " << k << " Vmod " << p.finalsubspace(subspace,j) << " state " << l << " values " << temp << " " << temp1 << std::endl;
                            }
                        }
                    }
                }
            }

        // saving matrix elements if desired 
        if (p.SAVE_ME) 
            {
            for(int i=0;i < p.IndexToPair[subspace].size();++i) 
                {
                for(int j=0;j < p.IndexToPair[subspace].size();++j) 
                    {
                    for (int k=0; k <= p.maxLev;++k) 
                        {
                        for (int l=0; l <= p.maxLev;++l) 
                            {
                            new1(k, l) = 0;
                            }
                        }

                    //if (fabsl(SC(pert, i, j) ) > .0000000000001) 
                    //    {
                        for(int k=0;k < st.num_chir[i];++k) 
                            {
                            int num_k = ( st.cstate_lev(i, k) > 0)?k-st.nstates_lev(i, st.cstate_lev(i,k)-1):k;
                            for(int l=0;l < st.num_chir[j];++l) 
                                {
                                int num_l = (st.cstate_lev(j,l) > 0)?l-st.nstates_lev(j, st.cstate_lev(j, l)-1):l;

                                snprintf(fname,100,"%s_v%d_%d_v%d_%d_l%d_l%d",fnames.MeChirName[pert].c_str(),  p.IndexToPair[subspace][i].first, p.IndexToPair[subspace][i].second, 
                                        p.IndexToPair[subspace][j].first, p.IndexToPair[subspace][j].second, st.cstate_lev(i, k), st.cstate_lev(j,l) );
                                if (new1(st.cstate_lev(i,k), st.cstate_lev(j, l)) == 0) 
                                    {
                                    new1(st.cstate_lev(i,k), st.cstate_lev(j,l) ) = 1;
                                    state_file = fopen(fname,"w");
                                    fprintf(state_file,"%d %d %3.24Lf\n %3.24Lf\n",num_k,num_l,me_chi(pert, i, j, k,l).first, me_chi(pert, i, j, k,l).second);
                                    fclose(state_file);
                                    }
                                else 
                                    {
                                    state_file = fopen(fname,"a");
                                    fprintf(state_file,"%d %d %3.24Lf\n %3.24Lf\n",num_k,num_l,me_chi(pert, i, j, k, l).first, me_chi(pert, i, j, k, l).second);
                                    fclose(state_file);
                                    }
                                }
                            }
                    //    }
                    }
                }
            }
        } // closing loop over pert
    

}

void chirME::calculate_chirME_for_ExpVal(parameters& p, states& st, int subspace)
{
    // calculation of chiral matrix element for the operator in the expectation values
    long double grade[MaxNumGrade], temp, temp1;

    int nm,mop;

    tensor4D< std::pair<long double, long double> >  me_chi(p.IndexToPair[subspace].size(), p.IndexToPair[subspace].size(), st.MaxChiralStates, st.MaxChiralStates ); 

    //tensor2D<long double> SC(p.subdim[subspace], p.subdim[subspace]);

    // Looping through chiral and antichiral components 
    //
    // chiralirity == 0 corresponds to left component ( chiral )
    // chiralirity == 1 corresponds to right component ( anti-chiral )
    for ( int chirality = 0; chirality < 2; ++chirality) 
        {
        int R_op;
        int S_op;
        double Del_op; 

        if ( chirality == 0 )
            {
            std::pair<int, int> op = p.fields[p.r_op].first;
            Del_op = p.Kac( op.first, op.second);
            }
        else
            {
            std::pair<int, int> op = p.fields[p.r_op].second;
            Del_op = p.Kac( op.first, op.second);
            }
  
        // copy relevant strcture constant 
       /* for(int i=0;i < p.subdim[subspace];++i)
            {
            for(int j=0;j < p.subdim[subspace];++j) 
                {
                int R1 = p.finalsubspace(subspace,i) / p.maxq;
                int S1 = p.finalsubspace(subspace,i) % p.maxq;
                int R2 = p.finalsubspace(subspace,j) / p.maxq;
                int S2 = p.finalsubspace(subspace,j) % p.maxq;
                SC(i,j)= std::sqrt(p.cpt(R1, S1, p.r_op, p.s_op, R2, S2));
                }
            } */



        for(int i=0; i < p.IndexToPair[subspace].size();++i) 
            {
            for(int j=0 ; j < p.IndexToPair[subspace].size();++j) 
                {
                //if (fabs(SC(i,j) ) > .0000001) 
                //    {
                    for(int l=0; l < st.num_chir[i];++l)  /*modes to left of op*/
                        {
                        for(int k=0; k < st.num_chir[j];++k)  /*modes to right of op*/
                            {

                            nm = 1 + st.cstate_nmode(i,l) + st.cstate_nmode(j,k) ;
                            mop = 1 + st.cstate_nmode(i,l);
                            for(int r=1; r <= st.cstate_nmode(i,l);++r)
                                {
                                grade[r] = -1.*st.cstate_mode(i,l,r);
                                }
                            grade[mop] = 0.;
                            for(int r=1; r<= st.cstate_nmode(j,k);++r)
                                grade[r+mop] = st.cstate_mode(j,k, st.cstate_nmode(j,k)-r+1);
                            temp = element_computation( p.C, grade,nm,mop,st.h[i],st.h[j],1.,Del_op)/sqrt(st.norm(i,l)*st.norm(j,k));
                            if (chirality == 0)
                                me_chi(i,j,l,k).first = temp; //SC(i,j)*temp;
                            else
                                me_chi(i,j,l,k).second = temp; //SC(i,j)*temp;
                            }
                        }
                //    }
                //else 
                //    {
                //    for(int l=0; l < st.num_chir[i];++l) 
                //        {
                //        for(int k=0; k < st.num_chir[j];++k) 
                //            {
                //            me_chi(i,j,l,k).first  = 0.;
                //            me_chi(i,j,l,k).second = 0.;
                //            }
                //        }
                //    }
                }
            }

        /*compute matrix elements of orthonormal states*/

        for(int i=0; i < p.IndexToPair[subspace].size();++i) 
            {
            for(int j=0; j < p.IndexToPair[subspace].size();++j) 
                {
                for(int k=0; k < st.num_chir[i];++k) 
                    {
                    for(int l=0; l < st.num_chir[j];++l) 
                        {
                        temp1 = 0.0;
                        for(int r=0; r < st.cstate_ortho_nmode(i,k);++r) 
                            {
                            for(int s=0; s < st.cstate_ortho_nmode(j,l);++s) 
                                {
                                if (chirality == 0)
                                    temp = me_chi(i, j, st.cstate_ortho_mode(i,k,r), st.cstate_ortho_mode(j,l,s)).first;
                                else
                                    temp = me_chi(i, j, st.cstate_ortho_mode(i,k,r), st.cstate_ortho_mode(j,l,s)).second;
                                temp *= st.cstate_ortho_modecoeff(i,k,r)*st.cstate_ortho_modecoeff(j,l,s);
                                temp1 += temp;
                                }
                            }
                        if (chirality == 0)
                            me_chi_ortho_op(i,j,k,l).first = temp1;
                        else
                            me_chi_ortho_op(i,j,k,l).second = temp1;
                        }
                    }
                }
            }
        }// closing loop over chirality


}


void chirME::calculate_ExpVal(parameters& p, states& st, double* coeff, int size, int Momentum, std::string fname, double R)
{
    //calculating Expectation values on perturbed eigenstates
    
    std::ofstream me_op_output;
    if ( Momentum == -1)
        {
        me_op_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
        }
    else
        {
        std::stringstream temp;
        temp << fname << "_Mom_" << Momentum;
        me_op_output.open(temp.str().c_str(), std::ofstream::out | std::ofstream::app);
        }


    if ( p.onlyDiagExpV)
        {
        for (int i=0; i < p.maxEigen; ++i)
            {
            double temp=0;
            for (int ii=0; ii < size; ++ii)
                {
                for (int jj=0;jj < size; ++jj)
                    {
                    /* int st_vmod = st.state_vmod_left[ii];
                    int st_left = st.state_left[ii];
                    int st_right = st.state_right[ii];
                    int st_vmodj = st.state_vmod_right[jj];
                    int jshift = jj;
                    std::pair <int,int> vmod_left = p.IndexToPair[p.currnetSubspace][ii];
                    std::pair <int,int> vmod_right = p.IndexToPair[p.currentSubspace][jj];
                    int left_right;
                    for (int f = 0; f < p.n_fields; f++)
                    {
                        if( p.fields[i] == std::make_pair(v_mod_left, vmod_right) ){
                            std::cout << "oh sisisisis" << std::endl;
                            left_right = f;
                        }
                    }
                    long double SCtmp = p.cpt(left_right, r_op, right_left); */
                    int st_vmod_left = st.state_vmod_left[ii];
                    int st_vmod_right = st.state_vmod_right[ii];
                    int st_left = st.state_left[ii];
                    int st_right = st.state_right[ii];
                    int st_vmodj_left = st.state_vmod_left[jj];
                    int st_vmodj_right = st.state_vmod_right[jj];
                    int jshift = jj;
                    std::pair <int,int> vmod_left_ii  = p.IndexToPair[p.currentSubspace][ii];
                    std::pair <int,int> vmod_right_ii = p.IndexToPair[p.currentSubspace][ii];
                    std::pair <int,int> vmod_left_jj  = p.IndexToPair[p.currentSubspace][jj];
                    std::pair <int,int> vmod_right_jj = p.IndexToPair[p.currentSubspace][jj];
                    int field_ii;
                    int field_jj;
                    for (int f = 0; f < p.n_fields; f++)
                    {
                        if( p.fields[i] == std::make_pair(vmod_left_ii, vmod_right_ii) ){
                            std::cout << "oh sisisisis" << std::endl;
                            field_ii = f;
                        }
                        if( p.fields[i] == std::make_pair(vmod_left_jj, vmod_right_jj) ){
                            std::cout << "oh sisisisis" << std::endl;
                            field_jj = f;
                        }
                    }
                    long double SCtmp = p.cptsmall(field_ii, p.r_op, field_jj);

                    double temp1 = me_chi_ortho_op(st_vmod_left, st_vmodj_left, st_left, st.state_left[jshift]).first;
                    temp1 *= me_chi_ortho_op(st_vmod_right, st_vmodj_right, st_right, st.state_right[jshift]).second;
                    temp += SCtmp * coeff[i*size + ii] * temp1 * coeff[i*size + jj];
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
                        int st_vmod_left = st.state_vmod_left[ii];
                        int st_vmod_right = st.state_vmod_right[ii];
                        int st_left = st.state_left[ii];
                        int st_right = st.state_right[ii];
                        int st_vmodj_left = st.state_vmod_left[jj];
                        int st_vmodj_right = st.state_vmod_right[jj];
                        int jshift = jj;
                        std::pair <int,int> vmod_left_ii  = p.IndexToPair[p.currentSubspace][ii];
                        std::pair <int,int> vmod_right_ii = p.IndexToPair[p.currentSubspace][ii];
                        std::pair <int,int> vmod_left_jj  = p.IndexToPair[p.currentSubspace][jj];
                        std::pair <int,int> vmod_right_jj = p.IndexToPair[p.currentSubspace][jj];
                        int field_ii;
                        int field_jj;
                        for (int f = 0; f < p.n_fields; f++)
                        {
                            if( p.fields[i] == std::make_pair(vmod_left_ii, vmod_right_ii) ){
                                std::cout << "oh sisisisis" << std::endl;
                                field_ii = f;
                            }
                            if( p.fields[i] == std::make_pair(vmod_left_jj, vmod_right_jj) ){
                                std::cout << "oh sisisisis" << std::endl;
                                field_jj = f;
                            }
                        }
                        long double SCtmp = p.cptsmall(field_ii, p.r_op, field_jj);

                        double temp1 = me_chi_ortho_op(st_vmod_left, st_vmodj_left, st_left, st.state_left[jshift]).first;
                        temp1 *= me_chi_ortho_op(st_vmod_right, st_vmodj_right, st_right, st.state_right[jshift]).second;
                        temp += SCtmp * coeff[i*size + ii] * temp1 * coeff[j*size + jj];
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

