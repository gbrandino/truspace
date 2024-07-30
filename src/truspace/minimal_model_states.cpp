/*
 * Class that contains and calculates chiral states by removing null vectors from each verma module
 *
 * It does so as follows:
 *
 * a) Determines the maximum possible set of chiral states for a given
 * module by first determining the set of partitions of integers (up to
 *  the maximum level, maxLev, to be considered)
 *
 * b) Having determined the levels via the partitions, it then determines
 * which are linearly independent by computing the Kac matrix.
 *
 * More specifically, suppose that there are N possible states at a given level
 * and that these states have been ordered from 1, ..., N.  To determine
 * the set of M < N linearly independent states, one begins by forming
 * the 1x1 Kac matrix with state 1.  If this matrix is non-zero then
 * state 1 is categorized as an admissible state.  Otherwise, we move on to state 2 and repeat
 * the procedure.  If state 1 is admissible, we then form a 2x2 Kac matrix
 * involving states 1 and 2.  If we diagonalize this 2x2 matrix and there
 * are no non-zero eigenvalues, state 2 is similarly admissible, and we
 * proceed to construct the 3x3 Kac matrix involving states 1, 2 and 3.  
 * Otherwise we form a 2x2 Kac matrix with states 1 and 3.
 *
 * This procedure repeats until we have formed an MxM Kac matrix involving 
 * the maximum number possible of M linearly independent states.  After this
 * final diagonalization, we save the information about the corresponding eigenvectors 
 * as these are the orthonormalized states we will use later on.
 *  This information is saved in the set of variables whose names involve "_ortho_" while 
 * the information on the non-orthogonal but linearly independent states are stored in cstate_mode, cstate_lev, 
 * and cstate_nmode (see setup1.c for the declaration of these variables).

 * To compute the Kac matrices, we employ the routine "element_computation"
 * to compute the overlap matrix elements.
 *
 * In computing the Kac matrices, we never recompute matrix elements.  The
 * code is written in such a way that previously computed matrix elements
 * are reused.
 *
 * The program gives one the ability to save both the (non-orthogonalized) states and their overlaps
 * (set the flags SAVE_STATES and SAVE_ME to 1). 
 *
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l. , Robert Konik - BNL
 * last revision: November 2016
 */

#include "minimal_model_filenames.h"
#include "minimal_model_states.h"
#include "characters.h"
#include "minimal_model_config.h"
#include "cpt.h"
#include "element_computation.h"
#include "diag.h"
#include <time.h>
#ifdef MPI_PAR
#include <mpi.h>
#endif

minimal_model_states::minimal_model_states()
{

}


minimal_model_states::minimal_model_states(minimal_model_parameters& p)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif

    int temp=0;
    for (int j=0;j<= p.maxLev;j++) /* it was lev+1 */
        {
        temp=temp+partitions(j);
        }
    MaxChiralStates = p.PairToIndex[p.currentSubspace].size()*temp+100;

    // we cannot use initialization list to allocate tensors like cstate_ortho_modecoeff, because MaxChiralStates is calculated in this functions. 
    // Calculate it outside this function seems conceptually strange
    tensor3D< long double > cstate_ortho_modecoeff_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates, MaxChiralStates); // probably enough cstate_ortho_modecoeff_tmp(p.PairToIndex[subspace].size(),MaxChiralStates, MaxChiralStates/2)
    cstate_ortho_modecoeff = cstate_ortho_modecoeff_tmp;

    tensor2D< long double > norm_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    norm = norm_tmp;

    tensor3D< long double > cstate_mode_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates, MaxNumGrade);
    cstate_mode = cstate_mode_tmp;

    tensor2D< long double > cstate_norm_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    cstate_norm =  cstate_norm_tmp;
   
    tensor2D< int > cstate_lev_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    cstate_lev =  cstate_lev_tmp;     

    tensor2D< int > cstate_nmode_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    cstate_nmode =  cstate_nmode_tmp;

    tensor2D< int > cstate_ortho_lev_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    cstate_ortho_lev =  cstate_ortho_lev_tmp;

    tensor2D< int > cstate_ortho_nmode_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates);
    cstate_ortho_nmode =  cstate_ortho_nmode_tmp;

    tensor3D< int > cstate_ortho_mode_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates,MaxChiralStates); // probably enough cstate_ortho_mode_tmp(p.PairToIndex[subspace].size(),MaxChiralStates, MaxChiralStates/2)
    cstate_ortho_mode =  cstate_ortho_mode_tmp; 

    tensor3D< long double > me_chi_0_tmp(p.PairToIndex[p.currentSubspace].size(),MaxChiralStates, MaxChiralStates);
    me_chi_0 = me_chi_0_tmp;

    tensor2D< int > nstates_lev_tmp(p.PairToIndex[p.currentSubspace].size(), p.maxLev + 1) ;
    nstates_lev = nstates_lev_tmp;
    
    tensor2D< int > nstates_lev_each_tmp(p.PairToIndex[p.currentSubspace].size(), p.maxLev + 1) ;
    nstates_lev_each = nstates_lev_each_tmp;   

    /*determine the conformal dimesions of the primary field labelling the verma module for each verma module in the subspace*/
    h= new long double[p.PairToIndex[p.currentSubspace].size()];
    num_chir= new int[p.PairToIndex[p.currentSubspace].size()];

    std::map < std::pair < int, int>, int >::iterator it;
    int counter = 0;
    //for(int j=0;j<p.PairToIndex[subspace].size();j++)
    for(it = p.PairToIndex[p.currentSubspace].begin(); it != p.PairToIndex[p.currentSubspace].end(); it++)
        {
        int indtemp1 = it->first.first;
        int indtemp2 = it->first.second;
        h[counter]=hrs((long double)indtemp1,(long double)indtemp2,p.p,p.q);
        /* if (std::div(p.finalsubspace(subspace,j),p.maxq).rem!=0)
            {
            indtemp1=std::div(p.finalsubspace(subspace,j),p.maxq).quot;
            indtemp2=std::div(p.finalsubspace(subspace,j),p.maxq).rem-1;
            indtemp1 = p.PairToIndex(subspace).
            h[j]=hrs((long double)indtemp1+1,(long double)indtemp2+1,p.p,p.q);
            }
        else
            {
            indtemp1=std::div(p.finalsubspace(subspace,j),p.maxq).quot-1;
            indtemp2=p.maxq-1;
            h[j]=hrs((long double)indtemp1+1,(long double)indtemp2+1,p.p,p.q);
            } */
        for (int k=0;k< p.maxLev+1 ;k++)
            {
            nstates_lev(counter,k) = p.vermadimcomul(indtemp1-1,indtemp2-1,k);
            nstates_lev_each(counter,k) = p.vermadim(indtemp1-1,indtemp2-1,k);
            }
        counter++;
        }

    MaxChiralVerma=0;
    for(int j=0;j< p.PairToIndex[p.currentSubspace].size();j++)
        {
        if(nstates_lev(j,p.maxLev) > MaxChiralVerma)
            {
            MaxChiralVerma = nstates_lev(j,p.maxLev);
            }
        }


        int tempdim = 0;
        MaxStates = 0;
            for(int j=0;j< p.subdim[p.currentSubspace];j++)
                {
                std::pair< std::pair<int, int>, std::pair<int, int> > field = p.fields[p.finalsubspace(p.currentSubspace,j)-1];
                int indtemp1 = field.first.first;
                int indtemp2 = field.first.second;
                int indtemp3 = field.second.first;
                int indtemp4 = field.second.second;
                for (int k=0;k<= p.maxLev;k++)
                    {
                    for(int k1=0;k1<= p.maxLev;k1++)
                        {
                        if (abs(k1-k)==p.momentum)
                            {
                            tempdim += p.vermadim(indtemp1-1,indtemp2-1,k)*p.vermadim(indtemp3-1,indtemp4-1,k1);
                            }
                        }
                    }
                }

                /* if (std::div(p.finalsubspace(subspace,j),p.maxq).rem!=0)
                    {
                    int indtemp1 = std::div(p.finalsubspace(subspace,j), p.maxq).quot;
                    int indtemp2 = std::div(p.finalsubspace(subspace,j), p.maxq).rem-1;
                    for (int k=0;k<= p.maxLev;k++)
                        {
                        for(int k1=0;k1<= p.maxLev;k1++)
                            {
                            if (abs(k1-k)==l)
                                {
                                tempdim += p.vermadim(indtemp1,indtemp2,k)*p.vermadim(indtemp1,indtemp2,k1);
                                }
                            }
                        }
                    }
                else
                    {
                    int indtemp1 = std::div(p.finalsubspace(subspace,j), p.maxq).quot-1;
                    int indtemp2 = p.maxq-1;
                    for (int k=0;k<= p.maxLev;k++)
                        {
                        for(int k1=0;k1<= p.maxLev;k1++)
                            {
                            if (abs(k1-k)==l)
                                {
                                tempdim += p.vermadim(indtemp1,indtemp2,k)*p.vermadim(indtemp1,indtemp2,k1);
                                }
                            }
                        }
                    }
                } */
            if(tempdim>MaxStates)
                {
                MaxStates=tempdim;
                }

    state_en=new double[MaxStates];
    state_eff_en=new double[MaxStates];
    state_left=new int[MaxStates];
    state_right=new int[MaxStates];
    state_vmod_left=new int[MaxStates];
    state_vmod_right=new int[MaxStates];
    state_field=new int[MaxStates];

    for (int j=0;j<MaxStates;j++)
        {
        state_en[j]=0.;
        state_eff_en[j]=0.;
        state_left[j]=0;
        state_right[j]=0;
        state_vmod_left[j]=0;
        state_vmod_right[j]=0;
	    state_field[j]=0;
        }


}

minimal_model_states::~minimal_model_states()
{
    delete [] state_en;
    delete [] state_eff_en;
    delete [] state_left;
    delete [] state_right;
    delete [] state_vmod_left;
    delete [] state_vmod_right;
    delete [] state_field;
    delete [] h;
    delete [] num_chir;
}

void minimal_model_states::calculate_chiral(minimal_model_parameters& p, int subspace, minimal_model_filenames& fnames)
{
    // MaxPart and MaxNumGrade are (at least for now) global constant 
    // defined in config.h

#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif

    char fname[100];
    FILE *state_file;
    diag diagonalizer;
    int MaxKac = partitions(p.maxLev) + 1;
    long double* eigen=new long double[MaxKac];
    long double* work=new long double[MaxKac];

    tensor2D<long double> Kac(MaxKac,MaxKac);
    tensor2D<long double> storeKac(MaxKac,MaxKac);
    tensor2D<long double> storeKac1(MaxKac,MaxKac);

    for ( int i=0; i<MaxKac*MaxKac; ++i)
        {
        Kac.data[i]=0.;
        storeKac.data[i]=0.;
        storeKac1.data[i]=0.;
        }

    /*determine partition of levels in order to denumerate chiral states*/

    tensor3D<int> results(p.maxLev+1,MaxPart,MaxNumGrade);
    tensor2D<int> results2(MaxPart,MaxNumGrade);
    
    for(int i=0;i< MaxPart;++i) 
        {
        for(int j=0;j< MaxNumGrade;++j) 
            {
            results(0,i,j) = 0;
            results2(i,j)  = 0;
            }
        }

    int* num_part = new int[p.maxLev+1];
    tensor2D<int> num_modes(p.maxLev+1,MaxChiralStates);
    
    num_part[0] = 1;
    num_modes(0,0) = 0;
    for(int k=1;k <= p.maxLev;++k) 
        {
        for(int i=0;i< MaxPart;++i) 
            {
            for(int j=0;j< MaxNumGrade;++j) 
                {
                results(k,i,j) = 0;
                }
            }

        num_part[k] = determine_partition(k,results2);

        for (int i=0;i< MaxPart;i++)
            {
            for (int j=0;j<MaxNumGrade;j++)
                {
                results(k,i,j)=results2(i,j);
                results2(i,j)=0;
                }
            }
        for(int i=0;i< num_part[k];++i) 
            {
            for(int j=0;j< MaxNumGrade;++j) 
                {
                if (results(k,i,j) == 0) 
                    {
                    num_modes(k,i) = j;
                    break;
                    }
                }
            }
        }
    

    /*construct list of potential chiral states*/
    
    int max_num_chir = 0;
    int* pot_cstate_lev= new int[MaxChiralStates];
    int* pot_cstate_nmode= new int[MaxChiralStates];
    tensor2D<long double> pot_cstate_mode(MaxChiralStates,MaxNumGrade);
    for (int i=0; i<MaxChiralStates; ++i)
        {
        pot_cstate_lev[i] = 0;
        pot_cstate_nmode[i] = 0;
        }

    for(int i=p.minLev;i<= p.maxLev;++i) 
        {
        for(int k=0;k< num_part[i];++k) 
            {
            pot_cstate_lev[max_num_chir] = i;
            pot_cstate_nmode[max_num_chir] = num_modes(i,k);
            for(int j=0;j< num_modes(i,k);++j) 
                {
                pot_cstate_mode(max_num_chir,j) = results(i,k,j);
                }
            max_num_chir++;
            }
        }

    /*determine norm of potential chiral states*/

    int nm;
    int mop;
 
    long double* grade= new long double[MaxNumGrade];
    for ( int i=0; i<MaxNumGrade; i++)
        grade[i]=0.0;

    tensor2D<long double> pot_norm(p.IndexToPair[subspace].size(),MaxChiralStates);
    for(int i=0;i< p.IndexToPair[subspace].size();++i) 
        {
        for(int l=0;l< max_num_chir;++l) 
            {
            nm = 1 + pot_cstate_nmode[l] + pot_cstate_nmode[l];
            mop = 1 + pot_cstate_nmode[l];
            for(int j=1;j<=pot_cstate_nmode[l];++j)
                grade[j] = -1.0*pot_cstate_mode(l,j-1);
            grade[mop] = 0.0;
            for(int j=1;j<=pot_cstate_nmode[l];++j)
                grade[j+mop] = pot_cstate_mode(l,pot_cstate_nmode[l]-j);
            pot_norm(i,l) = element_computation(p.C,grade,nm,mop,h[i],h[i],1.,0.);
            }
        }


    /*winnow potential chiral states to a set of linearly independent chiral states*/
    int sokm;
    //int non_zero;
    int none_zero;
    int count;
    int count1;
    int start_currlev;
    int currlev;
    //int old_sokm;
    int j_id;
    int k_id;
    //int s;
    int index;
    //int num;
    int ortho_flag;
    int pstate;
    long double temp_ld;
    long double err_margin;
    int num_k;
    int num_l; 
    int start_l;
    int end_l; 
    int* new0=new int[p.maxLev+1];

    for(int i=0;i< p.IndexToPair[subspace].size();++i) 
        {

        ortho_flag = 0;

        num_chir[i] = 0;

        pstate = 0;
        start_currlev = 0;
        currlev = 0;
        do {

        /*add next possible mode to be checked for lin. independence*/

            if (fabsl(pot_norm(i,pstate)) > .000001) 
                {
                cstate_lev(i, num_chir[i]) = pot_cstate_lev[pstate];
                if (pot_cstate_lev[pstate] > currlev) 
                    {
                    currlev = pot_cstate_lev[pstate];
                    start_currlev = num_chir[i];
                    ortho_flag = 0;
                    }
                cstate_nmode(i, num_chir[i]) = pot_cstate_nmode[pstate];
                cstate_norm(i ,num_chir[i] ) = pot_norm(i,pstate);

                for(int j=1;j<=pot_cstate_nmode[pstate];++j)
                    cstate_mode(i, num_chir[i], j) = pot_cstate_mode(pstate,j-1);

                /*compute Kac matrix*/

                sokm = num_chir[i]+1-start_currlev;

                /*don't recompute entire Kac matrix but use previous computation*/

                Kac(sokm,sokm) = 1.000;/*pot_norm[i][pstate];*/

                for(int j=0;j<sokm;++j) 
                    {
                    j_id = j+start_currlev;

                    if (j<(sokm-1)) 
                        {
                        for(int k=0;k<sokm;++k) 
                            {
                            if (k<(sokm-1)) 
                                {
                                Kac(j+1,k+1) = storeKac(j+1,k+1);
                                }
                            else 
                                {
                                k_id = k+start_currlev;
                                nm = 1 + cstate_nmode(i,j_id) + cstate_nmode(i,k_id);
                                mop = 1 + cstate_nmode(i,j_id);
               
                                 for(int l=1;l<=cstate_nmode(i,j_id);++l)
                                    grade[l] = -1.*cstate_mode(i,j_id,l);

                                grade[mop] = 0.;
                                for(int l=1;l<=cstate_nmode(i,k_id);++l)
                                      grade[l+mop] = cstate_mode(i, k_id, cstate_nmode(i,k_id)-l+1);

                                Kac(j+1,k+1) = element_computation( p.C, grade,nm,mop,h[i],h[i],1.,0.)/sqrtl(cstate_norm(i,j_id)*cstate_norm(i,k_id));
                                Kac(k+1,j+1) = Kac(j+1,k+1);
                                }
                            }
                        }
                    }

                for(int j=0;j<sokm;++j) 
                    {
                    for(int k=0;k<sokm;++k) 
                        {
                        storeKac(j+1,k+1) = Kac(j+1,k+1);
                        storeKac1(j+1,k+1) = Kac(j+1,k+1);
                        }
                    }
                //old_sokm = sokm;

                /*change norm of final state in Kac matrix to better bring out any zero eigenvalues*/

                for(int j=0;j<sokm;++j) 
                    {
                    Kac(sokm,j+1) *= .001;
                    }

                for(int j=0;j<sokm;++j) 
                    {
                    Kac(j+1,sokm) *= .001;
                    }


                /*Determine eigenvalues of Kac matrix*/

                diagonalizer.tred2_ld(Kac,sokm,eigen,work,0);
                diagonalizer.tqli_ld(eigen,work,sokm,Kac,0);

                /*Any zeros? If not, add state to list*/

                long double min = 1.;
                for(int j=0;j<sokm;++j) 
                    {
                    if (fabsl(eigen[j+1]) < min)
                        min = eigen[j+1];
                    }

                err_margin = 0.00000000000001;

                if ((cstate_lev(i,num_chir[i]) > 10) && (cstate_lev(i, num_chir[i]) < 16))
                    {
                    err_margin = 0.0000000000001;
                    }
                else if ((cstate_lev(i,num_chir[i])) > 15 && (i==0)) 
                    {
                    err_margin = 0.000000000000001;
                    }
                else if ((cstate_lev(i,num_chir[i])) > 15 && (i==1))
                    {
                    err_margin = 0.00000000000001;
                    }

                none_zero = 1;
                for(int j=0;j<sokm;++j) 
                    {
                    temp_ld = eigen[j+1];
                    if (fabsl(temp_ld) < err_margin) /*smaller -> more likely to accept state*/
                        { 
                        none_zero = 0;
                        break;
                        }
                    }

                if (none_zero == 1) 
                    {
                    num_chir[i]++;
                    }


                /*obtain orthogonal states*/

                if ((num_chir[i] == nstates_lev(i, currlev) ) && (ortho_flag == 0)) 
                    {

                    ortho_flag = 1;
                    index = currlev>0?nstates_lev(i, currlev-1):0;

                    for(int k=0;k<sokm;++k) 
                        {
                        for(int j=0;j<sokm;++j) 
                            {
                            me_chi_0(i, j+index, k+index) = storeKac(j+1,k+1);
                            }
                        }

                    diagonalizer.tred2_ld(storeKac1,sokm,eigen,work,1);
                    diagonalizer.tqli_ld(eigen,work,sokm,storeKac1,1);

                    for(int k=start_currlev;k<num_chir[i];k++) 
                        {
                        cstate_ortho_lev(i, k) = currlev;
                        cstate_ortho_nmode(i, k) = num_chir[i]-start_currlev;  //??
                        cstate_ortho_nmode(i, k) = sokm;                       //??

                        for(int j=0;j<cstate_ortho_nmode(i, k);++j) 
                            {
                            cstate_ortho_mode(i,k,j) = j+start_currlev;
                            cstate_ortho_modecoeff(i,k,j) = storeKac1(j+1, k-start_currlev+1);
                            }

                        /*normalize orthogonal states*/

                        long double norm_ortho = 0.;
                        for(int j=0;j<cstate_ortho_nmode(i,k);++j) 
                            {
                            for(int l=0;l<cstate_ortho_nmode(i,k);++l) 
                                {
                                norm_ortho += (cstate_ortho_modecoeff(i,k,j)*cstate_ortho_modecoeff(i,k,l)*storeKac(j+1,l+1) );
                                }
                            }

                        for(int j=0;j<cstate_ortho_nmode(i,k);++j) 
                            {
                            cstate_ortho_modecoeff(i,k,j) /= sqrtl(norm_ortho);
                            }

                        }
                    }
            }
            pstate++;
        } while (pstate < max_num_chir);

#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif
    std::cout << " Verma module (" << p.IndexToPair[subspace][i].first << "," << p.IndexToPair[subspace][i].second << ") has " << num_chir[i] << " states " << std::endl;
    }


    /*record norms of chiral states*/


    for(int i=0;i< p.PairToIndex[subspace].size();++i) 
        {
        for(int l=0;l<num_chir[i];++l) 
            {
            norm(i,l) = cstate_norm(i,l);
            }
        }

    if (p.SAVE_STATES ) 
        {
#ifdef MPI_PAR
    if (myrank_mpi == 0)
        {
#endif

        for(int i=0;i< p.PairToIndex[subspace].size();++i) 
            {
            for (int k=0;k<= p.maxLev;++k) 
                {
                new0[k] = 0;
                }
    
            for(int k=0;k<num_chir[i];++k) 
                {
                num_k = (cstate_lev(i,k) > 0)?k-nstates_lev(i,cstate_lev(i,k)-1):k;
                start_l = (cstate_lev(i,k) > 0)?nstates_lev(i,cstate_lev(i,k)-1):0;
                end_l = nstates_lev(i,cstate_lev(i,k) );

                for(int l=start_l;l<end_l;++l) 
                    {
                    num_l = (cstate_lev(i,l) > 0)?l-nstates_lev(i, cstate_lev(i,l)-1):l;

                    if (cstate_lev(i,l) == cstate_lev(i,k) ) 
                        {

                        // snprintf(fname,FNLEN,"ME/ME0_%d_%d_%d_%d_%d",(int)P,(int)P1,(int)R[i],(int)S[i],(int)cstate_lev[i][k]);
                        snprintf(fname,100,"%s_v%d_%d_l%d",fnames.OverlapName.c_str(),p.IndexToPair[subspace][i].first,p.IndexToPair[subspace][i].second,cstate_lev(i,k));
                        long double en = cstate_lev(i,k) + h[i] - p.C/24.;

                        if (new0[cstate_lev(i,k)] == 0) 
                            {
                            new0[cstate_lev(i,k)] = 1;
                            state_file = fopen(fname,"w");
                            fprintf(state_file,"%d %d %3.24Lf\n",num_k,num_l,((long double)(en*me_chi_0(i,k,l) )));
                            fclose(state_file);
                            }
                        else 
                            {
                            state_file = fopen(fname,"a");
                            fprintf(state_file,"%d %d %3.24Lf\n",num_k,num_l,((long double)(en*me_chi_0(i,k,l) )));
                            fclose(state_file);
                            }           
                        }
                    }
                }
            }
#ifdef MPI_PAR
        }
#endif
        }

#ifdef MPI_PAR
    MPI_Barrier(MPI_COMM_WORLD);
#endif


    /*write states to file if so desired*/

    if (p.SAVE_STATES ) 
        {
#ifdef MPI_PAR
    if (myrank_mpi == 0)
        {
#endif
        for(int i=0;i< p.PairToIndex[subspace].size();++i) 
            {
            count = 0;
            for(int j=p.minLev ;j<= p.maxLev;++j) 
                {
                //snprintf(fname,FNLEN,"States/States_%d_%d_%d_%d_%d",j,(int)P,(int)P1,(int)R[i],(int)S[i]);
                snprintf(fname,100,"%s_v%d_%d_l%d",fnames.StateName.c_str(),p.IndexToPair[subspace][i].first, p.IndexToPair[subspace][i].second,j);
                state_file = fopen(fname,"w");
                count1 = 0;
                while ((cstate_lev(i,count) == j) && (count < num_chir[i])) 
                    {
                    fprintf(state_file,"%d %d ",count1,cstate_nmode(i, count) );
                    for(int k=1;k<=cstate_nmode(i, count);++k)
                        fprintf(state_file,"%3.1Lf ",cstate_mode(i,count,k) );
                    fprintf(state_file,"%3.24Lf\n",norm(i,count) );
                    count++;
                    count1++;
                    };
                fprintf(state_file,"-1\n");
                fclose(state_file);
                }
            }
        }
#ifdef MPI_PAR
        }
#endif

#ifdef MPI_PAR
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    delete [] new0;
    delete [] num_part;
    delete [] pot_cstate_lev;
    delete [] pot_cstate_nmode;
    delete [] grade;
    delete [] eigen;
    delete [] work; 
}


void minimal_model_states::calculate_nonchiral(minimal_model_parameters& p, int subspace, int mom)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif

  /*delineate states of non-chiral states (with or up to given  momentum) up to a given energy*/

    double nonchir_state_mom; /*new line 13/04/24 rmk*/
    int nonchir_state_mom_int; /*new line 13/04/24 rmk*/
    if(mom == -1)
        {
        nstate = 0;
            //en_cutoff = 60.;
            for(int i=0; i<p.subdim[subspace];++i) 
                {
                int field_index = p.finalsubspace(subspace,i) - 1;
                std::pair<int,int>  chiral = p.fields[field_index].first;
                int chiral_index = p.PairToIndex[subspace][chiral];
                std::pair<int,int>  antichiral = p.fields[field_index].second;  
                int antichiral_index = p.PairToIndex[subspace][antichiral];
                for(int i_l=0; i_l< num_chir[chiral_index];++i_l) 
                    {
                    for(int i_r=0; i_r< num_chir[antichiral_index];++i_r) 
                        {
                        nonchir_state_mom = ((double)cstate_lev(chiral_index,i_l)) - ((double)cstate_lev(antichiral_index,i_r)) +  h[chiral_index] - h[antichiral_index]; /*new line 13/04/24 rmk*/
                        nonchir_state_mom_int = round(nonchir_state_mom); /*edits 13/04/24 rmk*/
                        //if (abs((cstate_lev(chiral_index,i_l) - cstate_lev(antichiral_index,i_r) )) == p.momentum) /*old line*/
                        if (abs(nonchir_state_mom_int) == p.momentum) /*left abs() in place, but why is it there? this seems to give momentum both positive and negative at the same time*/
                            {
                            //if ((cstate_lev[i][i_l]+cstate_lev[i][i_r]) < en_cutoff) 
                            //    {
                                state_en[nstate] = h[chiral_index] + h[antichiral_index] + cstate_lev(chiral_index,i_l)  + cstate_lev(antichiral_index,i_r) + nstate*.00000000000001;
                                state_left[nstate] = i_l;
                                state_right[nstate] = i_r;
                                state_vmod_left[nstate] = chiral_index;
                                state_vmod_right[nstate] = antichiral_index;
                                state_field[nstate] = field_index+1;
                                nstate++;
                            //    }
                            }
                        }
                    }
                }
       
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif  
        std::cout << " Number of states of momentum "<< p.momentum << " for the current subspace  "<< nstate << std::endl;
        }
    else
        {
        //en_cutoff = 60.;
        nstate = 0; 
        for(int i=0;i< p.PairToIndex[subspace].size();++i) 
            {
            int field_index = p.finalsubspace(subspace,i) - 1;
            std::pair<int,int>  chiral = p.fields[field_index].first;
            int chiral_index = p.PairToIndex[subspace][chiral];
            std::pair<int,int>  antichiral = p.fields[field_index].second;
            int antichiral_index = p.PairToIndex[subspace][antichiral];
            for(int i_l=0;i_l< num_chir[chiral_index];++i_l) 
                { 
                for(int i_r=0;i_r< num_chir[antichiral_index];++i_r) 
                    {
                    //edits 13/04/24
                    nonchir_state_mom = ((double)cstate_lev(chiral_index,i_l)) - ((double)cstate_lev(antichiral_index,i_r)) +  h[chiral_index] - h[antichiral_index]; /*new line 13/04/24 rmk*/
                    nonchir_state_mom_int = round(nonchir_state_mom); /*edits 13/04/24 rmk*/
                    // if (abs((cstate_lev(chiral_index,i_l)-cstate_lev(antichiral_index,i_r) )) == mom) 
                    if (abs(nonchir_state_mom_int) == mom)
                        { 
                        //if ((cstate_lev(i,i_l) + cstate_lev(i,i_r) ) < en_cutoff) 
                        //    {
                            state_en[nstate] = h[chiral_index] + h[antichiral_index] + cstate_lev(chiral_index,i_l) + cstate_lev(antichiral_index,i_r) +nstate*.00000000000001;
                            state_left[nstate] = i_l;
                            state_right[nstate] = i_r;
                            state_vmod_left[nstate] = chiral_index;
                            state_vmod_right[nstate] = antichiral_index;
                            state_field[nstate] = field_index+1; //rmmk edit 21/04/24
                            nstate++;
                        //    }
                        }
                    }
                }
            }
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif
        std::cout << " Number of states with momentum "<< p.momentum << " for the current subspace  "<< nstate << std::endl;
        }


    // order states in terms of energy  
    unsigned long j,inc;
    long double v;  //v1;
    int v2_l,v2_r,v3,v4,v5;
    inc = 1;
    do 
        {
        inc *= 3;
        inc++;
        } while (inc <= (unsigned)nstate);

    do 
        {
        inc /= 3;
        for(int i=inc+1;i<=nstate;i++) 
            {
            v=state_en[i-1];
            v2_l=state_vmod_left[i-1];
            v2_r=state_vmod_right[i-1];
            v3=state_left[i-1];
            v4=state_right[i-1];
            v5=state_field[i-1];

            j=i;
            while (state_en[j-1-inc] > v) 
                {
                state_en[j-1]=state_en[j-1-inc];
                state_vmod_left[j-1]=state_vmod_left[j-1-inc];
                state_vmod_right[j-1]=state_vmod_right[j-1-inc];
                state_left[j-1]=state_left[j-1-inc];
                state_right[j-1]=state_right[j-1-inc];
                state_field[j-1]=state_field[j-1-inc];
                j -= inc;
                if (j <= inc)
                    break;
                }
            state_en[j-1]=v;
            state_vmod_left[j-1]=v2_l;
            state_vmod_right[j-1]=v2_r;
            state_left[j-1]=v3;
            state_right[j-1]=v4;
            state_field[j-1]=v5;
            } 
        } while (inc > 1);

    // determine eff_energy as a function of state number
    double h_min = h[0];
    for(int i=1;i<p.PairToIndex[subspace].size() ;++i) 
        {
        if (h[i] < h_min) 
            {
            h_min = h[i];
            }
        }


    int fpt = 0;
    int bpt = 0;
    double start_en = 0.;
    double floor_en;
    double ceiling_en;
    double fin_en = 0.;
    int diff = 0;
    for(int i=0; i< nstate; ++i) 
        {
        if (i>=fpt and i < (nstate-1)) 
            {
            fin_en = (2.*p.maxLev + 2.*h_min);
            fpt = nstate-1;
            bpt = i;
            diff = fpt-bpt;
            start_en = ((double)(floor(state_en[i])));
            floor_en = ((double)(floor(state_en[i])));
            ceiling_en = floor_en + 1.;
            for(int j=i; j< nstate;++j) 
                {
                if ((state_en[j]-ceiling_en) > .000001) 
                    {
                    fpt = j;
                    fin_en = ceiling_en;
                    diff = fpt - bpt;
                    break;
                    }
                }
            }
        if (i < (nstate-1))
            {
            state_eff_en[i] = state_en[i] + (fin_en-start_en)*(i-bpt)/diff;
            }
        else
            {
            state_eff_en[i] = state_en[i];
            }
        }


    // the following tensors are used only in sweep calculations
    if ( p.nSweep > 0 )
        {
        /* the following part is correct only for zero momentum states and diagonal partion function 
        tensor3D < int > state_left_lev_tmp(p.PairToIndex[subspace].size(),p.maxLev+1,MaxChiralVerma);
        state_left_lev = state_left_lev_tmp;     

        tensor3D < int > state_right_lev_tmp(p.PairToIndex[subspace].size(),p.maxLev+1,MaxChiralVerma);
        state_right_lev = state_right_lev_tmp;

        //en_cutoff = 60.;
        //nstate = 0;
        int shiftl;
        int shiftr;
        for(int i=0;i< p.PairToIndex[subspace].size();++i)
            {
            for(int ilev=p.minLev; ilev<= p.maxLev;ilev++)
                {
                for(int i_l=0;i_l< nstates_lev_each(i,ilev);++i_l)
                    {
                    for(int i_r=0;i_r<nstates_lev_each(i,ilev);++i_r)
                        {
                        if (ilev==0)
                            {
                            shiftl=0;
                            shiftr=0;
                            }
                        else
                            {
                            shiftl=nstates_lev(i,ilev-1);
                            shiftr=nstates_lev(i,ilev-1);
                            }

                        // to be modified for higher momentum states in sweep
                        if ((cstate_lev(i,i_l+shiftl) - cstate_lev(i,i_r+shiftr) ) == 0)
                            {
                            //if ((cstate_lev(i,i_l+shiftl) + cstate_lev(i,i_r+shiftr) ) < en_cutoff)
                            //    {   
                                state_left_lev(i,ilev,i_l)=i_l+shiftl;
                                state_right_lev(i,ilev,i_r)= i_r+shiftr;
                            //    }
                            }
                        }
                    }
                }
            }

        tensor4D< int > state_nonchiral_lev_tmp(p.PairToIndex[subspace].size(), p.maxLev + 1, MaxChiralVerma, MaxChiralVerma);
        state_nonchiral_lev = state_nonchiral_lev_tmp;

        for (int i=0;i<nstate;i++)
            {   
            int vmod=state_vmod_left[i];
            std::cout << " correct here for sweep and nondiagonal models " << std::endl;
            int lev=cstate_lev(vmod,state_left[i]);
            int shift;

            if (lev==0)
                {
                shift=0;
                }
            else
                {
                shift=nstates_lev(vmod,lev-1);
                }

            int ii=state_left[i]+1-shift-1;
            int jj=state_right[i]+1-shift-1;

            state_nonchiral_lev(vmod,lev,ii,jj)=i;
            }
            */
        } // closing if ( p.nSweep > 0 )
        
}


void minimal_model_states::read_chiral(minimal_model_parameters& p, int subspace, minimal_model_filenames& fnames)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif
    char fname[100];
    FILE *state_file;
    diag diagonalizer;
    int MaxKac = partitions(p.maxLev);
    long double* eigen=new long double[MaxKac];
    long double* work=new long double[MaxKac];
    long double* grade= new long double[MaxNumGrade];

    tensor2D<long double> Kac(MaxKac,MaxKac);
    tensor2D<long double> storeKac(MaxKac,MaxKac);
    tensor2D< int > me0_exists(p.PairToIndex[subspace].size(), p.maxLev+1);

    int count;
    int num;
    int mode_num;
    int status;
    /*read states from files*/
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif
    std::cout << " Reading states and overlap from file " << std::endl;

    for (int i=0;i< p.PairToIndex[subspace].size();++i) 
        {
        count = 0;
        for(int j= p.minLev;j <= p.maxLev;++j) 
            {
            snprintf(fname,100,"%s_v%d_%d_l%d",fnames.StateName.c_str(),p.IndexToPair[subspace][i].first, p.IndexToPair[subspace][i].second,j);
            state_file = fopen(fname,"r");
            status = fscanf(state_file,"%d",&num);
            while (num != -1) 
                {
                status = fscanf(state_file,"%d",&mode_num);
                cstate_nmode(i,count) = mode_num;
                cstate_lev(i,count) = j;
                for(int k=1;k<=mode_num;++k) 
                    {
                    status = fscanf(state_file,"%Lf",&cstate_mode(i,count,k));
                    }
                status = fscanf(state_file,"%Lf",&norm(i,count));
                count++;
                status = fscanf(state_file,"%d",&num);
                };
            num_chir[i] = count;
            fclose(state_file);
            }
        //std::cout << " Verma module " << i+1 << "  has " << num_chir[i] << " states " << std::endl;
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif
        std::cout << " Verma module " << p.finalsubspace(subspace,i) << " has " << num_chir[i] << " states " << std::endl;
        }


    /*overlap of states (multiplied by energy) */
    int num_l;
    int num_lindex; 
    long double en;
    int eq_r, eq_s;
    long double temp;
    for(int i=0;i< p.PairToIndex[subspace].size();++i) 
        {
        for(int l=0;l<= p.maxLev;++l) 
            {
            snprintf(fname,100,"%s_v%d_%d_l%d",fnames.OverlapName.c_str(),p.IndexToPair[subspace][i].first,p.IndexToPair[subspace][i].second,l);
            state_file = fopen(fname,"r");
            if (state_file == NULL) 
                {
                me0_exists(i,l) = 0;
                }
            else 
                {
                me0_exists(i,l) = 1;
                num_l = (l==0)?1:nstates_lev(i,l)-nstates_lev(i,l-1);
                num_lindex = (l==0)?0:nstates_lev(i,l-1);
                en = l + h[i] - p.C/24.;
                for(int r=0;r<num_l;++r) 
                    {
                    for(int s=0;s<num_l;++s) 
                        {
                        status = fscanf(state_file,"%d %d %Lf",&eq_r,&eq_s,&temp);
                        me_chi_0(i,num_lindex+r,num_lindex+s) = temp/en;
                        }
                    }
                fclose(state_file);
                }
            }
        }

    int nm;
    int mop;
    for(int i=0;i< p.PairToIndex[subspace].size();++i) 
        {
        for(int l=0;l<num_chir[i];++l) /*modes to left of op*/
            {
            for(int k=0;k<num_chir[i];++k)  /*modes to right of op*/
                {
                /*if state doesn't exist, compute it*/
                if (me0_exists(i,cstate_lev(i,l) ) == 0) 
                    {
                    nm = 1 + cstate_nmode(i,l) + cstate_nmode(i,k);
                    mop = 1 + cstate_nmode(i,l);
                    en = h[i]-p.C/24.;
                    for(int j=1;j<= cstate_nmode(i,l);++j) 
                        {
                        grade[j] = -1.*cstate_mode(i,l,j);
                        en += cstate_mode(i,l,j);
                        }
                    grade[mop] = 0.;
                    for(int j=1;j<= cstate_nmode(i,k);++j)
                        grade[j+mop] = cstate_mode(i,k,cstate_nmode(i,k)-j+1);
                    me_chi_0(i,l,k) = element_computation( p.C, grade,nm,mop,h[i],h[i],1.,0.)/sqrt(norm(i,l)*norm(i,k));
                    if (cstate_lev(i,l) != cstate_lev(i,k) ) 
                        {
                        if (fabs(me_chi_0(i,l,k) > .0000001)) 
                            {
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif
                            printf(" Warning: Different levels but non-zero overlap %d %d\n",l,k);
                            for(int j=1;j<=nm;++j)
                                printf(" %d %3.8Le\n",j,grade[j]);
                            }
                        }
                    }
                }
            }
        }

    /*obtain orthonormal states*/
    int index;
    for(int i=0;i< p.PairToIndex[subspace].size();++i) 
        {
        for(int j=0;j<= p.maxLev;++j) 
            {
            num = (j>0)?nstates_lev(i,j) - nstates_lev(i,j-1):nstates_lev(i,j);
            index = (j>0)?nstates_lev(i,j-1):0;

            for(int l=0;l<num;++l) 
                {
                for(int r=0;r<num;++r) 
                    {
                    storeKac(l+1,r+1) = Kac(l+1,r+1) = me_chi_0(i,l+index,r+index);
                    }
                }

            if (num > 0) 
                {

                diagonalizer.tred2_ld(Kac,num,eigen,work,1);
                diagonalizer.tqli_ld(eigen,work,num,Kac,1);


                for(int k=0;k<num;++k) 
                    {
                    cstate_ortho_lev(i,k+index) = j;
                    cstate_ortho_nmode(i,k+index) = num;

                    for(int s=0;s<num;++s) 
                        {
                        cstate_ortho_mode(i,k+index,s) = s+index;
                        cstate_ortho_modecoeff(i,k+index,s) = Kac(s+1,k+1);
                        }


                    long double norm_ortho = 0.;
                    for(int s=0;s<num;++s) 
                        {
                        for(int l=0;l<num;++l) 
                            {
                            norm_ortho += (cstate_ortho_modecoeff(i,k+index,s)*cstate_ortho_modecoeff(i,k+index,l)*storeKac(s+1,l+1));
                            }
                        }
                    for(int s=0;s<num;++s) 
                        {
                        cstate_ortho_modecoeff(i,k+index,s) /= sqrtl(norm_ortho);
                        }
                    }
                }
            }
        }


    delete [] eigen;
    delete [] work;
    delete [] grade;

}




