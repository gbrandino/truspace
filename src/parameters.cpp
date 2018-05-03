/** Input parameter class. 
 *  Big fuzzy object that contains 
 *    - input paramenter
 *    - Kac table 
 *    - structure constant 
 *    - verma module in each subspace
 *  \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 *  \date November 2017
 */

#include "config.h"
#include "parameters.h"
#include "tensors.h"
#include "cpt.h"
#include "characters.h"
#include "block.h"
#include <iostream>
#include <fstream>
	
parameters::parameters()
{
    

}


parameters::~parameters()
{
    delete [] r;
    delete [] s;
    delete [] deltas;
    delete [] pConst;
    if ( nModes_l_op > 0)
        delete [] mode_l_op;
    if ( nModes_r_op > 0)
        delete [] mode_r_op;
    delete [] subdim;
    delete [] leveldim; 
}


parameters::parameters(const char* input_file_name)
{
    std::ifstream input_file;
    int tmp;

    input_file.open(input_file_name);

    // reading input paramters from file. Only first elements of each line is read, the rest of the line is ignored and can be used for comments        
    input_file>>p;  input_file.ignore(256, '\n');
    input_file>>q;  input_file.ignore(256, '\n'); 
    input_file>>numPert;  input_file.ignore(256, '\n');
    input_file>>maxLev;  input_file.ignore(256, '\n');
    input_file>>minLev;  input_file.ignore(256, '\n');
    input_file>>maxMomentum;  input_file.ignore(256, '\n');
    input_file>>allMomSubspace;  input_file.ignore(256, '\n');
    input_file>>maxEigen;  input_file.ignore(256, '\n');

    r = new int[numPert];
    s = new int[numPert];
    deltas=new double[numPert];
    pConst=new double[numPert];
    for (int i=0;i<numPert;i++)
        {
        input_file>>r[i];  input_file.ignore(256, '\n');
        input_file>>s[i];  input_file.ignore(256, '\n');
        input_file>>pConst[i];  input_file.ignore(256, '\n');
        }

    input_file>>rSteps;  input_file.ignore(256, '\n');
    input_file>>rInc; input_file.ignore(256, '\n');   
    input_file>>rStart; input_file.ignore(256, '\n'); 
    input_file>>Nsize; input_file.ignore(256, '\n');  
    input_file>>Nstep;  input_file.ignore(256, '\n');


    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0) 
        TCSA=false;
    else
        TCSA=true;

    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0)
        NRG=false;
    else
        NRG=true;

    input_file>>nSweep;  input_file.ignore(256, '\n');

    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0)
        SAVE_STATES=false;
    else
        SAVE_STATES=true;

    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0)
        READ_STATES=false;
    else
        READ_STATES=true;

    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0)
        SAVE_ME=false;
    else
        SAVE_ME=true;

    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0)
        READ_ME=false;
    else
        READ_ME=true;

    input_file>>calcExpVal;  input_file.ignore(256, '\n');
    input_file>>rl_op;  input_file.ignore(256, '\n');
    input_file>>sl_op;  input_file.ignore(256, '\n');
    input_file>>nModes_l_op;  input_file.ignore(256, '\n');
    if (nModes_l_op>0)
        {
        mode_l_op= new int[nModes_l_op];
        for (int i=0;i<nModes_l_op;i++)
            {
            input_file>>mode_l_op[i];  input_file.ignore(256, '\n');
            }
        }

    input_file>>rr_op;  input_file.ignore(256, '\n');
    input_file>>sr_op;  input_file.ignore(256, '\n');
    input_file>>nModes_r_op;  input_file.ignore(256, '\n');
    if (nModes_r_op>0)
        {
        mode_r_op= new int[nModes_r_op];
        for (int i=0;i<nModes_r_op;i++)
            {
            input_file>>mode_r_op[i];  input_file.ignore(256, '\n');
            }
        }

    input_file>>onlyDiagExpV;  input_file.ignore(256, '\n');
    input_file>>saveEigenvec;  input_file.ignore(256, '\n');

    input_file.close();

    // done reading input file. Calculating relevant quantities

    // central charge
    C=1.0-6.0*std::pow(p-q,2)/double(p*q);

    // size of Kac table
    if ( std::div(q-1,2).rem == 0 )
        {
        maxq = (q-1)/2;
        maxp = p-1;
        }
    else
        {
        maxq = q-1;
        maxp = (p-1)/2;
        }
    maxl = maxp*maxq;
    
    //temporary tensor, used to properly allocate Kac now that we know the sizes
    tensor2D<double> Kac_temp(maxp,maxq);
    Kac = Kac_temp;
    
    //temporary tensor, used to properly allocate cpt now that we know the sizes
    tensor6D< complex<long double>  > cpt_temp(maxp, maxq, maxp, maxq, maxp, maxq);
    cpt = cpt_temp; 

    //temporary tensor, used to properly allocate cptsmall now that we know the sizes
    tensor3D< complex<long double>  > cptsmall_temp(maxl, maxl, maxl);
    cptsmall = cptsmall_temp;

    //temporary tensor, used to properly allocate vermadim now that we know the sizes
    tensor3D< int > vermadim_temp(maxp, maxq, maxLev+1);
    vermadim = vermadim_temp; 
    
    //temporary tensor, used to properly allocate vermadimcomul now that we know the sizes
    tensor3D< int > vermadimcomul_temp(maxp, maxq, maxLev+1);
    vermadimcomul = vermadimcomul_temp;

    //temporary tensor, used to properly allocate finalsubspace now that we know the sizes
    tensor2D< int > finalsubspace_temp(maxl,maxl);
    finalsubspace=finalsubspace_temp; 

    subdim = new int[maxl];
    leveldim = new int[maxLev+1];

    // filling Kac table with primary fields weight
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            Kac(i-1,j-1)=hrs((long double)i,(long double)j,p,q);
            }
        }


    for (int i=0;i<numPert;i++)
        {
        deltas[i]= hrs( (long double)r[i],(long double)s[i],p,q);
        }

    complex<long double> tdot;

    // Filling cpt tensord with structure constants. 
    // The tensor is indexed in the usual way with 6 indexes, that are the r,s Kac index for each field 
    // We also check each of them is not NaN or infinite. 
    // Due to the symmetry of the tensor, if one elemenent is NaN or infinite, we compute the equivalent coefficient using the corresponding permutation
    // If no permutation gives a finite number, the code prints an error and quits.
    for (int i=0; i<maxp; i++)
        {
        for (int i1=0; i1<maxq; i1++)
            {
            for (int j=0; j<maxp; j++)
                {
                for (int j1=0; j1<maxq; j1++)
	                { 
	                for (int k=0; k<maxp; k++)
	                    {
	                    for (int k1=0; k1<maxq; k1++)
		                    {               
                            tdot=dot((int)p,(int)q,i1+1,i+1,j1+1,j+1,k1+1,k+1);//123
                            if (isnan(real(tdot)) or isinf(real(tdot)) or isnan(imag(tdot)) or isinf(imag(tdot)) )
                                {
                                tdot=dot((int)p,(int)q,i1+1,i+1,k1+1,k+1,j1+1,j+1);   //132
                                if (isnan(real(tdot)) or isinf(real(tdot)) or isnan(imag(tdot)) or isinf(imag(tdot)) )
                                   {
                                   tdot=dot((int)p,(int)q,j1+1,j+1,i1+1,i+1,k1+1,k+1);  //213
                                   if (isnan(real(tdot)) or isinf(real(tdot)) or isnan(imag(tdot)) or isinf(imag(tdot)) )
                                      {
                                      tdot=dot((int)p,(int)q,j1+1,j+1,k1+1,k+1,i1+1,i+1);  //231  
                                      if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) or isinf(imag(tdot)) )
                                         {
                                         tdot=dot((int)p,(int)q,k1+1,k+1,j1+1,j+1,i1+1,i+1); //321
                                         if (isnan(real(tdot)) or isinf(real(tdot)) or isnan(imag(tdot)) or isinf(imag(tdot)) )
                                            {
                                            tdot=dot((int)p,(int)q,k1+1,k+1,i1+1,i+1,j1+1,j+1); //312 
                                            if (isnan(real(tdot)) or isinf(real(tdot)) or isnan(imag(tdot)) or isinf(imag(tdot)) ) 
                                               {
                                               std::cerr << "  ERROR, structure constants are either NaN of infinite" << std::endl;
                                               std::exit(1);
                                               }
                                            }
                                         }
                                      }
                                   }
                                }
                            cpt(i, i1, j, j1, k, k1)=tdot;
		                    }
	                    }
	                }
                }
            }
        } 


    // Filling cptsmall tensord with structure constants. 
    // The tensor is indexed using 3 indexes, so each field is labelled by a single integer, that is its position in the Kac table, counting in row_major order
    // So, a field 2,3 on ai 2 by 3 table  will have index 4.  
    // We also check each of them is not NaN or infinite. 
    // Due to the symmetry of the tensor, if one elemenent is NaN or infinite, we compute the equivalent coefficient using the corresponding permutation
    // If no permutation gives a finite number, the code prints an error and quits.
    int ai=0;	 
    for (int i=0;i<maxp;i++)
        {
        for (int i1=0;i1<maxq;i1++)
            {
            int bi=0;        
            for (int j=0;j<maxp;j++)
                {
                for (int j1=0;j1<maxq;j1++)
                	{ 
                	int ci=0;
                	for (int k=0;k<maxp;k++)
                	    {
                	    for (int k1=0;k1<maxq;k1++)
                    		{
                            tdot=dot((int)p,(int)q,i1+1,i+1,j1+1,j+1,k1+1,k+1);//123
                            if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                {
                                tdot=dot((int)p,(int)q,i1+1,i+1,k1+1,k+1,j1+1,j+1);   //132
                                if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                    {
                                    tdot=dot((int)p,(int)q,j1+1,j+1,i1+1,i+1,k1+1,k+1);  //213
                                    if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                        {
                                        tdot=dot((int)p,(int)q,j1+1,j+1,k1+1,k+1,i1+1,i+1);  //231  
                                        if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                            {
                                            tdot=dot((int)p,(int)q,k1+1,k+1,j1+1,j+1,i1+1,i+1); //321
                                            if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                                {
                                                tdot=dot((int)p,(int)q,k1+1,k+1,i1+1,i+1,j1+1,j+1); //312 
                                                if (isnan(real(tdot)) || isinf(real(tdot)) || isnan(imag(tdot)) || isinf(imag(tdot)) )
                                                    {
                                                    std::cerr <<"WARNING, structure constants are either NaN of infinite" << std::endl;
                                                    std::exit(1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                            cptsmall(ai,bi,ci)=tdot; 
                    		ci=ci+1;
		                    }
	                    }
	                bi=bi+1;
	                }
                }
            ai=ai+1;
            }
        } 
}
 

void parameters::print() const
{
    // Printing input configuration, for easy review
    std::cout << "-------------------------------" << std::endl;
    std::cout << "  TruSpace - " << VERSION        << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << " Version " << VERSION << std::endl;
    std::cout << std::endl;
    std::cout << "   Minimal model " << std::endl;
    std::cout << std::endl;
    std::cout << "   p " << p << "  q " << q << std::endl;;
    std::cout << std::endl; 
    std::cout << "   Central charge " << C << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << " Reporting input parameter" << std::endl;
    std::cout << std::endl;
    std::cout << " Number of perturbations ---> " << numPert << std::endl;
    std::cout << " Maximum conformal level considered ---> " << maxLev << std::endl;
    std::cout << " Minimum conformal level considered ---> " << minLev << std::endl;

    if (allMomSubspace)
        {
        std::cout << " Calculating spectrum up to momentum " << maxMomentum << " in separeted subspaces " << std::endl;
        }
    else
        {
        std::cout << " Calculating spectrum up to momentum " << maxMomentum << " in 1 subspace " << std::endl;
        }

    // sweep cannot handle (for now) nonzero momentum states
    if(maxMomentum>0 and NRG )
        {
        std::cout << " Warning: sweep cannot handle nonzero momentum states, setting number of sweeps to 0 ( normal NRG )" << std::endl;
        }

    if (maxEigen==0)
        { 
        std::cout << " Printing all available eigenvalues" << std::endl;
        }
    else
        {
        std::cout << " Maximum number of eigenstates in eigen file(s) ---> " << maxEigen << std::endl;
        }

    for (int i=0;i<numPert;i++)
        {
        std::cout << " r index for perturbation " << i+1 << " ---> " << r[i] << std::endl;  
        std::cout << " s index for perturbation " << i+1 << " ---> " << s[i] << std::endl;
        std::cout << " Constant for perturbation "<< i+1 << " ---> " << pConst[i] << std::endl;;
        }

    std::cout << " Number of values of R ---> " << rSteps << std::endl;
    std::cout << " Incremet at each R step (in fraction of 2*pi) ---> " << rInc << std::endl;
    std::cout << " Starting value of R (in fraction of 2*pi) ---> " << rStart << std::endl;
    std::cout << " Number of states in the basic matrix for NRG ---> " << Nsize << std::endl;
    std::cout << " Number of states added in every NRG step ---> " << Nstep << std::endl;

    if (TCSA)
       {
       std::cout << " Full calculation including TCSA " << std::endl;
       }
    else 
        {
        std::cout <<" No TCSA "<<endl;
        }

    if (NRG)
        {
        std::cout << " NRG calcutation enabled" << std::endl;
        }
    else 
        {
        std::cout<<" No NRG calculation " << std::endl;
        }

    if (NRG and nSweep>0)
        {
        std::cout << " Performing " << nSweep << " NRG sweeps" << std::endl; 
        }

    if (SAVE_STATES)
        { 
        std::cout << " Saving states to file is enabled " << std::endl; 
        } 
    else 
        {
        std::cout << " Saving states to file is disabled" << std::endl;
        } 

    if (READ_STATES)
        { 
        std::cout << " Reading states from files is enabled" << std::endl; 
        } 
    else 
        {
        std::cout << " Reading states from files is disabled" << std::endl;
        } 

    if (SAVE_ME)
        { 
        std::cout << " Saving chiral matrix elements is enabled" << std::endl; 
        } 
    else 
        {
        std::cout << " Saving chiral matrix elements is disabled" << std::endl;
        }

    if (READ_ME)
        {
        std::cout << " Reading chiral matrix elements from file is enabled" << std::endl;
        }
    else
        {
        std::cout << " Reading chiral matrix elements from file is disabled" << std::endl;
        }


    if (calcExpVal)
        {
        std::cout << " Calculation of matrix elements of ";
        if (nModes_l_op >=1)
            {
            for(int i=0;i<nModes_l_op;i++) 
                {
                std::cout<<"L_("<<mode_l_op[i]<<")";
                }
            }
        std::cout << "phi_(" << rl_op << "," << sl_op << ")";
        if (nModes_r_op >= 1)
            {
            for(int i=0;i<nModes_r_op;i++)
                {
                std::cout << " {bar}L_(" << mode_r_op[i] << ")";
                }
            }
        std::cout << "{bar}phi_(" << rr_op << "," << sr_op << ") operator is enabled." << std::endl;
        }

    if(calcExpVal)
        {
        if (onlyDiagExpV)
            {
            std::cout << " Only diagonal matrix elements (<i|phi|i>) will be calculated" << std::endl;
            }
        else
            {
            std::cout<<" Full matrix elements calculation" << std::endl;
            }
        }

    if (saveEigenvec)
        {
        std::cout << " Saving Eigenvector to file is enable" << std::endl;
        } 
    std::cout << std::endl;
    std::cout << "-------------------------" << endl;
    std::cout << std::endl;
    std::cout << " Number of primary fields " << maxl << endl;
    std::cout << std::endl;
    std::cout << "   Reduced Kac table" << std::endl;
   
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            std::cout << "  " << hrs((long double)(i),(long double)(j),p,q);
            }
        std::cout << std::endl;
        }
    std::cout << std::endl;
    std::cout << " Indexes of fields (for indications like C(abc))" << std::endl;
    std::cout << std::endl;
    int index=1;
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            std::cout << "  " << index;
            index += 1;
            }
        std::cout << std::endl;
        }

    std::cout << " r max : " << maxp << std::endl;
    std::cout << " s max : "<< maxq << std::endl;
    std::cout << std::endl;
    std::cout << "        s     " << std::endl;
    std::cout << "  |---------->" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << " r|" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << "  v   " << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "   Structure constant tensor  C(r,s)(r',s')(r'',s'')" << std::endl;
    std::cout << std::endl;

    for (int i=0;i<maxp;i++)
        {
        for (int i1=0;i1<maxq;i1++)
            {
            for (int j=0;j<maxp;j++)
                {
                for (int j1=0;j1<maxq;j1++)
	                {
	                for (int k=0;k<maxp;k++)
	                    {
	                    for (int k1=0;k1<maxq;k1++)
		                    {
		                    std::cout << "C(" << i+1 << i1+1 << ")" << "(" << j+1 << j1+1 << ")" << "(" << k+1 << k1+1 << ")=" << cpt(i,i1,j,j1,k,k1) << "     ";
		                    }
                	        std::cout << std::endl;
	                        }
	
	                    }
                    std::cout << std::endl;
                }
            }
        }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "--  Perturbation of the CFT hamiltonian  ----" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<numPert;i++)
        {
        std::cout << " phi_(" << r[i] << "," << s[i] << ")  " << std::endl;
        }
    std::cout << std::endl;

    for (int i=0;i<numPert;i++)
        {
        std::cout << " Conformal dimension of perturbation "<<i+1<<" : "<<deltas[i] << std::endl;
        }
    std::cout << std::endl;

    std::cout<<"-- Decomposition of the hamiltonian --" << std::endl;


    std::cout << " Structure constants related to the perturbation of the hamiltonian" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<numPert;i++)
        {
        std::cout << "C(a,b,[" << r[i] << "," << s[i] << "])" << "+" ;
        }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << " Sketch of the interaction matrix, each element is a verma module " << std::endl;
    
    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            std::cout << intMatrix(i,j) << " "; 
            }
        std::cout << std::endl;
        }

    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << " Number of independet states at each level for each verma module"<< std::endl;
    std::cout << std::endl;
    std::cout << " Primary fields ";
    for (int i=0;i<maxp;i++)
        {
        for (int j=0;j<maxq;j++)
            {
            std::cout << "(" << i+1 << j+1 << ") ";
            }
        }
    std::cout << std::endl;

    for (int i=0;i<=maxLev;i++)
        {
        std::cout << " level   "<< i << "   ";
        for (int j=0;j<maxp;j++) 
            {
            for (int k=0;k<maxq;k++)
                {
                std::cout << "    " << vermadim(j,k,i);
                }
            }
        std::cout << std::endl;
        }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout <<"---------------------------------------------------------------" << std::endl;
    std::cout <<" Total number of independet states up to level n for each verma module" << std::endl;
    std::cout << std::endl;
    std::cout << " Primary fields ";
    for (int i=0;i<maxp;i++)
        {
        for (int j=0;j<maxq;j++)
            {
            std::cout << "(" << i+1 << j+1 << ") ";
            }
        }
    std::cout << std::endl;

    for (int i=0;i <= maxLev;i++)
        {
        std::cout << "level   " << i << "   ";
        for (int j=0;j<maxp;j++)
            {
            for (int k=0;k<maxq;k++)
                {
                std::cout << "    " << vermadimcomul(j,k,i);
                }
            }
        std::cout << std::endl;
        }

    std::cout << std::endl << std::endl;

    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << " Total number of independet states up to level n " << std::endl;
    std::cout << std::endl;
   
    int* MaxNumStates = new int[maxLev+1]; 
    for (int i=0;i<=maxLev;i++)
        {
        int temp=0;
        for (int j=0;j<maxp;j++)
            {
            for (int k=0;k<maxq;k++)
                {
                temp += vermadimcomul(j,k,i); 
                }
            }
        MaxNumStates[i]=temp;
        }

    for (int i=0;i<= maxLev;i++)
        {
        std::cout<< " level   " << i << "    "<< MaxNumStates[i] << std::endl;
        }
    delete [] MaxNumStates;
    std::cout << std::endl << std::endl;


    // Printing subspace found
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << " Number of independent subspaces found  " << numspace << std::endl;
    std::cout << std::endl;

    for (int i=0;i<maxl;i++)
        {
        cout<<" Subspace "<<i+1<<"      ";
        for (int j=0;j<maxl;j++)
            {
            if (finalsubspace(i,j) != 0)
                {
                std::cout << finalsubspace(i,j) <<"  ";
                }
            }
        std::cout << std::endl;
        }
    std::cout << std::endl;


    for (int i=0;i<maxl;i++)
        {
        std::cout<<" Subspace "<< i+1 << "      ";
        for (int j=0;j<maxl;j++)
            {
            if (finalsubspace(i,j) != 0)
                {
                if (div(finalsubspace(i,j),maxq).rem!=0)
                    {
                    std::cout << "(" << div(finalsubspace(i,j),maxq).quot+1 << div(finalsubspace(i,j),maxq).rem << ")" << "   ";
                    }
                else
                    {
                    std::cout << "(" << div(finalsubspace(i,j),maxq).quot<<maxq << ")" << "   ";
                    }
                }
            }
        std::cout << std::endl;
        }
    std::cout<< std::endl;
    std::cout << " Subspace dimensions (number of verma modules in each subspace) " << std::endl;

    for (int i=0;i<maxl;i++)
        {
        std::cout << " " << subdim[i];
        }
    std::cout << std::endl;

    std::cout << std::endl;


    if (allMomSubspace)
        {
        for (int i=0;i<numspace;i++)
            {
            std::cout << " Subspace "<<i+1<<" should include "<< subdimStates(i,0) <<" states "<<std::endl << std::endl;
            }
        }
    else
        {
        for (int i=0;i<numspace;i++)
            {
            for(int mom=0;mom<= maxMomentum;mom++)
                {
                std::cout << " Subspace "<<i+1<<", momentum "<<mom <<" should include "<< subdimStates(i,mom) <<" states "<< std::endl;
                }
            std::cout << std::endl;
            } 
        }

}

void parameters::analyse_subspaces()
{
    tensor2D<complex <long double> > a(maxl,maxl);
    // Creating interaction matrix for subspace analysis. The matrix has size maxl by maxl, ( maxl is the number of independent verma modules).
    // For each perturbation labelled by a combined single index k ( as in cptsmall above), the matrix a(i,j)  is filled with the corresponding strcture constant
    // c_(i,j,k). This will allow to understand if the are symmetry subspaces that can be diagonalized separately. 
    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            a(i,j)=0;
            }
        }

    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            for (int k=0;k<numPert;k++)
                {
                a(i,j)=a(i,j) + cptsmall(i,j,(r[k]-1)*maxq+s[k]-1);
                }
            }
        }


    // If the expectation of a certain operator is required, we should check that such operator does not connect different symmetry subspaces. If it does, 
    // its expecation values between different subspaces can be non-zero. In such a case, we should diagonalize together all the subspaces connected by 
    // the perturbation.
    if (calcExpVal)
        {
        for (int i=0;i<maxl;i++)
            {
            for (int j=0;j<maxl;j++)
                {
                a(i,j) = a(i,j) + cptsmall(i,j,(rl_op-1)*maxq+sl_op-1)  + cptsmall(i,j,(rr_op-1)*maxq+sr_op-1);
                }
            }
        }

    // Tensor intMatrix is just the 0-1 representation of the matrix a
    tensor2D<int> intMatrix_temp(maxl,maxl);
    intMatrix = intMatrix_temp;
    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            intMatrix(i,j)=0;
            }
        }

    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            if (std::abs(a(i,j))>0.000001 )
               {
               intMatrix(i,j)=1;
               }
            }
        }


    //Calculation of independent subspaces
    blockrec(maxl,a,finalsubspace,subdim,numspace);

    tensor2D<int>  subdimStates_temp(numspace,maxMomentum+1);
    subdimStates = subdimStates_temp;

    //calculation of the number of independet states at each level in each verma module*/
    for (int i=0;i<maxp;i++)
        {
        for (int j=0;j<maxq;j++)
            {
            characexpan(i+1,j+1,p,q,maxLev,leveldim );
            vermadim(i,j,0)=1;
            vermadimcomul(i,j,0)=1;
            for (int k=1;k <= maxLev;k++)
                {
                vermadim(i,j,k) = leveldim[k];
                vermadimcomul(i,j,k) = leveldim[k] + vermadimcomul(i,j,k-1);
                }
            }
        }

    // given the states in each verma module in a given subspace, we calculate the number of non-chiral states in each subspace
    if (allMomSubspace)
        {
        for (int i=0;i<numspace;i++)
            {
            int NumVMods=subdim[i];
            int tempdim=0;
            for(int mom=0;mom<=maxMomentum;mom++)
                {
                for(int j=0;j<NumVMods;j++)
                   {
                    if (std::div(finalsubspace(i,j),maxq).rem!=0)
                        {
                        int indtemp1 = std::div(finalsubspace(i,j),maxq).quot;
                        int indtemp2 = std::div(finalsubspace(i,j),maxq).rem-1;
                        for (int k=0;k<=maxLev;k++)
                            {
                            for (int k1=0;k1<=maxLev;k1++)
                                {
                                if(abs(k-k1) == mom)
                                    {
                                    tempdim=tempdim+vermadim(indtemp1,indtemp2,k)*vermadim(indtemp1,indtemp2,k1);
                                    }
                                }          
                            }

                        }
                    else
                        {
                        int indtemp1 = std::div(finalsubspace(i,j),maxq).quot-1;
                        int indtemp2 = maxq-1;
                        for (int k=0;k<= maxLev;k++)
                            {
                            for (int k1=0;k1<= maxLev;k1++)
                                {
                                if(abs(k-k1) == mom)
                                    {
                                    tempdim=tempdim+vermadim(indtemp1,indtemp2,k)*vermadim(indtemp1,indtemp2,k1);
                                    }
                                }
                            }
                        }
                    }
                }
            subdimStates(i,0) = tempdim;
            }
        }
    else
        {
        for (int i=0;i<numspace;i++)
            {
            int NumVMods=subdim[i];
            int tempdim=0;
            for(int mom=0;mom<= maxMomentum;mom++)
                {
                tempdim=0;
                for(int j=0;j<NumVMods;j++)
                    {
                    if (std::div(finalsubspace(i,j),maxq).rem!=0)
                        {
                        int indtemp1 = std::div(finalsubspace(i,j),maxq).quot;
                        int indtemp2 = std::div(finalsubspace(i,j),maxq).rem-1;
                        for (int k=0;k<= maxLev;k++)
                            {
                            for (int k1=0;k1<=maxLev;k1++)
                                {
                                if(abs(k-k1)==mom)
                                    {
                                    tempdim=tempdim+vermadim(indtemp1,indtemp2,k)*vermadim(indtemp1,indtemp2,k1);
                                    }
                                }
                            }

                        }
                    else
                        {
                        int indtemp1 = std::div(finalsubspace(i,j),maxq).quot-1;
                        int indtemp2 = maxq-1;
                        for (int k=0;k<= maxLev;k++)
                            {
                            for (int k1=0;k1<=maxLev;k1++)
                                {
                                if(abs(k-k1)==mom)
                                    {
                                    tempdim=tempdim+vermadim(indtemp1,indtemp2,k)*vermadim(indtemp1,indtemp2,k1);
                                    }
                                }
                            }
                        }
                    }
                subdimStates(i,mom)  = tempdim;
                }
            }
        }




}
