/** Input parameter class, derived from parameters, with model-specific data  
 *  Big fuzzy object that contains 
 *    - input paramenter
 *    - Kac table 
 *    - structure constant 
 *    - verma module in each subspace
 *  \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 *  \date November 2017
 */

#include "minimal_model_config.h"
#include "minimal_model_parameters.h"
#include "tensors.h"
#include "cpt.h"
#include "characters.h"
#include "block.h"
#include "sparse_tensor.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

minimal_model_parameters::minimal_model_parameters()
{
    

}


minimal_model_parameters::~minimal_model_parameters()
{
    delete [] r;
    delete [] s;
    delete [] deltas;
    delete [] pConst;
    delete [] subdim;
    delete [] leveldim; 
}


minimal_model_parameters::minimal_model_parameters(const char* model_input_file_name, const char* nrg_input_file_name)
:
parameters::parameters(nrg_input_file_name)
{
    std::ifstream model_input_file;
    int tmp;

    model_input_file.open(model_input_file_name);

    // reading input paramters from file. Only first elements of each line is read, the rest of the line is ignored and can be used for comments        
    model_input_file>>p;  model_input_file.ignore(256, '\n');
    model_input_file>>q;  model_input_file.ignore(256, '\n'); 
    model_input_file>>numPert;  model_input_file.ignore(256, '\n');
    model_input_file>>maxLev;  model_input_file.ignore(256, '\n');
    model_input_file>>minLev;  model_input_file.ignore(256, '\n');
    model_input_file>>momentum;  model_input_file.ignore(256, '\n');
    r = new int[numPert];
    s = new int[numPert];
    deltas=new std::pair<double, double>[numPert];
    pConst=new double[numPert];
    for (int i=0;i<numPert;i++)
        {
        model_input_file>>r[i];  model_input_file.ignore(256, '\n');
        //input_file>>s[i];  input_file.ignore(256, '\n');
        model_input_file>>pConst[i];  model_input_file.ignore(256, '\n');
        }

    model_input_file>>rSteps;  model_input_file.ignore(256, '\n');
    model_input_file>>rInc; model_input_file.ignore(256, '\n');   
    model_input_file>>rStart; model_input_file.ignore(256, '\n'); 


    model_input_file>>tmp; model_input_file.ignore(256, '\n');
    if (tmp == 0) 
        TCSA=false;
    else
        TCSA=true;

    model_input_file>>tmp; model_input_file.ignore(256, '\n');
    if (tmp == 0)
        SAVE_STATES=false;
    else
        SAVE_STATES=true;

    model_input_file>>tmp; model_input_file.ignore(256, '\n');
    if (tmp == 0)
        READ_STATES=false;
    else
        READ_STATES=true;

    model_input_file>>tmp; model_input_file.ignore(256, '\n');
    if (tmp == 0)
        SAVE_ME=false;
    else
        SAVE_ME=true;

    model_input_file>>tmp; model_input_file.ignore(256, '\n');
    if (tmp == 0)
        READ_ME=false;
    else
        READ_ME=true;


    model_input_file>>r_op;  model_input_file.ignore(256, '\n');
    //input_file>>s_op;  input_file.ignore(256, '\n');

    model_input_file>>fields_file;   model_input_file.ignore(256, '\n');
    model_input_file>>constants_file;   model_input_file.ignore(256, '\n');

    model_input_file.close();

    // done reading model input file. 
    // Calculating relevant quantities

    // central charge
    C=1.0-6.0*std::pow(p-q,2)/double(p*q);

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

    std::ifstream kac_table_file(fields_file);

    std::string line;
    std::vector<std::string>  fields_strings;
    int n_lines = 0;
    // read line by line, skip line that CONTAINS '#' and split the various cells separated by "|"
    while( std::getline(kac_table_file, line) ){
//        if (line.find("#")  != std::string::npos ){
//            continue;
//        }
        char * pch;
        char * str = &line[0u];
        pch = strtok (str,"#");
        if(pch){
            fields_strings.push_back(pch);            
        }
        
/*        while (pch != NULL)
            {
            fields_strings.push_back(pch);
            
            pch = strtok (NULL, "|");
            }*/
        n_lines++;
    }

    if (n_lines == 0){
        std::cout << " Fields file empty or not found " << std::endl;
        exit(0);
    }
    n_fields = n_lines;


    for (int i=0; i < n_fields; i++){
    //    std::cout << fields_strings[i] << std::endl;
        char * pch;
        std::string tmp = fields_strings[i];
        char * str = &tmp[0u];
        pch = strtok (str,",");
        std::pair<int,int> chiral;
        std::pair<int,int> antichiral;
        chiral.first = atoi(pch);
        int count = 1;
        while (pch != NULL)
            {
            pch = strtok (NULL, ",");
            if (pch) {
                switch (count) {
                    case 1 :
                        chiral.second = atoi(pch);
                        break;
                    case 2 :
                        antichiral.first = atoi(pch);
                        break;
                    case 3 :    
                        antichiral.second = atoi(pch);
                        break;
                    }
                count++;
                }
            }
        fields.push_back(std::make_pair(chiral, antichiral));
        fields_dimensions[i+1] = std::make_pair( hrs(chiral.first, chiral.second, p, q) , hrs(antichiral.first, antichiral.second, p, q) );
    }
    kac_table_file.close();

    sparse_tensor3D<double> CptNonChiral;
    std::ifstream structure_constant_file(constants_file);

    // set structure constant to zero, such that only the one different from zero needs to be in the input file
    for(int i = 1; i <= n_fields; i++)
        for(int j = 1; j <= n_fields; j++)
            for(int k = 1; k <= n_fields; k++)
                CptNonChiral.set(i, j, k) = 0.0;

    // read line by line, skip line that CONTAINS '#' and split the various cells separated by "|"
    n_lines = 0;
    while( std::getline(structure_constant_file, line) ){
        if (line.find("#")  != std::string::npos ){
            continue;
        }
        n_lines++;
        char * pch;
        char * str = &line[0u];
        pch = strtok (str,"|");
        std::string tmp[2];
        int count = 0;
        while (pch != NULL){
            tmp[count++] = pch;
            pch = strtok (NULL, "|");
            }
        int indexes[6];
        int countIndexes = 0;
        int valTmp;
        char * pch2;
        char * str2 = &tmp[0][0u];
        pch2 = strtok (str2,",");
        while (pch2 != NULL){
            indexes[countIndexes++] = atoi(pch2);
            pch2 = strtok (NULL, ",");
            }
        CptNonChiral.set(indexes[0], indexes[1], indexes[2]) = atof(tmp[1].c_str());
    }
    if (n_lines == 0){
       std::cout << " Constants file empty or not found " << std::endl; 
       exit(0);
    }
    //temporary tensor, used to properly allocate Kac now that we know the sizes
    tensor2D<double> Kac_temp(maxp,maxq);
    Kac = Kac_temp;
 
    //temporary tensor, used to properly allocate cpt now that we know the sizes
    //tensor6D< complex<long double>  > cpt_temp(maxp, maxq, maxp, maxq, maxp, maxq);
    //cpt = cpt_temp; 


    //temporary tensor, used to properly allocate cptsmall now that we know the sizes
    //tensor3D< complex<long double>  > cptsmall_temp(maxl, maxl, maxl);
    //cptsmall = cptsmall_temp;
    cptsmall = CptNonChiral;

    //temporary tensor, used to properly allocate vermadim now that we know the sizes
    tensor3D< int > vermadim_temp(maxp, maxq, maxLev+1);
    vermadim = vermadim_temp; 
    
    //temporary tensor, used to properly allocate vermadimcomul now that we know the sizes
    tensor3D< int > vermadimcomul_temp(maxp, maxq, maxLev+1);
    vermadimcomul = vermadimcomul_temp;

    //temporary tensor, used to properly allocate finalsubspace now that we know the sizes
    tensor2D< int > finalsubspace_temp(n_fields, n_fields);
    finalsubspace=finalsubspace_temp; 

    subdim = new int[n_fields];
    leveldim = new int[maxLev+1];

    for (int i=0;i<numPert;i++)
        {
        deltas[i] = fields_dimensions[r[i]];
        }

    
    // filling Kac table with primary fields weight
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            Kac(i-1,j-1)=hrs((long double)i,(long double)j,p,q);
            }
        }

}
 

void minimal_model_parameters::print() const
{
    // Printing input configuration, for easy review
    std::cout << "-------------------------------" << std::endl;
    std::cout << "  TruSpace - Non diagonal Potts - " << VERSION        << std::endl;
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

    std::cout << " Calculating spectrum for momentum subspace " << momentum << std::endl;

    // sweep cannot handle (for now) nonzero momentum states
    if(momentum>0 and NRG )
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
        std::cout << " index for perturbation " << i+1 << " ---> " << r[i] << std::endl;  
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
        std::cout << "phi_(" << r_op << "," << s_op << ")";
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
    std::cout << " Fields filename: " << fields_file << std::endl;
    std::cout << " Constants filename: " << constants_file << std::endl;
    std::cout << std::endl;

    std::cout << "-------------------------" << endl;
    std::cout << std::endl;
    std::cout << " Number of fields " << n_fields << endl;
    std::cout << std::endl;
    std::cout << "   Reduced Kac table" << std::endl;
   
    for (int i=maxp;i>=1;i--)
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
    for (int i=0;i< n_fields;i++)
        {
        int r_tmp = fields[i].first.first;
        int s_tmp = fields[i].first.second;
        int r1_tmp = fields[i].second.first;
        int s1_tmp = fields[i].second.second;
        std::cout << i+1 << "| ( " << r_tmp << ", " << s_tmp << " ) ( " << r1_tmp << ", " << s1_tmp << " ) | ( " <<  fields_dimensions.at(i+1).first << " , " << fields_dimensions.at(i+1).second << " ) " << std::endl;
        }
    std::cout << std::endl;
    std::cout << " r max : " << maxp << std::endl;
    std::cout << " s max : "<< maxq << std::endl;
    std::cout << std::endl;
/*    std::cout << "        s     " << std::endl;
    std::cout << "  |---------->" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << " r|" << std::endl;
    std::cout << "  |" << std::endl;
    std::cout << "  v   " << std::endl;
    std::cout << std::endl;*/
    std::cout << std::endl;
    std::cout << "--  Perturbation of the CFT hamiltonian  ----" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<numPert;i++)
        {
        std::cout << " phi_(" << r[i] << ")  " << std::endl;
        }
    std::cout << std::endl;

    for (int i=0;i<numPert;i++)
        {
        std::cout << " Conformal dimension of perturbation "<<i+1<<" : "<< deltas[i].first << " " << deltas[i].second << std::endl;
        }
    std::cout << std::endl;

    std::cout << " Sketch of the interaction matrix, each element is a verma module " << std::endl;
    
    for (int i=0;i<n_fields;i++)
        {
        for (int j=0;j<n_fields;j++)
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

    for (int i=0;i<numspace /*n_fields*/;i++)
        {
        cout<<" Subspace "<<i+1<<"      ";
        for (int j=0;j<n_fields;j++)
            {
            if (finalsubspace(i,j) != 0)
                {
                std::cout << finalsubspace(i,j) <<"  ";
                }
            }
        std::cout << std::endl;
        }
    std::cout << std::endl;


    std::cout << " Subspace dimensions (number of verma fields in each subspace) " << std::endl;

    for (int i=0;i<n_fields;i++)
        {
        std::cout << " " << subdim[i];
        }
    std::cout << std::endl;

    std::cout << std::endl;


    for (int i=0;i<numspace;i++)
        {
        std::cout << " Subspace "<<i+1<<", momentum "<< momentum <<" should include "<< subdimStates(i,0) <<" states "<< std::endl;
        std::cout << std::endl;
        } 

}

void minimal_model_parameters::analyse_subspaces()
{
    tensor2D<double> a(n_fields, n_fields);
    // Creating interaction matrix for subspace analysis. The matrix has size n_fields by n_fields, ( n_fileds is the number of fields in theory as from input file).
    // For each perturbation labelled by a combined single index k ( as in cptsmall above), the matrix a(i,j)  is filled with the corresponding strcture constant
    // c_(i,j,k). This will allow to understand if the are symmetry subspaces that can be diagonalized separately. 
    for (int i=0;i<n_fields;i++)
        {
        for (int j=0;j<n_fields; j++)
            {
            a(i,j)=0;
            }
        }

    for (int i=1;i <= n_fields ; i++)
        {
        for (int j=1;j <= n_fields; j++)
            {
            for (int k=0; k< numPert; k++)
                {
                a(i-1, j-1)=a(i-1, j-1) + cptsmall(i, r[k], j);
                }
            } 
        }

    // If the expectation of a certain operator is required, we should check that such operator does not connect different symmetry subspaces. If it does, 
    // its expecation values between different subspaces can be non-zero. In such a case, we should diagonalize together all the subspaces connected by 
    // the perturbation.
    if (calcExpVal)
        {
        for (int i=1;i <= n_fields ; i++)
            {
            for (int j=1;j <= n_fields; j++)
                {
                a(i-1,j-1) = a(i-1,j-1) + cptsmall(i, r_op, j);
                }
            }
        }

    // Tensor intMatrix is just the 0-1 representation of the matrix a
    tensor2D<int> intMatrix_temp(n_fields, n_fields);
    intMatrix = intMatrix_temp;
    for (int i=0;i<n_fields;i++)
        {
        for (int j=0;j<n_fields;j++)
            {
            intMatrix(i,j)=0;
            }
        }

    for (int i=0;i<n_fields;i++)
        {
        for (int j=0;j<n_fields;j++)
            {
            if (std::abs(a(i,j))>0.000001 )
               {
               intMatrix(i,j)=1;
               }
            }
        }

    //Calculation of independent subspaces
    blockrec(n_fields,a,finalsubspace,subdim,numspace);
 
    tensor2D<int>  subdimStates_temp(numspace,1);
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
    std::vector< std::map < std::pair<int, int>, int > > PairToIndex_tmp(numspace);
    PairToIndex = PairToIndex_tmp;
    std::vector< std::map < int, std::pair<int, int> > > IndexToPair_tmp(numspace);
    IndexToPair = IndexToPair_tmp;

        for (int i=0;i<numspace;i++)
            {
            int tempdim=0;
            int vermaCounter = 0;
                tempdim=0;
                for(int j=0;j<subdim[i];j++)
                    {
                    int r_tmp = fields[finalsubspace(i,j) - 1].first.first;
                    int s_tmp = fields[finalsubspace(i,j) - 1].first.second;
                    std::pair <int, int > tmp_pair(r_tmp, s_tmp);
                    if ( PairToIndex[i].find(tmp_pair) == PairToIndex[i].end() ) {
                        //PairToIndex[i][tmp_pair] = 1;
                        PairToIndex[i][tmp_pair] = vermaCounter++;
                    } 
                    //else {
                    //PairToIndex[i][tmp_pair] += 1;
                    //}


                    int r1_tmp = fields[finalsubspace(i,j) - 1].second.first;
                    int s1_tmp = fields[finalsubspace(i,j) - 1].second.second;
                    std::pair <int, int > tmp1_pair(r1_tmp, s1_tmp);
                    if ( PairToIndex[i].find(tmp1_pair) == PairToIndex[i].end() ) {
                        //PairToIndex[i][tmp1_pair] = 1;
                        PairToIndex[i][tmp_pair] = vermaCounter++;
                    } 
                    //else {
                    //   PairToIndex[i][tmp1_pair] += 1;
                    //}     



                    for (int k=0;k<= maxLev;k++)
                        {
                        for (int k1=0;k1<=maxLev;k1++)
                            {
                            if(abs(k-k1) == momentum)
                                {
                                tempdim=tempdim+vermadim(r_tmp-1, s_tmp-1, k)*vermadim(r1_tmp-1, s1_tmp-1, k1);
                                }
                            }
                        }
                    }
                subdimStates(i,0)  = tempdim;
            }

        for (int i=0;i<numspace;i++)
            {
            int countChiralFields = 0;
            std::map < pair<int,int> , int >::iterator it;
            //std::cout << " chiral size subspace " << PairToIndex[i].size() << std::endl;
            for (it = PairToIndex[i].begin(); it !=  PairToIndex[i].end(); it++){
                it->second = countChiralFields;
                std::pair<int, int> tempPair(it->first.first, it->first.second);
                //std::cout<< it->first.first << " " << it->first.second << " index " << it->second << std::endl;
                IndexToPair[i][countChiralFields] = tempPair;
                countChiralFields++;
            }
            //std::cout << " ----------------- " << std::endl;
            } 

        /* for (int i=0;i<numspace;i++)
            {
            std::map < int, pair<int,int> >::iterator it;
            std::cout << " chiral size subspace " << IndexToPair[i].size() << std::endl;
            for(int j=0; j<IndexToPair[i].size(); j++)
                {
                std::pair<int, int> tmp = IndexToPair[i][j];
                std::cout<< tmp.first << " " << tmp.second << " index " << j <<std::endl;
                }
            std::cout << " ----------------- " << std::endl;
            }*/
}
