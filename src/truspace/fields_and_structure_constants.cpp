#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include "cpt.h"
#include "tensors.h"
int main(int argc, char** argv)
{


    if (argc < 3)
        {
        std::cout<<" Usage: " << argv[0] << " [p] [q] " << std::endl;
        std::exit(0);
        }


    int p = atoi(argv[1]);
    int q = atoi(argv[2]);
    //int maxLev = atoi(argv[3]);
    //int minLev = atoi(argv[4]);

    double C=1.0-6.0*std::pow(p-q,2)/double(p*q);

    int maxp;
    int maxq;
    int maxl;
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


    tensor2D<double> Kac(maxp,maxq);
    
    //temporary tensor, used to properly allocate cpt now that we know the sizes
    tensor6D< std::complex<long double>  > cpt(maxp, maxq, maxp, maxq, maxp, maxq);

    //temporary tensor, used to properly allocate cptsmall now that we know the sizes
    tensor3D< std::complex<long double>  > cptsmall(maxl, maxl, maxl);

    std::stringstream fields_out_filename;
    fields_out_filename << "p_" << p << "_q_" << q << "_fields_diagonal.dat";
    std::ofstream fields_outstream(fields_out_filename.str().c_str(), std::ofstream::out);


    // filling Kac table with primary fields weight
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            Kac(i-1,j-1)=hrs((long double)i,(long double)j,p,q);
            fields_outstream << i << "," << j << "," << i << "," << j << " # " << Kac(i-1,j-1) << std::endl;
            }
        }
    fields_outstream.close();


    std::complex<long double> tdot;

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
                            if (std::isnan(real(tdot)) or std::isinf(real(tdot)) or std::isnan(imag(tdot)) or std::isinf(imag(tdot)) )
                                {
                                tdot=dot((int)p,(int)q,i1+1,i+1,k1+1,k+1,j1+1,j+1);   //132
                                if (std::isnan(real(tdot)) or std::isinf(real(tdot)) or std::isnan(imag(tdot)) or std::isinf(imag(tdot)) )
                                   {
                                   tdot=dot((int)p,(int)q,j1+1,j+1,i1+1,i+1,k1+1,k+1);  //213
                                   if (std::isnan(real(tdot)) or std::isinf(real(tdot)) or std::isnan(imag(tdot)) or std::isinf(imag(tdot)) )
                                      {
                                      tdot=dot((int)p,(int)q,j1+1,j+1,k1+1,k+1,i1+1,i+1);  //231  
                                      if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) or std::isinf(imag(tdot)) )
                                         {
                                         tdot=dot((int)p,(int)q,k1+1,k+1,j1+1,j+1,i1+1,i+1); //321
                                         if (std::isnan(real(tdot)) or std::isinf(real(tdot)) or std::isnan(imag(tdot)) or std::isinf(imag(tdot)) )
                                            {
                                            tdot=dot((int)p,(int)q,k1+1,k+1,i1+1,i+1,j1+1,j+1); //312 
                                            if (std::isnan(real(tdot)) or std::isinf(real(tdot)) or std::isnan(imag(tdot)) or std::isinf(imag(tdot)) ) 
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


    std::stringstream cptsmall_out_filename;
    cptsmall_out_filename << "p_" << p << "_q_" << q << "_constant.dat";
    std::ofstream cptsmall_outstream(cptsmall_out_filename.str().c_str(), std::ofstream::out);

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
                            if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
                                {
                                tdot=dot((int)p,(int)q,i1+1,i+1,k1+1,k+1,j1+1,j+1);   //132
                                if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
                                    {
                                    tdot=dot((int)p,(int)q,j1+1,j+1,i1+1,i+1,k1+1,k+1);  //213
                                    if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
                                        {
                                        tdot=dot((int)p,(int)q,j1+1,j+1,k1+1,k+1,i1+1,i+1);  //231  
                                        if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
                                            {
                                            tdot=dot((int)p,(int)q,k1+1,k+1,j1+1,j+1,i1+1,i+1); //321
                                            if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
                                                {
                                                tdot=dot((int)p,(int)q,k1+1,k+1,i1+1,i+1,j1+1,j+1); //312 
                                                if (std::isnan(real(tdot)) || std::isinf(real(tdot)) || std::isnan(imag(tdot)) || std::isinf(imag(tdot)) )
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
                            //std::cout << typeid(tdot.real).name() << std::endl;
                            //std::cout << tdot.real() << std::endl; 
                            if (std::abs(tdot.real()) > 1e-16 )
                                cptsmall_outstream << ai + 1  << "," << bi + 1 << "," << ci + 1 << "|" << std::real(tdot) << std::endl;
                            ci=ci+1;
                            }
                        }
                    bi=bi+1;
                    }
                }
            ai=ai+1;
            }
        } 

    cptsmall_outstream.close();

    std::cout << "-------------------------" << std::endl;
    std::cout << " Mininal model ( " << p <<  ", "  << q <<  " ) " << std::endl;
    std::cout << std::endl;
    std::cout << " Central charge " << C << std::endl;
    std::cout << std::endl;
    std::cout << " Number of primary fields " << maxl << std::endl;
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




} 
