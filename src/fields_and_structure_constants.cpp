#include <iostream>
#include "cpt.cpp"
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
    tensor6D< complex<long double>  > cpt(maxp, maxq, maxp, maxq, maxp, maxq);

    //temporary tensor, used to properly allocate cptsmall now that we know the sizes
    tensor3D< complex<long double>  > cptsmall(maxl, maxl, maxl);

    // filling Kac table with primary fields weight
    for (int i=1;i<=maxp;i++)
        {
        for (int j=1;j<=maxq;j++)
            {
            Kac(i-1,j-1)=hrs((long double)i,(long double)j,p,q);
            }
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




} 
