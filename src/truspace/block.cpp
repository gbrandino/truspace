/**\brief Helper function identify symmetry subspace by the analysis of the interaction matrix   
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */

#include "tensors.h" 
#include <complex>

void  blockrec(int maxl,tensor2D< double >& a, tensor2D< int>& finalsubspace,int *subdim,int &numspace)
{

    int flag1,flagc,flags;
    tensor2D<int> b(maxl,maxl), c(maxl,maxl) ,subspace(maxl,maxl);

    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            b(i,j) = 0;
	        c(i,j) = 0;
	        numspace = 0;
	        subdim[i] = 0;
	        subspace(i,j) = 0;
        	finalsubspace(i,j) = 0;
	        }   
        }

    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            if (i==j)
                {
                b(i,j) = j+1;
                }
            if (std::abs(a(i,j) )>0.000001  and i!=j)
                {
                b(i,j)=j+1;
                }  
            }
        }

    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
	        flag1=0;
            for (int k=0;k<maxl;k++)
                {
                for (int h=0;h<maxl;h++)
                    {
		            if(b(i,k) == b(j,h) and b(i,k) != 0 and b(j,h) != 0 )
		                {
		                flag1=1;
		                }
		            }

                } 
	        if (flag1==1)
	            {
		        for (int li=0;li<maxl;li++)
		            {
		            if (b(i,li) != 0 or b(j,li) != 0)
		                {
		                b(i,li) = li+1;
                        b(j,li) = li+1; 
		                }
	                }
                }
            }
        }


    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            subspace(i,j) = b(i,j);
	        }
       }

    for (int i=0;i<maxl;i++)
        {
        for (int j=i+1;j<maxl;j++)
            {
            flagc=0;
            for (int k=0;k<maxl;k++)
                {
                if (subspace(i,k) == subspace(j,k))
                    {
                    flagc=flagc+1;
                    }
                }
            if (flagc==maxl)
                {
                for (int k=0;k<maxl;k++)
                    {
                    subspace(j,k)=0;
                    }
                }    
            }
        }

    int h=0;
    for (int i=0;i<maxl;i++)
        {
        int k=0;
        flags=0;
        for (int j=0;j<maxl;j++)
            {
	        if (subspace(i,j) != 0)
	            {
	            c(h,k) = subspace(i,j);
	            k=k+1;
                flags=1;
	            }
	        }
        if (flags==1)
            {
            h=h+1;
            }
        }
	

    for (int i=0;i<maxl;i++)
        {
        int k=0;
        for (int j=0;j<maxl;j++)
            {
	        if (c(i,j) != 0)
	            {
                finalsubspace(i,k) = c(i,j);
                k=k+1;
                }
            }
        }
     
    for (int i=0;i<maxl;i++)
        {
        for (int j=0;j<maxl;j++)
            {
            if (finalsubspace(i,j) != 0)
	            {
	            subdim[i]=subdim[i]+1;
	            }
	        }
        }


    for (int i=0;i<maxl;i++)
        {
        if (subdim[i] != 0)
            {
            numspace=numspace+1;
            }
        }     

}
