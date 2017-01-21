/* Set of helper functions manage ginac computation 
 * 
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l.
 * last revision: October 2016
 *          
 */

#include "characters.h"

ex invphi(int k,symbol x)
    {
    ex a=1;
    ex b;
    for (int i=1;i<=k;i++)
        {
        a=a/(1-pow(x,i));
        }
    b=a.series(x==0,k);
    return series_to_poly(b);
    }

ex ksymb1(int r,int s, int p, int q,int k,symbol x)
    {
    ex sum=0;
    for (int n=-k;n<=k;n++)
        {
        sum=sum+pow(x,p*q*n*n+n*(r*q-s*p));
        }
    return sum;
    }

ex ksymb2(int r,int s, int p, int q,int k,symbol x)
    {
    ex sum=0;
    for (int n=-k;n<=k;n++)
        {
        sum=sum+pow(x,p*q*n*n+n*(r*q+s*p)+s*r);
        }
    return sum;
    }



ex chisymb(int r,int s,int p, int q,int k,symbol x)
    {
    return expand((ksymb1(r,s,p,q,k,x)-ksymb2(r,s,p,q,k,x))*invphi(k,x));
    }



void characexpan(int r,int s, double p,  double q,int l,int* leveldim ) 
    {
    symbol x("x");
    for (int k=0;k<=l;k++)
        {
        leveldim[k]=to_int(ex_to<numeric>(chisymb(r,s,(int)p,(int)q,l+3,x).coeff(x,k)));
        }
    }  


int partitions(int k)
    {
    ex a;
    symbol x("x");
    a=invphi(k+3,x);
    return to_int(ex_to<numeric>(a.coeff(x,k)));
    }

/* 
 * routine to determine partition of integers 
 * the routine accepts a given integer (i.e. "level" and returns
 * the partitions in results
 * 
 * This routine can probably be optimize or omp parallized, but with deep rewriting
 *
 */
int determine_partition(int level, tensor2D<int>&  results)
{

    int num,  flag, total;
    int maxpart = results.n1;
    int maxgrade = results.n2;
    tensor2D<int> results1(maxpart,maxgrade);

    for(int i=0;i< maxpart;++i) 
        {
        for(int j=0;j< maxgrade;++j) 
            {
            results1(i,j) = 0;
            }
        }

    if (level==0) 
        {
        return 0;
        }

    if (level==1) 
        {
        results(0,0) = 1;
        return 1;
        }

    total = 1;
    results(0,0) = level;
    for(int i=level-1;i>=1;i--) 
        {
        num = determine_partition(level-i,results1);
        for(int j=0;j<num;++j)      
            {
            flag = 1;
            for(int k=0;k<maxgrade;++k) 
                {
                if (results1(j,k) > i) 
                    {
                    flag = 0;
                    break;
                    }
                }
            if (flag == 1) 
                {
                results(total,0) = i;
                for(int k=1;k<maxgrade;++k)
                    results(total,k) = results1(j,k-1);
                total++;
                }
            }   
        }

    return total;
}
