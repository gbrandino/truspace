/**  Set of helper function to calculate structure constant of Minimal CFTs.  
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */

#include "cpt.h"
#include "gamma.h"
#include <complex>
#include <cstdlib>

using namespace std;

long double hrs(long double r,long double s,long double p, long double q) 
    { 
    return (pow((r*q-s*p),2)- pow((p-q),2))/(4*p*q);
    }



long double dotmu(int ll,int l,long double ro,long double ro1)
    {
    long double prod1,prod2,prod3;
    prod1=1;
    prod2=1;
    prod3=1;
//    if (ll-1>0 && l-1>0)
//        {
    for (int i=1;i<=ll-1;i++)
        {
        for (int j=1;j<=l-1;j++)
            {
            prod1=prod1*pow((long double)i-((long double)j)*ro,-2);                
            }
        }
        for (int i=1;i<=ll-1;i++)
            {
            prod2=prod2*gamma(((long double)i)*ro1)/gamma(((long double)1)-((long double)i)*ro1);
            }
        for (int j=1;j<=l-1;j++)
            {
            prod3=prod3*gamma(((long double)j)*ro)/gamma(((long double)1)-((long double)j)*ro);
            }
       return pow(ro,(long double)(4*(ll-1)*(l-1)))*prod1*prod2*prod3;
//      }
//      else 
//          {
//          return 0;
//          }
    }


long double dotc(int ll,int l,long double ro, long double ro1, long double kn1, long double kn, long double ks1, long double ks, long double kp1, long double kp )
    {
    long double prod1,prod2,prod3;
    prod1=1;
    prod2=1;
    prod3=1;
//    if (ll-2>=0 && l-2>=0)
//      { 
    for (int i=0;i<=ll-2;i++)
        {
        for (int j=0;j<=l-2;j++)
            {
            prod1=prod1*pow(pow((long double)(ks1-1-i)-ro*((long double)(ks-1-j)) ,2)*pow((long double)(kn1-1-i)-ro*((long double)(kn-1-j)) ,2)*pow((long double)(kp1+1+i)-ro*((long double)(kp+1+j)) ,2),-1);
            }
        }
        for (int i=0;i<=ll-2;i++)
            {
            prod2=prod2*gamma(-ro1*((long double)(ks1-1-i))+(long double)ks)*gamma(-ro1*((long double)(kn1-1-i))+(long double)kn)*gamma(ro1*((long double)(kp1+1+i))-(long double)kp)/(gamma(1+ro1*((long double)(ks1-1-i))-(long double)ks)*gamma(1+ro1*((long double)(kn1-1-i))-(long double)kn)*gamma(1-ro1*((long double)(kp1+1+i))+(long double)kp)) ;
            }
        for (int j=0;j<=l-2;j++)
            {
            prod3=prod3*gamma(-ro*((long double)(ks-1-j))+(long double)ks1)*gamma(-ro*((long double)(kn-1-j))+(long double)kn1)*gamma(ro*((long double)(kp+1+j))-(long double)kp1)/(gamma(1+ro*((long double)(ks-1-j))-(long double)ks1)*gamma(1+ro*((long double)(kn-1-j))-(long double)kn1)*gamma(1-ro*((long double)(kp+1+j))+(long double)kp1)) ;
            }
    return dotmu(ll,l,ro,ro1)*prod1*prod2*prod3;
//      }
//      else 
//          {
//          return 0;
//          }
    }

long double dota(int kn1,int kn,long double ro,long double ro1)
    {
    long double prod1,prod2,prod3;
    prod1=1;
    prod2=1;
    prod3=1;
    for (int i=1;i<=kn1-1;i++)
        {
        for (int j=1;j<=kn-1;j++)
            {
            prod1=prod1*pow(((long double)(1+i))-ro*((long double)(1+j)),2)/pow((long double)i-((long double)j)*ro,2);
            }
        }
    for (int i=1;i<=kn1-1;i++)
        {
        prod2=prod2*gamma(((long double)i)*ro1)*gamma(2-ro1*((long double)(1+i)))/(gamma(((long double)1)-((long double)i)*ro1)*gamma(-1+ro1*((long double)(1+i))));
        }
    for (int j=1;j<=kn-1;j++)
        {
        prod3=prod3*gamma(((long double)j)*ro)*gamma(2-ro*((long double)(1+j)))/(gamma(((long double)1)-((long double)j)*ro)*gamma(-1+ro*((long double)(1+j))));
        }
    return prod1*prod2*prod3;
    }
 






complex<long double> dot(int v,int v1,int n1,int n, int s1,int s,int p1,int p)
       {
       if (div(n1+s1-p1+1,2).rem==0 && div(n+s-p+1,2).rem==0 && (n1+s1-p1+1)>0 && (n+s-p+1)>0)
       {
       
       return
      
dotc(div(n1+s1-p1+1,2).quot,div(n+s-p+1,2).quot,((long double)v1)/((long double)v),((long double)v)/((long double)v1),n1,n,s1,s,p1,p)*pow((complex<long double>)(dota(n1,n,((long double)v1)/((long double)v),((long double)v)/((long double)v1))*dota(s1,s,((long double)v1)/((long double)v),((long double)v)/((long double)v1))/dota(p1,p,((long double)v1)/((long double)v),((long double)v)/((long double)v1))),(long double)0.5);
       }
       else if  (div(n1+s1+p1-v1+1,2).rem==0 && div(n+s+p-v+1,2).rem==0 && (n1+s1+p1-v1+1)>0 && (n+s+p-v+1)>0)  
            {
             return
	    
dotc(div(n1+s1+p1-v1+1,2).quot,div(n+s+p-v+1,2).quot,((long double)v1)/((long double)v),((long double)v)/((long double)v1),n1,n,s1,s,v1-p1,v-p)*pow((complex<long double>)(dota(n1,n,((long double)v1)/((long double)v),((long double)v)/((long double)v1))*dota(s1,s,((long double)v1)/((long double)v),((long double)v)/((long double)v1))/dota(v1-p1,v-p,((long double)v1)/((long double)v),((long double)v)/((long double)v1))),(long double)0.5);
         
            }
            else 
            {
            return 0;
            }
       }





