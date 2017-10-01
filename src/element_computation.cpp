/*
 * Routine to compute matrix elements between two chiral states of an arbitrary perturbation.
 * The perturbation is characterized solely by its conformal weight and the central charge C
 * of the theory.  The input to the routine takes the following form:
 * 
 * grade[MaxNumGrade]: contains the string of modes + operator; 
 * n = total number of modes plus operator
 * m = position of operator
 * hr = conformal weight of right vacuum state
 * hl = conformal weight of left vacuum state
 * Del = conformal weight of perturbing operator \phi_Del
 *
 * for example <hr|L_{k1}L_{k2}\phi_Del(0)L_{l2}L_{l1}|hl> would be encoded as
 *
 * n = 5
 * m = 3
 * grade[1] = l1
 * grade[2] = l2
 * grade[3] = 0
 * grade[4] = k2
 * grade[5] = k1
 * 
 * If one want to compute the overlap of two states, one simple sets Del = 0, i.e. the identity
 * operator.
 *
 * The routine works by recursive calls.  At each call it performs an elementary operation on
 * the matrix element (i.e. some commutation of modes) and either passes the results (simpler matrix
 * elements) onto a subsequent call of element_computation or, if the resulting matrix element is simple enough,
 * directly evaluates it.
 *
 * For a given matrix element, element_computation first attempts to move all modes to right of the operator.
 * If this is already done, it will find the rightmost mode with positive grade and commute it with
 * its right neighbour.  

 * The routine tries to deal with special cases in order to speed up its performance.  For example,
 * if there is a string of modes of the form L_{-n}L_{n}^k, it will compute not [L_{-n},L_{n}]
 * but the results of commuting L_{-n} through L_{n}^k in one go.
 *
 * authors: Robert Konik - BNL, Giuseppe Piero Brandino - eXact-lab s.r.l. 
 * last revision: November 2016
 */

#include "element_computation.h"
#include <cmath>
#include <iostream>    

long double element_computation(long double C, long double* grade,int n,int m,long double hr,long double hl,long double coeff,long double Del)
{
    
  int i, j, rightmode; //k, length, flag, neg_mode[MaxNumGrade], pos_mode[MaxNumGrade];
  
  long double grade_res[MaxNumGrade], coeff_res, c, ct, temp, temp1, tgrade; //term, t1grade, temp2;
  
  long double Gam[MaxNumGrade], coeffs[MaxNumGrade], wr, wl, diff, gr; //Sum_IGam[MaxNumGrade]
  
  int len, pt; // len_r, len_l, min;     
  
  long double nterms, gam[MaxNumGrade][MaxNumGrade], alpha[MaxNumGrade][MaxNumGrade];
  int len_m1, len_1;
  

 
  if (n>(MaxNumGrade - 2))
    std::cout << " Warning: element too long" << std::endl;
  
  if (fabsl(coeff) < 0.000000000000001)
    return(0.);
  
  tgrade = 0.;
  for(i=1;i<=(m-1);++i) {
    tgrade += grade[i];			
    if (tgrade > .00001)
      return(0.);
  }
  
  
  /*move modes to the left of the operator to the right*/
  
  if (n > m) {
    
    /*determine length of string of modes with the same grade
      to the left of the operator*/
    
    len = 1;
    pt = m+2;
    diff = fabsl(grade[pt]-grade[m+1]);
    while ((diff<.00000001) && (pt<=n)) {
      pt++;
      len++;
      diff = fabsl(grade[pt]-grade[m+1]);
    }

    /*define coefficients*/
    
    wr = hr;
    for(i=1;i<m;++i) 
      wr -= grade[i];
    wl = hl;
    for(i=m+len+1;i<=n;++i) 
      wl += grade[i];
    
    for(i=1;i<=len;++i) {
      Gam[i] = wl-wr+grade[m+1]*Del+(len-i)*grade[m+1];
    }
    
    /*now commute L^(len)_k through the operator; set up resulting matrix elements*/
    
    temp = 0.;
    for(i=0;i<=len;++i) {
      coeffs[i] = 1.;
      for(j=1;j<=(len-i);++j) {
	coeffs[i] *= Gam[j]; 
      }
      coeffs[i] *= fac(len)/fac(len-i)/fac(i);
      
      for(j=1;j<m;++j) {
	grade_res[j] = grade[j];
      }
      for(j=1;j<=i;++j) {
	grade_res[m-1+j] = grade[m+1];
      }
      grade_res[m+i] = 0.;
      
      for(j=1;j<=(n-m-len);++j) {
	grade_res[m+i+j] = grade[m+len+j];
      }
      temp += element_computation( C, grade_res,n-len+i,m+i,hr,hl,coeff*coeffs[i],Del);
    }
    return(temp);
  }
  
  /* all modes are to the right of the operator*/
  
  if (n == m) {
    
    if (n == 1)
      return(coeff);
    
    /* find most rightward mode with non-negative grade */
    
    rightmode = 0;
    for(i=1;i<n;++i) {
      if (grade[i] > -.5) {
	rightmode = i;
	break;
      }
    }
    
    /*check if sum of grades of non-negative modes exceeds that of negative modes*/
    
    tgrade = 0.;
    for(i=1;i<n;++i) {
      tgrade += grade[i];
    }
    if (tgrade > 0.0001)
      return(0.);
    
    
    /* all modes have negative grade */
    
    if (rightmode == 0) {
      coeff_res = coeff;
      for(i=n-1;i>=1;--i) {
	c = -grade[i]*Del+hr-hl;
	for(j=1;j<i;++j) {
	  c -= grade[j];
	}
	coeff_res *= c;
      }
      return(coeff_res);
    }
    
    /* positive grade mode is adjacent to highest weight state*/
    
    if (rightmode == 1) {
      return(0.);
    }
    
    /*or it isn't*/
    
    /*special case: commute L_{n} through a string of L_{-n}'s*/
    
    else if ((rightmode > 1) && (fabsl(grade[rightmode]+grade[rightmode-1]) < .0001)) {
      
      gr = grade[rightmode];
      ct = (powl(gr,3.)-gr)*C/12.;

      /*determine length of string of L_{-n}'s*/
      
      len_m1 = 0;
      for(i=rightmode-1;i>=1;i--) {
	if (fabsl(grade[i]+gr) < .001) 
	  len_m1++;
	else
	  break;
      }
      
      /*determine length of string of L_n's*/
      
      len_1 = 0;
      for(i=rightmode;i<m;i++) {
	if (fabsl(grade[i]-gr) < .001) 
	  len_1++;
	else
	  break;
      }
      
      /*determine weight of state to the right of all L_{-n}s*/
      
      wr = hr;
      for(i=1;i<(rightmode-len_m1);i++) {
	wr -= grade[i];
      }
      
      /*compute coefficients of gam[k][j] = sum^k_l=1 2*(wr+l-1-j)*/
      
      gam[1][0] = 2.*gr*wr+ct;
      for(i=2;i<=len_m1;++i) 
	gam[i][0] = gam[i-1][0]+2.*gr*(wr+gr*(i-1.))+ct;
      
      for(j=0;j<=(len_1-1);++j) {
	for(i=1;i<=len_m1;++i) 
	  gam[i][j] = gam[i][0]-2.*gr*gr*j*i;
      }
      
      
      /*number of terms results from commutations*/
      
      nterms = (len_m1 < len_1)?len_m1:len_1;
      
      for(i=1;i<(rightmode-len_m1);i++) 
	grade_res[i] = grade[i];
      for(i=rightmode-len_m1;i<(rightmode-len_m1+len_1);++i)
	grade_res[i] = gr;
      for(i=rightmode-len_m1+len_1;i<(rightmode+len_1);++i)
	grade_res[i] = -gr;
      for(i=rightmode+len_1;i<=n;++i)
	grade_res[i] = grade[i];
      temp1 = element_computation( C, grade_res,n,m,hr,hl,coeff,Del);	      
      
      
      for(j=0;j<=len_1;++j) 
	alpha[0][j] = 1.;
      
      for(i=1;i<=nterms;++i) {
	for(j=0;j<=(len_1-i);++j) {
	  if (j==0)
	    alpha[i][j] = gam[len_m1-i+1][j]*alpha[i-1][j];
	  else if (j>0)
	    alpha[i][j] = gam[len_m1-i+1][j]*alpha[i-1][j]+alpha[i][j-1];
	}
	
	coeff_res = alpha[i][len_1-i];
	
	for(j=1;j<(rightmode-len_m1);j++) 
	  grade_res[j] = grade[j];
	for(j=rightmode-len_m1;j<(rightmode-len_m1+len_1-i);++j)
	  grade_res[j] = gr;
	for(j=rightmode-len_m1+len_1-i;j<(rightmode+len_1-2*i);++j)
	  grade_res[j] = -gr;
	for(j=rightmode+len_1-2*i;j<=(n-2*i);++j)
	  grade_res[j] = grade[(j-rightmode-len_1+2*i)+rightmode+len_1];
	
	temp1 += element_computation( C, grade_res,n-2*i,m-2*i,hr,hl,coeff_res*coeff,Del);	      
	
      }
      return(temp1);
      
    }
    
    /*general case*/
    
    else if (rightmode > 1) {
      
      /*determine length of string of modes with the same grade
	as rightmode*/
      
      len = 1;
      pt = rightmode+1;
      diff = grade[pt]-grade[rightmode];
      diff *= (diff<0.)?-1.:1.;
      while ((diff<.00000001) && (pt<=n)) {
	pt++;
	len++;
	diff = grade[pt]-grade[rightmode];
	diff = (diff<0.0)?-1.*diff:diff;
      }
      
      /*now commute string of modes beginning with rightmode through rightmode-1*/
      
      temp = temp1 = 0.;
      
      /*string of modes annihilate rightmode-1 through the commutation*/

      if (fabsl(grade[rightmode]+grade[rightmode-1])<.0000001) {
	
	/*string commutes throught mode at rightmode-1*/
	
	for(i=1;i<MaxNumGrade;++i) {
	  grade_res[i] = grade[i];
	}
	grade_res[rightmode+len-1] = grade[rightmode-1];
	grade_res[rightmode-1] = grade[rightmode];
	
	temp = element_computation( C, grade_res,n,m,hr,hl,coeff,Del);	      
	
	
	/*mode at rightmode-1 is annihilated in commuting through string*/
	
	wr = hr;
	for(i=1;i<(rightmode-1);++i) {
	  wr -= grade[i];
	}
	c = 0.;
	ct = C/12.*(powl(grade[rightmode],3.)-grade[rightmode]);
	for(i=0;i<len;++i) {
	  c += (2.*grade[rightmode]*(wr-grade[rightmode]*i)+ct);
	}
	
	for(i=1;i<(rightmode-1);++i) {
	  grade_res[i] = grade[i];
	}
	for(i=rightmode+1;i<MaxNumGrade;++i) {
	  grade_res[i-2] = grade[i];
	}		
	
	temp1 = element_computation( C, grade_res,n-2,m-2,hr,hl,c*coeff,Del);	      
	
      }
      
      /*string at rightmode and rightmode-1 do not annihilate one another*/
      
      else if (fabs(grade[rightmode]+grade[rightmode-1])>.0000001) {	    
	
	temp1 = 0.;
	for(i=0;i<=len;++i) {
	  
	  coeffs[i] = 1.;
	  for(j=1;j<=i;++j) {
	    coeffs[i] *= (grade[rightmode]*(2.-j)-grade[rightmode-1]); 
	  }
	  coeffs[i] *= fac(len)/fac(len-i)/fac(i);
	  
	  wr = hr;
	  for(j=1;j<(rightmode-1);++j) {
	    grade_res[j] = grade[j];
	    wr -= grade[j];
	  }
	  for(j=1;j<=(len-i);++j) {
	    grade_res[rightmode+j-2] = grade[rightmode];
	    wr -= grade[rightmode];
	  }
	  
	  diff = grade[rightmode-1]+i*grade[rightmode];
	  diff *= (diff<0.)?-1.:1.;
	  
	  if (diff > 0.0001) {
	    grade_res[rightmode+len-i-1] = grade[rightmode-1]+i*grade[rightmode];
	    for(j=(rightmode+len);j<=n;++j) {
	      grade_res[rightmode+len-i+(j-rightmode-len)] = grade[j];
	    }
	    temp1 += element_computation( C, grade_res,n-i,m-i,hr,hl,coeff*coeffs[i],Del);
	  }
	  else if (diff < .0001) {
	    for(j=(rightmode+len);j<=n;++j) {
	      grade_res[rightmode+len-i-1+(j-rightmode-len)] = grade[j];
	    }
	    ct = C/12.*(grade[rightmode]*grade[rightmode]*grade[rightmode]-grade[rightmode])/(2.*grade[rightmode]);
	    temp1 += element_computation( C, grade_res,n-i-1,m-i-1,hr,hl,coeff*(wr*coeffs[i]+coeffs[i]*ct),Del);
	  }
	  
	}
      }
      return(temp+temp1);
    }
  }
  return 0; //CHECK
}

/*compute n factorial*/
long double fac(int n)
{
    long double temp;

    temp = 1.;
    for (int i=1;i<=n;++i) 
        {
        temp = temp*((long double)i);
        }
  return(temp);
}
