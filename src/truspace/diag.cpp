/** Class containg a set of c-programs that diagonalize matrices based on code found in Numerical Recipes, for both double and long double. 
 * \author Robert Konik - BNL, Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2016
 */

#include "diag.h"
#include "tensors.h"
#include <math.h>
#include "minimal_model_config.h"

void diag::shell_ld(unsigned long n,long double d[],long double* r[])
{

  unsigned long i,j,l,inc;
  long double v, ve[n+1];


  inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {
    inc /= 3;
    for(i=inc+1;i<=n;i++) {
      v=d[i];
      for(l=0;l<n;++l) {
	  ve[l+1] = r[l+1][i];
      }
      j=i;
      while (d[j-inc] > v) {
	d[j]=d[j-inc];
	for(l=0;l<n;++l) {
	    r[l+1][j] = r[l+1][j-inc];
	}
	j -= inc;
	if (j <= inc) 
	  break;
      }
      d[j]=v;
      for(l=0;l<n;++l) {
	  r[l+1][j] = ve[l+1];
      }
    }
  } while (inc > 1);
}

void diag::tqli_ld(long double* d, long double*  e,int  n,tensor2D<long double >& z,int  ev_flag)
{
  int m,l,iter,i,k;
  long double s,r,p,g,f,dd,c,b,temp,temp1;

  for(i=2;i<=n;i++)
    e[i-1] = e[i];
  e[n] = 0.0;
  for(l=1;l<=n;l++) {
    iter = 0;
    do {
      for(m=l;m<=n-1;m++) {
	dd=fabsl(d[m])+fabsl(d[m+1]);
	temp = fabsl(e[m]+dd);
	temp1 = fabsl(temp-dd);
	if (temp1 < .000000000000000001)
	  break;
      }
      if (m != l) {
	if (iter++ == 30)
	  printf("Too many iterations in tqli\n");
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = sqrtl(g*g+1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	
	for(i=m-1;i>=l;i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r = sqrtl(f*f+g*g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i]-g)*s+2.0*c*b;
	  d[i+1] = g + (p=s*r);
	  g = c*r - b;
	  
	  if (ev_flag == 1) {
	    for (k=1;k<=n;++k) {
	      f = z(k,i+1);
	      z(k,i+1) = s*z(k,i) + c*f;
	      z(k,i) = c*z(k,i) - s*f;
	    }
	  }
	}
	if ((r == 0.0) && (i >= l))
	    continue;
	
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}	



void diag::tred2_ld(tensor2D< long double >& a,int  n,long double* d ,long double* e ,int ev_flag)
{
  int l,k,j,i;
  long double scale,hh,h, g,f;


  for(i=n;i>=2;i--) {
    l = i-1;
    h = scale = 0.0;
    if (l>1) {
      for(k=1;k<=l;k++) 
	scale += fabsl(a(i,k));
      if (scale == 0.0)
	e[i] = a(i,l);
      else {
	for(k=1;k<=l;k++) {
	  a(i,k) /= scale;
	  h += a(i,k)*a(i,k);
	}
	f = a(i,l);
	g = (f >= 0.0 ? -sqrtl(h):sqrtl(h));
	e[i] = scale*g;
	h -= f*g;
	a(i,l) = f - g;
	f = 0.0;
	for(j=1;j<=l;j++) {
	  if (ev_flag == 1)
	    a(j,i) = a(i,j)/h;
	  g = 0.0;
	  for(k=1;k<=j;k++) 
	    g += a(j,k)*a(i,k);
	  for(k=j+1;k<=l;k++)
	    g += a(k,j)*a(i,k);
	  e[j] = g/h;
	  f += e[j]*a(i,j);
	}
	hh = f/(h+h);
	for(j=1;j<=l;j++) {
	  f = a(i,j);
	  e[j] = g = e[j] - hh*f;
	  for (k=1;k<=j;k++) 
	    a(j,k) -= (f*e[k] + g*a(i,k));
	}
      }
    }
    else
      e[i] = a(i,l);
    d[i] = h;
  }

  d[1] = 0.0;
  e[1] = 0.0;

  for(i=1;i<=n;i++) {
    if (ev_flag == 1) {
      l = i-1;
      if (d[i]) {
	for(j=1;j<=l;j++) {
	  g = 0.0;
	  for(k=1;k<=l;k++) 
	    g += a(i,k)*a(k,j);
	  for(k=1;k<=l;k++)
	    a(k,j) -= g*a(k,i);
	}
      }
    }
    d[i] = a(i,i);
    if (ev_flag == 1) {
      a(i,i) = 1.;
      for(j=1;j<=l;j++) {
	a(j,i)=a(i,j) = 0.0;
      }
    }
  }
}


void diag::shell_double(unsigned long n,double d[] , double *r[])
{

  unsigned long i,j,l,inc;
  double v, ve[n+1];


  inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {
    inc /= 3;
    for(i=inc+1;i<=n;i++) {
      v=d[i];
      for(l=0;l<n;++l) {
	  ve[l+1] = r[l+1][i];
      }
      j=i;
      while (d[j-inc] > v) {
	d[j]=d[j-inc];
	for(l=0;l<n;++l) {
	    r[l+1][j] = r[l+1][j-inc];
	}
	j -= inc;
	if (j <= inc) 
	  break;
      }
      d[j]=v;
      for(l=0;l<n;++l) {
	  r[l+1][j] = ve[l+1];
      }
    }
  } while (inc > 1);
}

void diag::tqli_double(double d[],double e[] ,int n,double *z[] ,int ev_flag)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b,temp,temp1;

  for(i=2;i<=n;i++)
    e[i-1] = e[i];
  e[n] = 0.0;
  for(l=1;l<=n;l++) {
    iter = 0;
    do {
      for(m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	temp = fabs(e[m]+dd);
	temp1 = fabs(temp-dd);
	if (temp1 < .0000000000000000000001)
	  break;
      }
      if (m != l) {
	if (iter++ == 30)
	  printf("Too many iterations in tqli\n");
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = sqrt(g*g+1.0);
	g = d[m]-d[l]+e[l]/(g+SIGNd(r,g));
	s = c = 1.0;
	p = 0.0;
	
	for(i=m-1;i>=l;i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r = sqrt(f*f+g*g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i]-g)*s+2.0*c*b;
	  d[i+1] = g + (p=s*r);
	  g = c*r - b;
	  
	  if (ev_flag == 1) {
	    for (k=1;k<=n;++k) {
	      f = z[k][i+1];
	      z[k][i+1] = s*z[k][i] + c*f;
	      z[k][i] = c*z[k][i] - s*f;
	    }
	  }
	}
	if ((r == 0.0) && (i >= l))
	    continue;
	
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}	



void diag::tred2_double(double *a[],int n,double d[], double e[],int ev_flag)
{
  int l,k,j,i;
  double scale,hh,h, g,f;

  for(i=n;i>=2;i--) {
      /*printf("tred2: %d\n",i);fflush(stdout);*/
    l = i-1;
    h = scale = 0.0;
    if (l>1) {
      for(k=1;k<=l;k++) 
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i] = a[i][l];
      else {
	for(k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f = a[i][l];
	g = (f >= 0.0 ? -sqrt(h):sqrt(h));
	e[i] = scale*g;
	h -= f*g;
	a[i][l] = f - g;
	f = 0.0;
	for(j=1;j<=l;j++) {
	  if (ev_flag == 1)
	    a[j][i] = a[i][j]/h;
	  g = 0.0;
	  for(k=1;k<=j;k++) 
	    g += a[j][k]*a[i][k];
	  for(k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j] = g/h;
	  f += e[j]*a[i][j];
	}
	hh = f/(h+h);
	for(j=1;j<=l;j++) {
	  f = a[i][j];
	  e[j] = g = e[j] - hh*f;
	  for (k=1;k<=j;k++) 
	    a[j][k] -= (f*e[k] + g*a[i][k]);
	}
      }
    }
    else
      e[i] = a[i][l];
    d[i] = h;
  }

  d[1] = 0.0;
  e[1] = 0.0;

  for(i=1;i<=n;i++) {
    if (ev_flag == 1) {
      l = i-1;
      if (d[i]) {
	for(j=1;j<=l;j++) {
	  g = 0.0;
	  for(k=1;k<=l;k++) 
	    g += a[i][k]*a[k][j];
	  for(k=1;k<=l;k++)
	    a[k][j] -= g*a[k][i];
	}
      }
    }
    d[i] = a[i][i];
    if (ev_flag == 1) {
      a[i][i] = 1.;
      for(j=1;j<=l;j++) {
	a[j][i]=a[i][j] = 0.;
      }
    }
  }
}

