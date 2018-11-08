/*===========================================================================*/
/* runfunc - running window functions                                        */
/* Adapted from Jarek Tuszynski functions in caTools (archived on CRAN)      */
/* Distributed under GNU General Public License version 3                    */
/*===========================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <float.h>

#include <R.h>
#include <Rinternals.h>

#define notNaN(x)   ((x)==(x))
#define isNaN(x)  (!((x)==(x)))
#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))
#define SQR(x) ((x)*(x))

/* SumErr - macro calculating error of the summing operation */
#define SumErr(a,b,ab) ((((a)>(b)) == ((a)>-(b))) ?  (b) - ((ab)-(a)) : (a) - ((ab)-(b)) )
/* SUM_1 - macro for calculating Sum+=x; Num+=n; Which is NaN aware and have minimal (single number) overflow error correction */
#define SUM_1(x,n, Sum, Err, Num)   if (R_finite(x)){ y=Sum; Err+=x; Sum+=Err; Num+=n; Err=SumErr(y,Err,Sum);  } 
#define mpartial 1024	

void SUM_N(double x, int n, double *partial, int *npartial, int *Num) 
{
  if (R_finite(x)){ 
    int j, i;
    double hi, lo, y;
    for (i=j=0; j<*npartial; j++) {
      y  = partial[j];
      hi = y + x;
      lo = SumErr(x,y,hi); 
      if (lo==0 && i<mpartial) partial[i++] = lo;
      x = hi; 
    }
    partial[i] = x; 
    *npartial   = i+1;
    *Num+=n; 
  }
}  

/*==================================================================================*/
/* Mean function applied to (running) window. All additions performed using         */
/* addition algorithm which tracks and corrects addition round-off errors (see      */  
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps)*/
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runmean(double *In, double *Out, const int *nIn, const int *nWin)
{ /* medium size version with NaN's and edge calculation, but only one level of round-off correction*/
  int i, k2, Num, n=*nIn, m=*nWin;
  double *in, y, *out, Err, Sum;
  double NaN = (0.0/0.0);
  k2  = m>>1;         /* right half of window size */
  in=In; out=Out; 
  Sum = 0;           /* we need to calculate initial 'Sum' */
  Err = 0;
  Num = 0;
  /* step 1 - find mean of elements 0:(k2-1) */      
  for(i=0; i<k2; i++) {
    SUM_1(in[i], 1, Sum, Err, Num)
  }
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m; i++, out++) {
    SUM_1(in[i], 1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
  /* step 3: runsum of the rest of the vector. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m); but with round of error correction */
  for(i=m; i<n; i++, out++, in++) {
    SUM_1(in[m] ,  1, Sum, Err, Num)
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++, out++, in++) {
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
}


/*==================================================================*/
/* minimum function applied to moving (running) window              */ 
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void runmin(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k2, n=*nIn, m=*nWin;
  double ptOut, Min, *in, *out, CST = DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;               /* right half of window size */  
  in  = In;
  out = Out;
  /* --- step 1 - find min of elements 0:(k2-1) */      
  Min=CST;                  /* we need to calculate  initial 'Min' */
  for(i=0; i<k2; i++) Min = MIN(Min,in[i]);  /* find minimum over a window of length k2 */
  /* --- step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m-1; i++) {
    Min=MIN(Min,in[i]);     /* cumulative min */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 3 - the inner section - window of constant size is moving  */      
  ptOut=CST;
  for(i=m-1; i<n; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's min than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m; j++) 
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    } else                  /* if point comining out of the window was NOT window min than min of ... */
      Min=MIN(Min,in[m-1]); /* ... window's first m-1 points is still 'Min', so we have to add a single point */
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's extreme than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m-i-1; j++) 
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    } 
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min);  /* and fill the space with window extreme and move window */
  }
}

#undef MIN
#undef MAX
#undef SQR
#undef SUM_1
#undef SumErr
#undef mpartial
