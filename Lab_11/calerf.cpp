#include "stdafx.h"
#include "math.h"
#include "calerf.h"

// -------- netlib package by J. W. Cody ---------
// Translated into C++ by L. Bieniasz
// -----------------------------------------------


double calerf(const double arg, const int jint)
{
//-----------------------------------------------------------------
//
//  This packet evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
//  for a real argument x.  It contains four double functions:
//  erf(), erfc(), and erex(), and calerf().
//
//  The calling statements for the primary entries are:
//
//         y = erf(x),
//
//         y = erfc(x),
//  and
//         y = erex(x).
//
//  The function calerf() is intended for internal packet use only,
//  all computations within the packet being concentrated in this
//  routine. The functions invoke calerf() with the
//  statement
//
//         calerf(arg,jint)
//
//  where the parameter usage is as follows
//
//      Function                     Parameters for calerf()
//       call              arg            jint
//
//     erf(arg)     ANY REAL ARGUMENT      0
//     erfc(arg)    fabs(arg)  < XBIG      1
//     erex(arg)    XNEG < arg < XMAX      2
//
//  The main computation evaluates near-minimax approximations
//  from "Rational Chebyshev approximations for the error function"
//  by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
//  transportable program uses rational functions that theoretically
//  approximate  erf(x) and  erfc(x) to at least 18 significant
//  decimal digits. The accuracy achieved depends on the arithmetic
//  system, the compiler, the intrinsic functions, and proper
//  selection of the machine-dependent constants.
//
//******************************************************************
//
//  Explanation of machine-dependent constants
//
//   XMIN   = the smallest positive floating-point number.
//   XINF   = the largest positive finite floating-point number.
//   XNEG   = the largest negative argument acceptable to erex();
//            the negative of the solution to the equation
//            2*exp(x*x) = XINF.
//   XSMALL = argument below which erf(x) may be represented by
//            2*x/sqrt(pi) and above which  x*x  will not underflow.
//            A conservative value is the largest machine number x
//            such that   1.0 + x = 1.0   to machine precision.
//   XBIG   = largest argument acceptable to erfc();  solution to
//            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
//            W(x) = exp(-x*x)/[x*sqrt(pi)].
//   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
//            machine precision.  A conservative value is
//            1/[2*sqrt(XSMALL)]
//   XMAX   = largest acceptable argument to erex(); the minimum
//            of XINF and 1/[sqrt(pi)*XMIN].
//
//  Approximate values for some important machines are:
//
//                          XMIN       XINF        XNEG     XSMALL
//
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
//  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
//  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
//  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
//  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
//
//
//                          XBIG       XHUGE       XMAX
//
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
//  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
//  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
//  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
//  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
//
//******************************************************************
//
//  Error returns
//
//  The program returns  erfc = 0      for  arg >= XBIG;
//
//                       erex = XINF   for  arg <  XNEG;
//      and
//                       erex = 0      for  arg >= XMAX.
//
//
//  Intrinsic functions required are:
//
//     fabs(), exp()
//
//-----------------------------------------------------------------
//  Based on the netlib FORTRAN package by W. J. Cody,
//  Mathematics and Computer Science Division
//  Argonne National Laboratory
//  Argonne, IL 60439
//
//  Latest modification of the above package: March 19, 1990
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//  Mathematical constants
//-----------------------------------------------------------------
const double ZERO    =  0.0e0;
const double HALF    =  0.5e0;
const double ONE     =  1.0e0;
const double TWO     =  2.0e0;
const double FOUR    =  4.0e0;
const double SIXTEEN = 16.0e0;

static const double SQRPI  = 5.6418958354775628695e-1;
static const double THRESH = 0.46875e0;

//-----------------------------------------------------------------
//  Machine-dependent constants (for IBM PC)
//-----------------------------------------------------------------
static const double XINF   =    1.79e308;
static const double XNEG   = -26.628e0;
static const double XSMALL =    1.11e-16;
static const double XBIG   =  26.543e0;
static const double XHUGE  =    6.71e7;
static const double XMAX   =    2.53e307;

//-----------------------------------------------------------------
//  Coefficients for approximation to  erf  in first interval
//-----------------------------------------------------------------
static const double A[5] = {
                           3.16112374387056560e00,
                           1.13864154151050156e02,
                           3.77485237685302021e02,
                           3.20937758913846947e03,
                           1.85777706184603153e-1
                           };

static const double B[4] = {
                           2.36012909523441209e01,
                           2.44024637934444173e02,
                           1.28261652607737228e03,
                           2.84423683343917062e03
                           };
//-----------------------------------------------------------------
//  Coefficients for approximation to  erfc  in second interval
//-----------------------------------------------------------------
static const double C[9] = {
                           5.64188496988670089e-1,
                           8.88314979438837594e0,
                           6.61191906371416295e01,
                           2.98635138197400131e02,
                           8.81952221241769090e02,
                           1.71204761263407058e03,
                           2.05107837782607147e03,
                           1.23033935479799725e03,
                           2.15311535474403846e-8
                           };

static const double D[8] = {
                           1.57449261107098347e01,
                           1.17693950891312499e02,
                           5.37181101862009858e02,
                           1.62138957456669019e03,
                           3.29079923573345963e03,
                           4.36261909014324716e03,
                           3.43936767414372164e03,
                           1.23033935480374942e03
                           };
//-----------------------------------------------------------------
//  Coefficients for approximation to  erfc  in third interval
//-----------------------------------------------------------------
static const double P[6] = {
                           3.05326634961232344e-1,
                           3.60344899949804439e-1,
                           1.25781726111229246e-1,
                           1.60837851487422766e-2,
                           6.58749161529837803e-4,
                           1.63153871373020978e-2
                           };

static const double Q[5] = {
                           2.56852019228982242e00,
                           1.87295284992346047e00,
                           5.27905102951428412e-1,
                           6.05183413124413191e-2,
                           2.33520497626869185e-3
                           };
//-----------------------------------------------------------------

register int i;
int iysq;
double del,x,xden,xnum,y,ysq;
double result;

x = arg;
y = fabs(x);
if(y <= THRESH)
  {
  //------------------------------------
  //  Evaluate  erf  for  |x| <= 0.46875
  //------------------------------------
  ysq = ZERO;
  if(y > XSMALL)ysq = y * y;
  xnum = A[4]*ysq;
  xden = ysq;
  for(i=0; i<3; i++)
     {
     xnum = (xnum + A[i])*ysq;
     xden = (xden + B[i])*ysq;
     }
  result = x * (xnum + A[3])/(xden + B[3]);
  if(jint != 0)result = ONE - result;
  if(jint == 2)result = exp(ysq)*result;

  return result;
  }
//-------------------------------------------
//  Evaluate  erfc  for 0.46875 <= |x| <= 4.0
//-------------------------------------------
else
if(y <= FOUR)
  {
  xnum = C[8]*y;
  xden = y;
  for(i=0; i<7; i++)
     {
     xnum = (xnum + C[i])*y;
     xden = (xden + D[i])*y;
     }
  result = (xnum + C[7])/(xden + D[7]);
  if(jint != 2)
    {
    iysq = y*SIXTEEN;
    ysq = ((double)iysq)/SIXTEEN;
    del = (y-ysq)*(y+ysq);
    result = exp(-ysq*ysq)*exp(-del)*result;
    }
  }
//-------------------------------
//  Evaluate  erfc  for |x| > 4.0
//-------------------------------
else{
    result = ZERO;
    if(y >= XBIG)
      {
      if((jint != 2) || (y >= XMAX))goto LAB300;
      if(y >= XHUGE)
        {
        result = SQRPI/y;
        goto LAB300;
        }
      }
    ysq = ONE/(y * y);
    xnum = P[5]*ysq;
    xden = ysq;
    for(i=0; i<4; i++)
       {
       xnum = (xnum + P[i])*ysq;
       xden = (xden + Q[i])*ysq;
       }
    result = ysq *(xnum + P[4])/(xden + Q[4]);
    result = (SQRPI - result)/y;
    if(jint != 2)
      {
      iysq = y*SIXTEEN;
      ysq = ((double)iysq)/SIXTEEN;
      del = (y-ysq)*(y+ysq);
      result = exp(-ysq*ysq)*exp(-del)*result;
      }
    }

//-----------------------------------------
//  Fix up for negative argument, erf, etc.
//-----------------------------------------
LAB300:
if(jint == 0)
  {
  result = (HALF - result)+HALF;
  if(x < ZERO)result = -result;
  }
else
if(jint == 1)
  {
  if(x < ZERO)result = TWO - result;
  }
else{ // jint == 2

    if(x < ZERO)
      {
      if(x < XNEG)result = XINF;
      else{
          iysq = x*SIXTEEN;
          ysq = ((double)iysq)/SIXTEEN;
          del = (x-ysq)*(x+ysq);
          y = exp(ysq*ysq)*exp(del);
          result = (y+y)-result;
          }
      }
    }

return result;
}





double erf(const double x)
{
//-------------------------------------------------------------------
//  This function computes approximate values for erf(x).
//  (see comments heading calerf()).
//
//  Based on the netlib package by W. J. Cody, January 8, 1985
//-------------------------------------------------------------------
return calerf(x,0);
}






double calerf_erfc( const double x )
{
//-------------------------------------------------------------------
//  This function computes approximate values for erfc(x).
//  (see comments heading calerf()).
//
//  Based on the netlib package by W. J. Cody, January 8, 1985
//-------------------------------------------------------------------
return calerf(x,1);
}





double erex(const double x)
{
//-----------------------------------------------------------------
//  This function computes approximate values for
//  exp(x*x) * erfc(x).
//  (see comments heading calerf()).
//
//  Based on the netlib package by W. J. Cody, March 30, 1987
//-----------------------------------------------------------------
return calerf(x,2);
}











