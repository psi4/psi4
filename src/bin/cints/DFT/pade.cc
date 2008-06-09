/*! \file pade.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <memory.h>
#include <cstdlib>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"physconst.h"


namespace psi { namespace CINTS {

/*! Implementation of the VWN functional extrapolation
  functions as given in Can J. Phys., 58, 1980
   pp. 1200 */
double Pade_int(double p, double x0, double b, double c, double A, double Q){

    double x,X,X0,invQ;
    double txpb,Qd2xpb;
    double term1,term2,mult1,term3,term4;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7;
    double ec;
  
    x = sqrt(p);
    txpb = 2.0*x+b;
    Qd2xpb = Q/txpb;
    invQ = 1.0/Q;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);

    temp1 = x*x*X;
    term1 = log(temp1);

    temp3 = atan(Qd2xpb);
    term2 = 2.0*b*invQ*temp3;

    mult1 = b*x0*X0;

    temp4 = (x-x0)*(x-x0)*X;
    term3 = log(temp4);

    temp5 = 2.0*(2.0*x0+b)*invQ;
    term4 = temp5*atan(Qd2xpb);

    ec = A*(term1+term2-mult1*(term3+term4));
    return ec;

}

double d_Pade_int(double p, double x0, double b, double c, double A, double Q){
    
    double x,Q2,X,X0,invQ;
    double txpb;
    double Qsptxpbs;
    double term1,term2,term3,term4,term5,term6;
    double mult1;
    double temp1,temp2,temp3,temp4,temp5;
    double dec;

    x=sqrt(p);
    A=A/(2*x);
    Q2 = Q*Q;
    txpb = 2.0*x+b;
    Qsptxpbs = txpb*txpb+Q2;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);
    
    term1 = 2.0/x;
   
    term2 = -txpb*X;

    term3 = -4.0*b/Qsptxpbs;

    mult1 = -b*x0*X0;

    temp3 = x-x0;
    term4 = 2.0/temp3;

    temp4 = 2.0*x0+b;
    term5 = -4.0*temp4/Qsptxpbs;

    dec = A*(term1+term2+term3+mult1*(term4+term2+term5));
    return dec;
}                          
};}
