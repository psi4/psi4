/*! \file functional.cc
  \ingroup CINTS
  functionals go here

  The way the nomenclature workes is
  
  name_e just computes the energy
  name_ed computes both the energy and the potential
  name_ed2 includes second derivatives of the functional
  
  this is done for efficiency
  
  ----------------------------------------- */ 

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
#include"pade.h"

namespace psi { namespace cints {
  struct fun_info_s slater_e(struct den_info_s den_info){
    double Cx; /*= -0.930525736349100;*/
    struct fun_info_s exch_info;
    
    /***************************************************/
    Cx = -(9.0/4.0)*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));
    /***************************************************/

    exch_info.eval = 2.0*Cx*pow(den_info.den,4.0/3.0);
    return exch_info;
  }


  struct fun_info_s slater_ed(struct den_info_s den_info){
    double Cx /*= -0.930525736349100*/;
    double dCx /*= -1.240700981798800*/;
    struct fun_info_s exch_info;
    
    /**********************************************/
    Cx = -(9.0/4.0)*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));
    dCx = -3.0*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));
    /**********************************************/
    
    exch_info.eval = 2.0*Cx*pow(den_info.den,4.0/3.0);
    exch_info.dpval = dCx*pow(den_info.den,1.0/3.0);
    
    return exch_info;
  }


  /*struct density(struct den_info_s den_info){
    return den_info.den;
    }
  */

  struct fun_info_s no_funct(struct den_info_s den_info){
    struct fun_info_s fun_info;
    
    fun_info.eval = 0.0;
    fun_info.dval = 0.0;
    fun_info.dpval = 0.0;
    fun_info.dgval = 0.0;
    
    return fun_info;
  }


  struct fun_info_s VWN5_e(struct den_info_s den_info){
    double A = 0.0621814/2.0;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double Q = 6.1519908198;
    double threed4pi = 0.2387324146;
    
    double p;
    double ec;
    double x;
    double temp;
    
    struct fun_info_s corr_info;

    p = 2.0*den_info.den;
    temp = threed4pi/p;
    x = pow(temp,0.3333333333333);

    ec = Pade_int(x,x0,b,c,A,Q);

    corr_info.eval = p*ec;
    
    return corr_info;
  }


  struct fun_info_s VWN5_ed(struct den_info_s den_info){

    double A = 0.0621814/2.0;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double Q = 6.1519908198;
    double threed4pi = 0.2387324146;
    double x, temp1;
    double dxdp;
    double ec,dec;
    double ec_sum,dec_sum;
    double p;
    
    struct fun_info_s corr_info;

    p = 2.0*den_info.den;

    temp1 = threed4pi/p;
    x = pow(temp1,0.3333333333);
    dxdp = -x/(3.0*p);
    ec = Pade_int(x,x0,b,c,A,Q);;
    dec = d_Pade_int(x,x0,b,c,A,Q);
    
    corr_info.eval = p*ec;
    corr_info.dpval = ec+p*dxdp*dec;
    
   return corr_info;
  }

  /* This is the functional in which Gaussian uses */

  struct fun_info_s VWN4_e(struct den_info_s den_info){
    double A = 0.0621814/2.0;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;
    double Q = 0.04489988864;
    double threed4pi = 0.2387324146;
    
    double p;
    double ec;
    double x;
    double temp;
    
    struct fun_info_s corr_info;

    p = 2.0*den_info.den;
    temp = threed4pi/p;
    x = pow(temp,0.3333333333333);

    ec = Pade_int(x,x0,b,c,A,Q);

    corr_info.eval = p*ec;
    
    return corr_info;
  }

  struct fun_info_s VWN4_ed(struct den_info_s den_info){

    double A = 0.0621814/2.0;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;
    double Q = 0.04489988864;
    double threed4pi = 0.2387324146;
    double x, temp1;
    double dxdp;
    double ec,dec;
    double ec_sum,dec_sum;
    double p;
    
    struct fun_info_s corr_info;

    p = 2.0*den_info.den;

    temp1 = threed4pi/p;
    x = pow(temp1,0.3333333333);
    dxdp = -x/(3.0*p);
    ec = Pade_int(x,x0,b,c,A,Q);;
    dec = d_Pade_int(x,x0,b,c,A,Q);
    
    corr_info.eval = p*ec;
    corr_info.dpval = ec+p*dxdp*dec;
    
    return corr_info;
  }

  struct fun_info_s Becke88_e(struct den_info_s den_info){

    double b = 0.0042;
    double gamma;
    double x;
    double p,p43;
    double fourthirds = 4.0/3.0;
    double threehalves = 3.0/2.0;
    double threefourthspi = 3.0/(4.0*_pi);
    double onethirds = 1.0/3.0;
    double sxp1;
    double xpsxp1;
    double grad;
    double gx;
    double Cx;
    double arcsinhx;
    
    struct fun_info_s exch_info;

    p = den_info.den;
    grad = sqrt(den_info.gamma);
    
    Cx = -threehalves*pow(threefourthspi,onethirds);
    /*fprintf(outfile,"\nCx = %10.15lf",Cx);*/
    p43 = pow(p,fourthirds);
    /*fprintf(outfile,"\np43 = %10.15lf",p43);
    fprintf(outfile,"\ngrad = %10.10lf",grad);*/
    x = grad/p43;
    
    sxp1 = sqrt(x*x+1.0);
    xpsxp1 = x+sxp1;
    arcsinhx = log(xpsxp1);
    
    /*fprintf(outfile,"\nx = %10.15lf arcsinhx = %10.15lf",x,arcsinhx);*/
    gx = Cx - (b*x*x)/(1.0+6.0*b*x*arcsinhx);
    
    exch_info.eval = 2.0*p43*gx;
    
    return exch_info;
    
  }
  struct fun_info_s Becke88_ed(struct den_info_s den_info){

    double b = 0.0042;
    double gamma;
    double x,x2;
    double p,p13,p43;
    double fourthirds = 4.0/3.0;
    double threehalves = 3.0/2.0;
    double threefourthspi = 3.0/(4.0*_pi);
    double onethirds = 1.0/3.0;
    double sxp1;
    double xpsxp1;
    double xdsxp1;
    double bx,bbxx6;
    double bx6arcsin;
    double onepbx6arcsin;
    double grad;
    double gx,gpx;
    double Cx;
    double arcsinhx;
    
    struct fun_info_s exch_info;

    p = den_info.den;
    grad = sqrt(den_info.gamma);
    
    Cx = -threehalves*pow(threefourthspi,onethirds);
    /*fprintf(outfile,"\nCx = %10.15lf",Cx);*/
    p13 = pow(p,onethirds);
    p43 = p13*p;
    /*fprintf(outfile,"\np43 = %10.15lf",p43);
    fprintf(outfile,"\ngrad = %10.10lf",grad);*/
    x = grad/p43;
    bx = b*x;
    x2 = x*x;
    bbxx6 = bx*bx*6.0;
    sxp1 = sqrt(x2+1.0);
    xpsxp1 = x+sxp1;
    xdsxp1 = x/sxp1;
    arcsinhx = log(xpsxp1);
    bx6arcsin = 6.0*bx*arcsinhx;
    onepbx6arcsin = 1.0+bx6arcsin;
    
    /*fprintf(outfile,"\nx = %10.15lf arcsinhx = %10.15lf",x,arcsinhx);*/
    gx = Cx - (b*x2)/onepbx6arcsin;
    gpx = (bbxx6*(xdsxp1-arcsinhx)-2.0*bx)/(onepbx6arcsin*onepbx6arcsin);
    /*fprintf(outfile,"\ngpx = %10.10lf",gpx);*/
    exch_info.eval = 2.0*p43*gx;
    exch_info.dpval = fourthirds*p13*(gx-x*gpx);
    exch_info.dgval = 0.5*gpx/grad;
    fprintf(outfile,"\ndgval = %10.10lf",exch_info.dgval);
    return exch_info;
    
  }
}}
