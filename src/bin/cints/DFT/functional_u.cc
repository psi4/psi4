/*! \file functional_u.cc
    \ingroup CINTS
    functionals go here
    
    ----------------------------------------- */ 

#include <cmath>
#include <cstring>
#include <cstdio>
#include <memory.h>
#include <cstdlib>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>
#include<pthread.h>
#include<libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"physconst.h"
#include"pade.h"

namespace psi { namespace CINTS {
  
  struct fun_info_s slater_u_e(struct den_info_s den_info){
    double Cx = -0.930525736349100;
    struct fun_info_s exch_info;

    /*Cx = -(9.0/4.0)*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));*/    
    
    exch_info.eval = Cx*(pow(den_info.dena,4.0/3.0)
			 +pow(den_info.denb,4.0/3.0));
    return exch_info;
  }

  struct fun_info_s slater_u_ed(struct den_info_s den_info){
    double Cx = -0.930525736349100;
    double dCx = -1.240700981798800;
    struct fun_info_s exch_info;
  
    /*dCx = -3.0*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));*/
    
    exch_info.eval = Cx*(pow(den_info.dena,4.0/3.0)
			 +pow(den_info.denb,4.0/3.0));
    
    exch_info.dvala = dCx*pow(den_info.dena,1.0/3.0);
    exch_info.dvalb = dCx*pow(den_info.denb,1.0/3.0);
    
    return exch_info;
  }
          
  struct fun_info_s no_funct_u(struct den_info_s den_info){
    struct fun_info_s fun_info;
    
    fun_info.eval = 0.0;
    fun_info.dvala = 0.0;
    fun_info.dvalb = 0.0;
    fun_info.ddvala = 0.0;
    fun_info.ddvalb = 0.0;
    
    return fun_info;
  }



  /* This is the functional in which Gaussian uses */

  struct fun_info_s VWN4_u_e(struct den_info_s den_info){
    
    /* paramagnetic case*/
    double Ap = 0.0621814/2.0;
    double x0p = -0.409286;
    double bp = 13.0720;
    double cp = 42.7198;
    double Qp = 0.044899888641577;
    
    /* ferromagnetic case*/
    double Af = 0.0621814/4.0;
    double x0f = -0.743294;
    double bf = 20.1231;
    double cf = 101.578;
    double Qf = 1.171685277708971;
    
    /* spin stiffness*/ 
    double Aa = -0.016886863940390; /*-1.0/(6.0*_pi*_pi)*/
    double x0a = -0.228344;
    double ba = 1.06835;
    double ca = 11.4813;
    double Qa = 6.692072046645942;
    
    double p,pa,pb;
    double ea,ep,ef,beta;
    double ec;
    double eta,geta,eta2,eta4;
    double onepluseta;
    double onemineta;
    double d2f0 = 1.709920934161365;  /* 4/9*(2^1/3-1) */
    double nined8 = 9.0/8.0;
    double fourthirds = 4.0/3.0;
    double onethirds = 1.0/3.0;
    double x;
    double temp;
    
    struct fun_info_s corr_info;
    
    pa = den_info.dena;
    pb = den_info.denb;
    p = pa+pb;
    eta = (pa-pb)/p;
    eta4 = eta*eta*eta*eta;
    onepluseta = pow(1+eta,fourthirds);
    onemineta = pow(1-eta,fourthirds);
    geta = nined8*(onepluseta+onemineta-2.0);
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,onethirds);
    ep = Pade_int(x,x0p,bp,cp,Ap,Qp);
    ef = Pade_int(x,x0f,bf,cf,Af,Qf);
    ea = Pade_int(x,x0a,ba,ca,Aa,Qa);
    beta = ((d2f0*(ef-ep))/ea)-1.0;
    
    ec = ep+ea*geta*(1.0+(beta*eta4));
    corr_info.eval = p*ec;
    
    return corr_info;
    
  }

  struct fun_info_s VWN4_u_ed(struct den_info_s den_info){
    
    /* paramagnetic case*/
    double Ap = 0.0621814/2.0;
    double x0p = -0.409286;
    double bp = 13.0720;
    double cp = 42.7198;
    double Qp = 0.044899888641577;
    
    /* ferromagnetic case*/
    double Af = 0.0621814/4.0;
    double x0f = -0.743294;
    double bf = 20.1231;
    double cf = 101.578;
    double Qf = 1.171685277708971;
    
    /* spin stiffness*/ 
    double Aa = -0.016886863940390; /*-1.0/(6.0*_pi*_pi)*/
    double x0a = -0.228344;
    double ba = 1.06835;
    double ca = 11.4813;
    double Qa = 6.692072046645942;
    
    double p,pa,pb;
    double ea,ep,ef;
    double ec,dec,deca,decb;
    double beta,dbeta,betaplus;
    double dep,def,dea;
 
    double eta,detaa,detab,eta2,eta3,eta4;
    double geta,dgeta;
    double onepluseta,onemineta;
    double onepluseta2,onemineta2;
    double d2f0 = 1.709920934161365;  /* 4/9*(2^1/3-1)*/
    double nineeighths = 9.0/8.0;
    double fourthirds = 4.0/3.0;
    double onethirds = 1.0/3.0;
    double threehalves = 3.0/2.0;
    double x, temp;
    double dxdp;
    double term2,term3,term4a,term4b;
    
    struct fun_info_s corr_info;
    
    pa = den_info.dena;
    pb = den_info.denb;
    p = pa+pb;
    eta = (pa-pb)/p;
    detaa = (1-eta)/p;
    detab = -(1+eta)/p;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    eta4 = eta2*eta2;
    onepluseta = pow(1+eta,fourthirds);
    onemineta = pow(1-eta,fourthirds);
    geta = nineeighths*(onepluseta+onemineta-2.0);
    onepluseta2 = pow(1+eta,onethirds);
    onemineta2 = pow(1-eta,onethirds);
    dgeta = threehalves*(onepluseta2-onemineta2);
        
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,onethirds);
    dxdp = -x/(3.0*p);
    ep = Pade_int(x,x0p,bp,cp,Ap,Qp);
    ef = Pade_int(x,x0f,bf,cf,Af,Qf);
    ea = Pade_int(x,x0a,ba,ca,Aa,Qa);
    dep = d_Pade_int(x,x0p,bp,cp,Ap,Qp);
    def = d_Pade_int(x,x0f,bf,cf,Af,Qf);
    dea = d_Pade_int(x,x0a,ba,ca,Aa,Qa);
    /*fprintf(outfile,"\nep = %e ef = %e ea = %e",ep,ef,ea);
    fprintf(outfile,"\ndep = %e def = %e dea = %e",dep,def,dea);*/
    beta = ((d2f0*(ef-ep))/ea)-1.0;
    dbeta = (d2f0/ea)*(def-dep-((ef-ep)/ea)*dea);
    betaplus = 1.0+(beta*eta4);
    
    ec = ep+ea*geta*(1.0+(beta*eta4));
    
    corr_info.eval = p*ec;
    
    term2 = dea*geta*betaplus;
    term3 = ea*geta*dbeta*eta4;
    term4a = (dgeta*betaplus+4.0*geta*beta*eta3)*detaa;
    term4b = (dgeta*betaplus+4.0*geta*beta*eta3)*detab;
    
    dec = dxdp*(dep+term2+term3);
    deca = dec+ea*term4a;
    decb = dec+ea*term4b;
    corr_info.dvala = ec+p*deca;
    corr_info.dvalb = ec+p*decb;
    
    return corr_info; 
  }

  struct fun_info_s VWN5_u_e(struct den_info_s den_info){

    /* paramagnetic case*/
    double Ap = 0.0621814/2.0;
    double x0p = -0.10498;
    double bp = 3.72744;
    double cp = 12.9352;
    double Qp = 6.151990819759080;
    
    /* ferromagnetic case*/
    double Af = 0.0621814/4.0;
    double x0f = -0.32500;
    double bf = 7.06042;
    double cf = 18.0578;
    double Qf = 4.730926909560114;

    /* spin stiffness */
    double Aa = -0.016886863940390; /*-1.0/(6.0*_pi*_pi)*/
    double x0a = -0.00475840;
    double ba = 1.13107;
    double ca = 13.0045;
    double Qa = 7.123108917818118;

    double p,pa,pb;
    double ea,ep,ef,beta;
    double ec;
    double eta,geta,eta2,eta4;
    double onepluseta;
    double onemineta;
    double d2f0 = 1.709920934161365;  /* 4/9*(2^1/3-1) */
    double nined8 = 9.0/8.0;
    double fourthirds = 4.0/3.0;
    double onethirds = 1.0/3.0;
    double x;
    double temp;
    
    struct fun_info_s corr_info;
    
    pa = den_info.dena;
    pb = den_info.denb;
    p = pa+pb;
    eta = (pa-pb)/p;
    eta4 = eta*eta*eta*eta;
    onepluseta = pow(1+eta,fourthirds);
    onemineta = pow(1-eta,fourthirds);
    geta = nined8*(onepluseta+onemineta-2.0);
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,onethirds);
    ep = Pade_int(x,x0p,bp,cp,Ap,Qp);
    ef = Pade_int(x,x0f,bf,cf,Af,Qf);
    ea = Pade_int(x,x0a,ba,ca,Aa,Qa);
    beta = ((d2f0*(ef-ep))/ea)-1.0;
    
    ec = ep+ea*geta*(1.0+(beta*eta4));
    corr_info.eval = p*ec;
    
    return corr_info;
  }

  struct fun_info_s VWN5_u_ed(struct den_info_s den_info){
    
    /* paramagnetic case*/
    double Ap = 0.0621814/2.0;
    double x0p = -0.10498;
    double bp = 3.72744;
    double cp = 12.9352;
    double Qp = 6.151990819759080;
    
    /* ferromagnetic case*/
    double Af = 0.0621814/4.0;
    double x0f = -0.32500;
    double bf = 7.06042;
    double cf = 18.0578;
    double Qf = 4.730926909560114;

    /* spin stiffness */
    double Aa = -0.016886863940390; /*-1.0/(6.0*_pi*_pi)*/
    double x0a = -0.00475840;
    double ba = 1.13107;
    double ca = 13.0045;
    double Qa = 7.123108917818118;
        
    double p,pa,pb;
    double ea,ep,ef;
    double ec,dec,deca,decb;
    double beta,dbeta,betaplus;
    double dep,def,dea;
 
    double eta,detaa,detab,eta2,eta3,eta4;
    double geta,dgeta;
    double onepluseta,onemineta;
    double onepluseta2,onemineta2;
    double d2f0 = 1.709920934161365;  /* 4/9*(2^1/3-1)*/
    double nineeighths = 9.0/8.0;
    double fourthirds = 4.0/3.0;
    double onethirds = 1.0/3.0;
    double threehalves = 3.0/2.0;
    double x, temp;
    double dxdp;
    double term2,term3,term4a,term4b;
    
    struct fun_info_s corr_info;
   
    pa = den_info.dena;
    pb = den_info.denb;
    p = pa+pb;
    eta = (pa-pb)/p;
    detaa = (1-eta)/p;
    detab = -(1+eta)/p;
    eta2 = eta*eta;
    eta3 = eta2*eta;
    eta4 = eta2*eta2;
    onepluseta = pow(1+eta,fourthirds);
    onemineta = pow(1-eta,fourthirds);
    geta = nineeighths*(onepluseta+onemineta-2.0);
    onepluseta2 = pow(1+eta,onethirds);
    onemineta2 = pow(1-eta,onethirds);
    dgeta = threehalves*(onepluseta2-onemineta2);
        
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,0.333333333);
    dxdp = -x/(3.0*p);
    ep = Pade_int(x,x0p,bp,cp,Ap,Qp);
    ef = Pade_int(x,x0f,bf,cf,Af,Qf);
    ea = Pade_int(x,x0a,ba,ca,Aa,Qa);
    dep = d_Pade_int(x,x0p,bp,cp,Ap,Qp);
    def = d_Pade_int(x,x0f,bf,cf,Af,Qf);
    dea = d_Pade_int(x,x0a,ba,ca,Aa,Qa);
    /*fprintf(outfile,"\nep = %e ef = %e ea = %e",ep,ef,ea);
    fprintf(outfile,"\ndep = %e def = %e dea = %e",dep,def,dea);*/
    beta = ((d2f0*(ef-ep))/ea)-1.0;
    dbeta = (d2f0/ea)*(def-dep-((ef-ep)/ea)*dea);
    betaplus = 1.0+(beta*eta4);
    
    ec = ep+ea*geta*(1.0+(beta*eta4));
    
    corr_info.eval = p*ec;
    
    term2 = dea*geta*betaplus;
    term3 = ea*geta*dbeta*eta4;
    term4a = (dgeta*betaplus+4.0*geta*beta*eta3)*detaa;
    term4b = (dgeta*betaplus+4.0*geta*beta*eta3)*detab;
    
    dec = dxdp*(dep+term2+term3);
    deca = dec+ea*term4a;
    decb = dec+ea*term4b;
    corr_info.dvala = ec+p*deca;
    corr_info.dvalb = ec+p*decb;
    
    return corr_info;   
  }                         
};};
