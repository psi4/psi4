/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
  
  $Id$
*/

//\file wavef.cc
//\brief The hydrogenic bound and continuum states
/************************************************************************
 * Here is a madness representation of the hydrogenic wave functions.
 * The bound states come from the Gnu Scientific Library. The unbound
 * states are generated with the confluent hypergeometric function which
 * uses gmp and mpfr for extended precision
 * 
 * Using: Gnu Scientific Library          http://www.gnu.org/software/gsl/
 *        GNU Multiple Precision library  http://gmplib.org/
 *        mpfr                            http://www.mpfr.org/
 * By:    Nick Vence
 ************************************************************************/


#include "wavef.h"
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <float.h>
#include <math.h>

using namespace madness;

//MPI printing macros
time_t before, after;
//#define END_TIMER_C(msg,cplx) tt=cpu_time()-tt; std::cout << "Timer: " << msg << "(" << real(cplx) << " + " << imag(cplx) << "I) took " << tt << " seconds" << std::endl
//#define END_TIMER(msg) tt=cpu_time()-tt; printf("timer: %24.24s    took%8.2f seconds\n", msg, tt)
//Defining static members
const double baseWF::PI = M_PI;
const complexd baseWF::I(0.0,1.0);

/**********************************************************************
 * Philk (angular momentum resolved scattering states)
 * See Freidrich, Theoretical Atomic Physics
 * Appendix 5
 **********************************************************************/
Phikl::Phikl(const double Z, const double k, const int l, const double cutoff) 
    :ScatteringWF(Z, cutoff) 
    ,l_(l) {}
Phikl::Phikl(World& world, const double Z, const double k, const int l, const double cutoff) 
    :ScatteringWF(world, Z, cutoff) 
  ,l_(l) {}
//NOT IMPLEMENTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complexd Phikl::operator()(const vector3D& rVec) const {
    if( fabs(rVec[0])<cutoff_ && fabs(rVec[1])<cutoff_ && fabs(rVec[2])<cutoff_ ) {
        return 0.01;
    } else {
        return 0.0;
    }
}
//NOT IMPLEMENTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///f11 determines when to use the series conhyp and when to use the asymptotic aForm
///This particular f11 belongs to PhiK. I haven't tested that for Phikl and the C++ 
///design needs to be modified to support it.
complexd Phikl::f11(double xx) const {
    complexd ZZ(0.0,-xx);
    //The transition point was done by finding the minimum difference between
    //conhyp(k,r) - aForm(k,r) for different k values
    if(xx <= 4.1*Z_*Z_/(k_*k_) + 21.6) return conhyp(-I*Z_/k_, one, ZZ);
    else return aForm(ZZ);  //??? Should I implemnt aForm in ScatteringWF or in PhiK and Phikl ???
}
double Phikl::getk() const  { return k_; }
complexd Phikl::setAA() { return l_ + 1.0 + I*Z_/k_;}
complexd Phikl::setBB() { return complexd(2.0*l_ + 2.0, 0.0); }

/**********************************************************************
 * PhiK (directional scattering states)
 * See Landau and Lifshitz Quantum Mechanics Volume 3
 * Third Edition Formula (136.9)
 **********************************************************************/
PhiK::PhiK(const double Z, const vector3D& kVec, const double cutoff)
    :ScatteringWF(Z, cutoff)
    ,kVec_(kVec) {}
PhiK::PhiK(World& world, const double Z, const vector3D& kVec, const double cutoff)
    :ScatteringWF(world, Z, cutoff)
    ,kVec_(kVec) {}
double PhiK::getk() const {
    return sqrt(kVec_[0]*kVec_[0] + kVec_[1]*kVec_[1] + kVec_[2]*kVec_[2]);
}
complexd PhiK::setAA() { return -I*Z_/k_; }
complexd PhiK::setBB() { return complexd(1.0, 0.0); }
complexd PhiK::operator()(const vector3D& rVec) const {
    if( fabs(rVec[0])<cutoff_ && fabs(rVec[1])<cutoff_ && fabs(rVec[2])<cutoff_ ) {
        double kDOTr =    kVec_[0]*rVec[0] + kVec_[1]*rVec[1] + kVec_[2]*rVec[2];
        double r     = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);
        return 0.0634936359342 //  = (2PI)^-(3/2)
               * expPIZ_2kXgamma1pIZ_k
               * exp(I*kDOTr)
               * fit1F1(k_*r + kDOTr);
    } else {
        return 0.0;
    }
}
///f11(double) determines when to use the series conhyp and when to use the asymptotic aForm
complexd PhiK::f11(double xx) const {
    complexd ZZ(0.0,-xx);
    //The cutoff was done by finding the r value for which the precision (error)
    //of aForm was approximately that of it's asymptotic value: (1e-15, 1e-14)
    
    //conhyp's precision is strictly less than 1e-15, but large values of r take
    //much longer to evaluate and conhyp is set to throw an insufficient precision
    //error (which we don't handle) if convergence isn't reached in 20000 terms.
    
    //Convergence is different for each value of Z
    switch( int(Z_) ) {
    case 1:
        if( xx > 40.3 +  0.412/(k_*k_) - 10.5*k_ + 3*k_*k_) return aForm(ZZ);
        else return conhyp(AA,BB,ZZ);
    case 2:
        if( xx > 49.3 + 1.32/(k_*k_) - 16.3*k_ + 3.7*k_*k_) return aForm(ZZ);
        else return conhyp(AA,BB,ZZ);
    case 3:
        if( xx > 45.9 + 2.93/(k_*k_) - 10.2*k_ + 2*k_*k_ ) return aForm(ZZ);
        else return conhyp(AA,BB,ZZ);
    }
    throw "PhiK::f11 Z != {1, 2, 3}";
    return  conhyp(AA,BB,ZZ);
}

/****************************************************
 * The Scattering Wave Function
 * An abstract base class for the anglar momentum 
 * resolved basis or the vector basis scattering states
 ****************************************************/
ScatteringWF::ScatteringWF(World& world, const double Z, const double cutoff) : Z_(Z), cutoff_(cutoff) {}
ScatteringWF::ScatteringWF(const double Z, const double cutoff) : Z_(Z), cutoff_(cutoff) {}
void ScatteringWF::Init(World& world) {
    one = complexd(1.0, 0.0);
    dx = 4e-3;   //Mesh spacing <- OPTIMIZE
    k_ = getk();
    AA = setAA();
    BB = setBB();
    gammaBB = gamma(BB);
    expmIPIAArGammaBBmAA = std::exp(-I*PI*AA)/gamma(BB-AA);
    rGammaAA = 1.0/gamma(AA);
    AAmBB = AA-BB;
    mAA = -AA;
    expPIZ_2kXgamma1pIZ_k = std::exp(PI*Z_/(2*k_))*gamma(1.0+I*Z_/k_);
/**********************************************************************
 * How far must we tabulate our 1F1 to cover the domain?
 * V^(1/3) gives us the length of the box
 * Our spherical function needs only one octant of the box    V^(1/3)/2
 * sqrt(3) allows us to reach the corner of the cube  sqrt(3)*V^(1/3)/2
 * kr + kDOTr brings along another factor of 2k     k*sqrt(3)*V^(1/3)
 **********************************************************************/
    domain = k_*sqrt(3)*pow(FunctionDefaults<NDIM>::get_cell_volume(),1.0/3.0);    
    n = floor(domain/dx +1);
    MemberFuncPtr p1F1(this); //this level of wrapping now seems redundant
    //World is needed for timing the length of the CubicInterpolationTable
    fit1F1 = CubicInterpolationTable<complexd>(world, 0.0, domain, n, p1F1);
}
/****************************************************************
 * The asymptotic form of the hypergeometric function given by
 * Abramowitz and Stegun 13.5.1
 * **************************************************************/
complexd ScatteringWF::aForm(complexd ZZ) const {
    //complexd cA1 = expPIZ_k; 
    complexd ZZPmAA = pow(ZZ,mAA);
    complexd cA  = ZZPmAA*expmIPIAArGammaBBmAA;
    complexd expZZ = exp(ZZ);
    complexd ZZPAAmBB = pow(ZZ, AAmBB);
    complexd cB  = expZZ*ZZPAAmBB*rGammaAA;
    complexd termA(0,0);
    complexd termB(0,0);
    const int maxTerms = 24;
    double nFact = 1.0;        // 0! = 1
    complexd  zr = 1.0/ZZ;     
    complexd         zrn(1.0,0.0);   //(1/z)^0
    complexd        mzrn(1.0,0.0);   //(-1/z)^0
    complexd      pochAA(1.0,0.0);   //Pochhammer is the counting up factorial (A)_0 = 1
    complexd poch1pAAmBB(1.0,0.0);
    complexd   pochBBmAA(1.0,0.0);
    complexd    poch1mAA(1.0,0.0);
    
    for(int n=0; n<=maxTerms; n++) {
        //Suming the n'th term
        complexd contribA = pochAA*poch1pAAmBB*mzrn/nFact;
        termA += contribA;
        complexd contribB = pochBBmAA*poch1mAA*zrn/nFact;
        termB += contribB;
        //Calculating the (n+1)'th term
        zrn         *=  zr;         //   z^-n
        mzrn        *= -zr;         // (-z)^-n
        nFact       *= n+1;
        pochAA      *=        AA + 1.0*n;  // (x)_n = x(x+1)(x+2)..(x+n-1)
        poch1pAAmBB *= 1.0+AA-BB + 1.0*n;
        pochBBmAA   *=   BB - AA + 1.0*n;
        poch1mAA    *=  1.0 - AA + 1.0*n;
    }
    return gammaBB*(cA*termA + cB*termB);
}

complexd ScatteringWF::gamma(double re, double im) {
    gsl_sf_result lnr;
    gsl_sf_result arg;
    int status = gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
    if(status != 0) throw "Error: gsl_sf_lngamma: " + status;
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}
complexd ScatteringWF::gamma(complexd AA) {
    gsl_sf_result lnr;
    gsl_sf_result arg;
    int status = gsl_sf_lngamma_complex_e(real(AA), imag(AA), &lnr, &arg);
    if(status != 0) throw "Error: gsl_sf_lngamma: " + status;
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}
// void testGamma(World& world) {
//     if(world.rank() == 0) std::cout << "Testing Gamma:================================================" << std::endl;
//     if(world.rank() == 0) std::cout << "gamma(3.0,0.0) = " << gamma(3.0,0.0) << std::endl;
//     if(world.rank() == 0) std::cout << "gamma(0.0,3.0) = " << gamma(0.0,3.0) << std::endl;
//     if(world.rank() == 0) std::cout << "gamma(3.0,1.0) = " << gamma(3.0,1.0) << std::endl;
//     if(world.rank() == 0) std::cout << "gamma(1.0,3.0) = " << gamma(1.0,3.0) << std::endl;
// }


/******************************************
 * BoundWF
 ******************************************/
BoundWF::BoundWF(double Z, int nn, int ll, int mm ) : Z(Z) {
    gsl_set_error_handler_off();
    if(nn < 1) {
        std::cerr << "Thou shalt not have negative n!" << std::endl;
	exit(1);
    }
    if(ll<0 || ll>=nn) {
    std::cerr << "n = " << nn << "\tl = " << ll << std::endl;
	std::cerr << "l has broken the quantum commandments!" << std::endl;
	exit(1);
    }
    if(abs(mm) > ll) {
    std::cerr << "n = " << nn << "\tl = " << ll << "\tm = " << mm << std::endl;
	std::cerr << "m out of bounds error!" << std::endl;
	exit(1);
    }
    n=nn;
    l=ll;
    m=mm;
}
complexd BoundWF::operator()(const vector3D& rVec) const {
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += rVec[i]*rVec[i]; }
    double r = sqrt(sum);
    double cosTH;
    if(r==0.0) 
        cosTH = 1.0;
    else
        cosTH = rVec[2]/r;
    gsl_sf_result Rnl;
    int status = gsl_sf_hydrogenicR_e(n, l, Z, r, &Rnl);
    //    gsl_set_error_handler(NULL);     //Turns on the default error handler
    if(status == GSL_EUNDRFLW) { return complexd(0,0); }
    else if(status != 0)            MADNESS_EXCEPTION("gsl_ERROR: ",status);
    if(m==0) { return complexd(Rnl.val * gsl_sf_legendre_sphPlm(l, m, cosTH), 0.0); }
    else {
	gsl_sf_result rPhi;
	gsl_sf_result phi ; 
	gsl_sf_rect_to_polar(rVec[0], rVec[1], &rPhi, &phi);
	return complexd(   Rnl.val 
                         * gsl_sf_legendre_sphPlm(l, abs(m), cosTH)
                         * exp(complexd(0,m*phi.val))
                       );
    }
}


/*****************************************
 *Exp[ I*(k.r) ]
 *****************************************/
Expikr::Expikr( const vector3D& kVec) : kVec(kVec) {
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += kVec[i]*kVec[i]; }
    k = sqrt(sum);
    costhK = kVec[2]/k;
}
complexd Expikr::operator()(const vector3D& rVec) const {
    double kDOTr = 0.0;
    for(int i=0; i<NDIM; i++) {
        kDOTr += kVec[i]*rVec[i];
    }
    return exp(I*kDOTr);
}

/*****************************************
 *Exp[ -r^2 ]
 *****************************************/
Gaussian::Gaussian(const double aa) {
    if( a>0 ) a = aa;
    else a = 1.0;
}
complexd Gaussian::operator()(const vector3D& rVec) const {
    double r2 = 0.0;
    for(int i=0; i<NDIM; i++) {
        r2 += rVec[i]*rVec[i];
    }
    return complexd(exp(-a*r2),0);
}


Yl0::Yl0( int l=0 ) : l_(l)  { }

double Yl0::operator()(const vector3D& r) const {
    double cosTH = r[2]/std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    return  gsl_sf_legendre_sphPlm(l_, 0, cosTH);
}
