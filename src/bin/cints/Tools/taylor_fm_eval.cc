/*! \file taylor_fm_eval.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

void free_Taylor_Fm_Eval();

static void compute_fm(double *F, double T, unsigned int l) {
  return;
}

/*------------------------------------------------------
  Initialize Taylor_Fm_Eval object (computes incomplete
  gamma function via Taylor interpolation)
 ------------------------------------------------------*/
void init_Taylor_Fm_Eval(unsigned int mmax, double epsilon)
{
  int i, m;
  int T_idx;
  double T, T_new;
  double egamma, func, dfuncdT;
  double term, sum, denom, rel_error;

/*  if (mmax == 0)
    throw std::domain_error("Called init_Taylor_Fm_Eval with argument of zero");*/

  Taylor_Fm_Eval.cutoff = epsilon;
  /*---------------------------------------
    We are doing Taylor interpolation with
    n=TAYLOR_ORDER terms here:
    error <= delT^n/(n+1)!
   ---------------------------------------*/
  Taylor_Fm_Eval.order_interp = TAYLOR_ORDER;
  Taylor_Fm_Eval.delT = 2.0*pow(Taylor_Fm_Eval.cutoff*fac[Taylor_Fm_Eval.order_interp+1],
			1.0/Taylor_Fm_Eval.order_interp);
  Taylor_Fm_Eval.max_m = mmax + Taylor_Fm_Eval.order_interp - 1;

  /*------------------------------------------------
     Check if Taylor_Fm_Eval has been initialized with
     the same mmax before:
     2) yes  - re-initialize again
     3) no - initialize
   ------------------------------------------------*/
  if (Taylor_Fm_Eval.grid != NULL || Taylor_Fm_Eval.T_crit != NULL) {
    free_Taylor_Fm_Eval();
  }
  
  Taylor_Fm_Eval.T_crit = init_array(Taylor_Fm_Eval.max_m + 1);   /*--- m=0 is included! ---*/
  Taylor_Fm_Eval.max_T = 0;
  /*--- Figure out T_crit for each m and put into the T_crit ---*/
  for(m=Taylor_Fm_Eval.max_m;m>=0;m--) {
    /*------------------------------------------
      Damped Newton-Raphson method to solve
      T^{m-0.5}*exp(-T) = epsilon*Gamma(m+0.5)
      The solution is the max T for which to do
      the interpolation
     ------------------------------------------*/
    T = -log(epsilon);
    egamma = epsilon*sqrt(M_PI)*df[2*m]/pow(2,m);
    T_new = T;
    do {
      const double damping_factor = 0.2;
      double deltaT, max_deltaT, sign_deltaT;
      T = T_new;
      /* f(T) = the difference between LHS and RHS of the equation above */
      func = pow(T,m-0.5) * exp(-T) - egamma;
      dfuncdT = ((m-0.5) * pow(T,m-1.5) - pow(T,m-0.5)) * exp(-T);
      /* f(T) has 2 roots and has a maximum in between. If f'(T) > 0 we are to the left of the hump. Make a big step to the right. */
      if (dfuncdT > 0.0)
	T_new *= 2.0;
      else {
	/* damp the step */
	deltaT = -func/dfuncdT;
	sign_deltaT = (deltaT > 0.0) ? 1.0 : -1.0;
	max_deltaT = damping_factor * T;
	if (fabs(deltaT) > max_deltaT)
	  deltaT = sign_deltaT * max_deltaT;
	T_new = T + deltaT;
      }
    } while (fabs(func/egamma) >= SOFT_ZERO);
    Taylor_Fm_Eval.T_crit[m] = T_new;
    T_idx = (int) floor(T_new/Taylor_Fm_Eval.delT);
    if (T_idx > Taylor_Fm_Eval.max_T)
      Taylor_Fm_Eval.max_T = T_idx;
  }

  /*-------------------------------------------------------
    Tabulate the gamma function from t=delT to T_crit[m]:
    1) include T=0 though the table is empty for T=0 since
       Fm(0) is simple to compute
    2) modified MacLaurin series converges fastest for
       the largest m -> use it to compute Fmmax(T)
       see JPC 94, 5564 (1990).
    3) then either use the series to compute the rest
       of the row or maybe use downward recursion
   -------------------------------------------------------*/
  Taylor_Fm_Eval.grid = block_matrix(Taylor_Fm_Eval.max_T+1,Taylor_Fm_Eval.max_m+1);
  /*--- do the mmax first ---*/
  for(m=0;m<=Taylor_Fm_Eval.max_m;m++)
  for(T_idx = Taylor_Fm_Eval.max_T;
      T_idx >= 0;
      T_idx--) {
    T = T_idx*Taylor_Fm_Eval.delT;
    denom = (m+0.5);
    term = 0.5*exp(-T)/denom;
    sum = term;
    do {
      denom += 1.0;
      term *= T/denom;
      sum += term;
      rel_error = term/sum;
    } while (rel_error >= Taylor_Fm_Eval.cutoff);

    Taylor_Fm_Eval.grid[T_idx][m] = sum;
  }

  /*--- final touch - set the pointer to the function to compute Fm(T) ---*/
  Taylor_Fm_Eval.compute_Fm = compute_fm;

  return;
}


void free_Taylor_Fm_Eval() {

  free(Taylor_Fm_Eval.T_crit);
  Taylor_Fm_Eval.T_crit = NULL;
  free_block(Taylor_Fm_Eval.grid);
  Taylor_Fm_Eval.grid = NULL;
  
  return;
}

void taylor_compute_fm(double *F, double T, unsigned int l) {

  int m;
  unsigned int T_ind;
  double T_crit, two_T, exp_mT, h, F_m, F_mp1;
  double *F_row;

#define STATIC_OO2NP1
#define STATIC_OON
#include "static.h"

  T_crit = Taylor_Fm_Eval.T_crit[l];
  two_T = 2.0*T;

  /*------------------------
    First compute Fl(T) ...
   ------------------------*/
  if (T > T_crit) {
    /*--- Asymptotic formula ---*/
    F[l] = df[2*l]*sqrt(M_PI/2)/pow(two_T,l+0.5);
  }
  else {
    /*--- Taylor interpolation ---*/
    T_ind = (unsigned int) floor(0.5+T/Taylor_Fm_Eval.delT);
    h = T_ind*Taylor_Fm_Eval.delT - T;
    F_row = Taylor_Fm_Eval.grid[T_ind] + l;
    F[l] =          F_row[0] +
	         h*(F_row[1] +
	  oon[2]*h*(F_row[2] +
	  oon[3]*h*(F_row[3] +
	  oon[4]*h*(F_row[4] +
	  oon[5]*h*(F_row[5])))));
  }

  /*------------------------------------
    And then do downward recursion in m
   ------------------------------------*/
  if (l > 0) {
    F_mp1 = F[l];
    exp_mT = exp(-T);
    for(m=l-1;m>=0;m--) {
      F_m = (exp_mT + two_T*F_mp1)*oo2np1[m];
      F[m] = F_m;
      F_mp1 = F_m;
    }
  }

  return;
}

  


};};
