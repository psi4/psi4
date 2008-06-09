/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstring>
#include<cstdio>
#include<cstdlib>
#include<memory.h>
#include<libint/libint.h>
#include<libderiv/libderiv.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif

namespace psi {
  namespace CINTS {
/*!--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives
 --------------------------------------------------------------------------------*/
void deriv1_quartet_data(prim_data *Data, double_array_t *fjt_table, double AB2, double CD2,
			 struct shell_pair *sp1, struct shell_pair *sp2, 
			 int am, int pi, int pj, int pk, int pl, double scale)
{
#define STATIC_OO2NP1
#include "static.h"

  /*----------------
    Local variables
   ----------------*/
  struct coordinates PQ, W;
#ifdef USE_TAYLOR_FM
  double F[2*CINTS_MAX_AM+1];
#endif
  int i;
  double small_T = UserOptions.cutoff;       /*--- Use only one term in Taylor expansion of Fj(T) if T < small_T ---*/
  double T;
  double coef1;
  double PQ2;
  double oozn;
  double zeta, eta, rho;

  zeta = sp1->gamma[pi][pj];
  eta = sp2->gamma[pk][pl];
  oozn = 1.0/(zeta+eta);
  Data->twozeta_a = 2.0*sp1->a1[pi];
  Data->twozeta_b = 2.0*sp1->a2[pj];
  Data->twozeta_c = 2.0*sp2->a1[pk];
  Data->twozeta_d = 2.0*sp2->a2[pl];
  Data->poz = eta*oozn;
  rho = zeta*Data->poz;
  coef1 = 2.0*sqrt(rho*M_1_PI)*scale*sp1->Sovlp[pi][pj]*sp2->Sovlp[pk][pl];
  PQ.x = sp1->P[pi][pj][0] - sp2->P[pk][pl][0];
  PQ.y = sp1->P[pi][pj][1] - sp2->P[pk][pl][1];
  PQ.z = sp1->P[pi][pj][2] - sp2->P[pk][pl][2];
  PQ2 = PQ.x*PQ.x;
  PQ2 += PQ.y*PQ.y;
  PQ2 += PQ.z*PQ.z;
  T = rho*PQ2;

  Data->oo2zn = 0.5*oozn;
  Data->pon = zeta*oozn;
  Data->oo2z = 0.5/zeta;
  Data->oo2n = 0.5/eta;
  W.x = (sp1->P[pi][pj][0]*zeta+sp2->P[pk][pl][0]*eta)*oozn;
  W.y = (sp1->P[pi][pj][1]*zeta+sp2->P[pk][pl][1]*eta)*oozn;
  W.z = (sp1->P[pi][pj][2]*zeta+sp2->P[pk][pl][2]*eta)*oozn;

  if(fabs(PQ2)<ZERO){ 
    for(i=0; i<=am+DERIV_LVL; i++) 
      Data->F[i] = oo2np1[i]*coef1;
    }
  else {
#ifdef USE_TAYLOR_FM
      taylor_compute_fm(F,T,am+DERIV_LVL);
      for(i=0;i<=am+DERIV_LVL;i++)
	Data->F[i] = F[i]*coef1;
#else
      int_fjt(fjt_table,am+DERIV_LVL,T);
      for(i=0;i<=am+DERIV_LVL;i++)
	Data->F[i] = fjt_table->d[i]*coef1;
#endif
    }

  /* PA */
  Data->U[0][0] = sp1->PA[pi][pj][0];
  Data->U[0][1] = sp1->PA[pi][pj][1];
  Data->U[0][2] = sp1->PA[pi][pj][2];
  /* PB */
  Data->U[1][0] = sp1->PB[pi][pj][0];
  Data->U[1][1] = sp1->PB[pi][pj][1];
  Data->U[1][2] = sp1->PB[pi][pj][2];
  /* QC */
  Data->U[2][0] = sp2->PA[pk][pl][0];
  Data->U[2][1] = sp2->PA[pk][pl][1];
  Data->U[2][2] = sp2->PA[pk][pl][2];
  /* QD */
  Data->U[3][0] = sp2->PB[pk][pl][0];
  Data->U[3][1] = sp2->PB[pk][pl][1];
  Data->U[3][2] = sp2->PB[pk][pl][2];
  /* WP */
  Data->U[4][0] = W.x - sp1->P[pi][pj][0];
  Data->U[4][1] = W.y - sp1->P[pi][pj][1];
  Data->U[4][2] = W.z - sp1->P[pi][pj][2];
  /* WQ */
  Data->U[5][0] = W.x - sp2->P[pk][pl][0];
  Data->U[5][1] = W.y - sp2->P[pk][pl][1];
  Data->U[5][2] = W.z - sp2->P[pk][pl][2];

  return;
}
};};
