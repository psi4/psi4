/*! \file quartet_data.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<memory.h>
#include<stdlib.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif

namespace psi { namespace CINTS {

/*!--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives
 --------------------------------------------------------------------------------*/
void quartet_data(prim_data *Data, double_array_t *fjt_table, double AB2, double CD2,
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
  double F[4*CINTS_MAX_AM];
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

  if (!am) {
#ifdef USE_TAYLOR_FM
    taylor_compute_fm(F,T,0);
    Data->F[0] = F[0]*coef1;
#else
    int_fjt(fjt_table,0,T);
    Data->F[0] = fjt_table->d[0]*coef1;
#endif
  }
  else {
    Data->oo2zn = 0.5*oozn;
    Data->pon = zeta*oozn;
    Data->oo2z = 0.5/zeta;
    Data->oo2n = 0.5/eta;
    W.x = (sp1->P[pi][pj][0]*zeta+sp2->P[pk][pl][0]*eta)*oozn;
    W.y = (sp1->P[pi][pj][1]*zeta+sp2->P[pk][pl][1]*eta)*oozn;
    W.z = (sp1->P[pi][pj][2]*zeta+sp2->P[pk][pl][2]*eta)*oozn;

    if(T < small_T){ 
      for(i=0; i<=am; i++) 
	Data->F[i] = oo2np1[i]*coef1;
    }
    else {
#ifdef USE_TAYLOR_FM
      taylor_compute_fm(F,T,am);
      for(i=0;i<=am;i++)
	Data->F[i] = F[i]*coef1;
#else
      int_fjt(fjt_table,am,T);
      for(i=0;i<=am;i++)
	Data->F[i] = fjt_table->d[i]*coef1;
#endif
    }

    /* PA */
    Data->U[0][0] = sp1->PA[pi][pj][0];
    Data->U[0][1] = sp1->PA[pi][pj][1];
    Data->U[0][2] = sp1->PA[pi][pj][2];
    /* QC */
    Data->U[2][0] = sp2->PA[pk][pl][0];
    Data->U[2][1] = sp2->PA[pk][pl][1];
    Data->U[2][2] = sp2->PA[pk][pl][2];
    /* WP */
    Data->U[4][0] = W.x - sp1->P[pi][pj][0];
    Data->U[4][1] = W.y - sp1->P[pi][pj][1];
    Data->U[4][2] = W.z - sp1->P[pi][pj][2];
    /* WQ */
    Data->U[5][0] = W.x - sp2->P[pk][pl][0];
    Data->U[5][1] = W.y - sp2->P[pk][pl][1];
    Data->U[5][2] = W.z - sp2->P[pk][pl][2];
  }

  return;
}


};};
