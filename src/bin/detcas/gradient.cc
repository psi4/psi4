/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

/*
** calc_grad_1
**
** Calculates the orbital gradient from the MO Lagrangian read in
** from the CLAG program
**
** An alternative method is available using calc_grad_2()
**/
void calc_grad_1(int npair, int *ppair, int *qpair, double **lag, double *grad)
{
  int pair, p, q;
  double value;

  /* this def of the gradient is consistent with orbital rotations defined
   * in the same sense as our VBD paper...also, what we would call the
   * Lagrangian in the VBD paper is actually twice the Lagrangian computed
   * by CLAG, so double the contribution.
   */
  for (pair=0; pair<npair; pair++) {
    p = *ppair++;
    q = *qpair++;
    value = 2.0 * (lag[q][p] - lag[p][q]);
    grad[pair] = value;
  }

}




/*
** calc_grad_2
**
** Calculates the orbital gradient from the F_core and F_act quantities,
** not from CLAG
**
** I am assuming that the pairs (p,q) are always given such that p>=q
**
*/
void calc_grad_2(int npairs, int *ppair, int *qpair, double *F_core, 
                 double *tei, double **opdm, double *tpdm, double *F_act, 
                 int firstact, int lastact, double *grad)
{

  int pair, p, q, pq, pp, qq;
  int i,ii,a,aa,t,tt,ti,u,tu,au,iu,v,w,vw,tuvw,auvw,iuvw;
  double value;

  /* this def of the gradient is consistent with orbital rotations defined
   * in the same sense as our VBD paper...we have to reverse the sign
   * on all their terms
   */

  /* loop over the independent pairs */
  for (pair=0; pair<npairs; pair++) {
    p = ppair[pair];
    q = qpair[pair];
    pq = ioff[p] + q;
    pp = ioff[p] + p;
    qq = ioff[q] + q;
  
    /* g_{ai}, i.e., inactive virt/inactive occ */
    if (p >= lastact && q < firstact) {
      grad[pair] = -4.0 * (F_core[pq] + F_act[pq]);
    }

    /* g_{at}, i.e., inactive virt with active orb */
    else if (p >= lastact && q >= firstact) {
      a = p;  t = q;
      aa = ioff[a] + a;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        au = INDEX(a,u);
        value += opdm[t][u] * F_core[au];
      }
      grad[pair] = -2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        au = INDEX(a,u); 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            auvw = INDEX(au,vw);
            value += tpdm[tuvw] * tei[auvw];
          }
        }
      }
      grad[pair] -= 2.0 * value; 
       
    } 

    /* g_{ti}, i.e., active orb with inactive occ */
    else if (p >= firstact && q < firstact) {
      t = p;  i = q;
      tt = ioff[t] + t;
      ii = ioff[i] + i; 
      ti = ioff[t] + i;
      
      grad[pair] = -4.0 * (F_core[ti] + F_act[ti]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        iu = INDEX(i,u);
        value += opdm[t][u] * F_core[iu];
      }
      grad[pair] += 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        iu = INDEX(i,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            iuvw = INDEX(iu,vw);
            value += tpdm[tuvw] * tei[iuvw];
          }
        }
         
      }
      grad[pair] += 2.0 * value; 
    }

    else {
      fprintf(outfile, 
             "(calc_grad_2): Error, unrecognized class of indep pair\n");
    }

  } /* end loop over pairs */

}

}} // end namespace psi::detcas

