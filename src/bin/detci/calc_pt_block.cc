/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#define EXTERN

#include <cstdio>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

extern int *ioff;
extern struct stringwr **alplist;
extern struct stringwr **betlist;
                                                                                
                                                                                
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

/*
** calc_pt_block(): Function calculates the 2nd order perturbation theory
** correction to a block of the CI vector assuming all corrections 
** involve excitations to the external orbitals only
**
** C. David Sherrill
** Feburary 2004
*/
void calc_pt_block(struct stringwr *alplist_local, struct stringwr
  *betlist_local, double **C, int nas, int nbs)
{
  int na;                /* number of alpha electrons */
  int nb;                /* number of beta electrons */
  int acnt, bcnt         /* row (acnt) and column (bcnt) of block */
  int first_external;    /* first external orbital for PT */ 
  int last_external;     /* last external orbital for PT */
  int a1, a2, b1, b2;    /* string indices */
  int i, j;              /* internal orbitals */
  int A, B;              /* external orbitals */
  int Ai, Aj, Bi, Bj;    /* integral lookup indices */
  int AiBj, AjBi;        /* more integral lookup indices */
  struct stringwr *betlist0; /* to reset the beta list for every alpha row */
  double Cval;           /* CI coefficient */

  betlist0 = betlist_local;

  /* loop over all the elements of the block */
  for (acnt=0; acnt<nas; acnt++) {
    for (bcnt=0,betlist_local=bestlist0; bcnt<nbs; bcnt++) {

      Cval = C[acnt][bcnt];

      /* here we should compute the denominator using code like that
         in calc_hd_block(), but need to replace one or two of the orbitals
         with external ones A and B, which is done in loops below...
         loops need to be coupled.  I can see that orbital energies
         make the denominators easier to compute... */

      /* double alpha excitations */
      for (a1=0; a1<na; a1++) {
        i = (int) alplist_local->occs[a1];
        for (a2=0; a2<a1; a2++) {
          j = (int) alplist_local->occs[a2];
          /* i > j */
          for (A=first_external; A<=last_external; A++) {
            Ai = ioff[A]+i;
            Aj = ioff[A]+j;
            for (B=first_external; B<A; B++) {
              /* A > B */
              /* matrix element <AB||ij>, will square so sign no matter */
              /* <AB||ij> = <AB|ij>-<AB|ji> = (Ai|Bj)-(Aj|Bi) */ 
              Bj = ioff[B]+j;
              Bi = ioff[B]+i;
              AiBj = ioff[Ai]+Bj;
              AjBi = ioff[Aj]+Bi;
              value = tei[AiBj] - tei[AjBi];
              value *= value;
              value *= Cval;
              value /= denom;
            } /* end loop over B */
          } /* end loop over A */
        } /* end loop over j */
      } /* end loop over i */

    betlist_local++;
    } /* end loop over beta list */

  alplist_local++;
  } /* end loop over alpha list */

}

}} // namespace psi::detci

