/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Code to compute the sigma_3 part of sigma

    \sigma_3(Ia,Ib) = \sum_{Ja,Jb} \sum_{ijkl} 
                      <Jb|E^b_{ij}|Ib> <Ja|E^a_{kl}|Ia> (ij|kl) C(Ja,Jb)
*/

#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/detci/structs.h"

namespace psi {
namespace detci {

int form_ilist(struct stringwr *alplist, int Ja_list, int nas, int kl, int *L,
               int *R, double *Sgn);
int form_ilist_rotf(int *Cnt, int **Ridx, signed char **Sn, int **Ij, int nas,
                    int kl, int *L, int *R, double *Sgn);

#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))

/*
** S3_BLOCK_VDIAG()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
** 
** alplist and betlist refer to the alpha and beta lists for a particular
** alpha and beta codes (the combination of which specifies the current
** block of sigma)
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/
void s3_block_vdiag(struct stringwr *alplist, struct stringwr *betlist,
                    double **C, double **S, double *tei, int nas, int nbs,
                    int cnas, int Ib_list, int Ja_list, int Jb_list, int Ib_sym,
                    int Jb_sym, double **Cprime, double *F, double *V,
                    double *Sgn, int *L, int *R, int norbs, int *orbsym) {
  struct stringwr *Ia;
  unsigned int Ia_ex;
  int ij, i, j, t, kl, I, J, RJ;
  double tval, VS, *CprimeI0, *CI0;
  int jlen, Jacnt, *Iaij, Ia_idx;
  unsigned int *Iaridx;
  signed char *Iasgn;
  double *Tptr;
  int npthreads, rc, status;

  /* loop over i, j */
  for (i = 0; i < norbs; i++) {
    for (j = 0; j <= i; j++) {
      if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
      ij = ioff[i] + j;
      jlen = form_ilist(betlist, Jb_list, nbs, ij, L, R, Sgn);

      if (!jlen) continue;
      /*  outfile->Printf("S3_BLOCK_VDIAG: ij = %d\t jlen = %d\n", ij, jlen);
       */
      Tptr = tei + ioff[ij];

      /* gather operation */
      for (I = 0; I < cnas; I++) {
        CprimeI0 = Cprime[I];
        CI0 = C[I];
        for (J = 0; J < jlen; J++) {
          tval = Sgn[J];
          CprimeI0[J] = CI0[L[J]] * tval;
        }
      }

      for (Ia = alplist, Ia_idx = 0; Ia_idx < nas; Ia_idx++, Ia++) {
        /* loop over excitations E^a_{kl} from |A(I_a)> */
        Jacnt = Ia->cnt[Ja_list];
        Iaridx = Ia->ridx[Ja_list];
        Iasgn = Ia->sgn[Ja_list];
        Iaij = Ia->ij[Ja_list];

        zero_arr(V, jlen);
        for (Ia_ex = 0; Ia_ex < Jacnt && (kl = *Iaij++) <= ij; Ia_ex++) {
          I = *Iaridx++;
          tval = *Iasgn++;
          if (ij == kl) tval *= 0.5;
          VS = Tptr[kl] * tval;
          CprimeI0 = Cprime[I];

#ifdef USE_BS
          C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
          for (J = 0; J < jlen; J++) {
            V[J] += VS * CprimeI0[J];
          }
#endif
        }

        /* scatter */
        for (J = 0; J < jlen; J++) {
          RJ = R[J];
          S[Ia_idx][RJ] += V[J];
        }

      } /* end loop over Ia */

    } /* end loop over j */
  }   /* end loop over i */
}

/*
** S3_BLOCK_V()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For non-diagonal blocks of s3
**
*/
void s3_block_v(struct stringwr *alplist, struct stringwr *betlist, double **C,
                double **S, double *tei, int nas, int nbs, int cnas,
                int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
                double **Cprime, double *F, double *V, double *Sgn, int *L,
                int *R, int norbs, int *orbsym) {
  struct stringwr *Ia;
  unsigned int Ia_ex;
  int ij, i, j, kl, ijkl, I, J, RJ;
  double tval, VS, *CprimeI0, *CI0;
  int jlen, Ia_idx, Jacnt, *Iaij;
  unsigned int *Iaridx;
  signed char *Iasgn;
  double *Tptr;

  /* loop over i, j */
  for (i = 0; i < norbs; i++) {
    for (j = 0; j <= i; j++) {
      if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
      ij = ioff[i] + j;
      jlen = form_ilist(betlist, Jb_list, nbs, ij, L, R, Sgn);

      if (!jlen) continue;

      Tptr = tei + ioff[ij];

      /* gather operation */
      for (I = 0; I < cnas; I++) {
        CprimeI0 = Cprime[I];
        CI0 = C[I];
        for (J = 0; J < jlen; J++) {
          tval = Sgn[J];
          CprimeI0[J] = CI0[L[J]] * tval;
        }
      }

      timer_on("CIWave: s3_mt");
      for (Ia = alplist, Ia_idx = 0; Ia_idx < nas; Ia_idx++, Ia++) {
        /* loop over excitations E^a_{kl} from |A(I_a)> */
        Jacnt = Ia->cnt[Ja_list];
        Iaridx = Ia->ridx[Ja_list];
        Iasgn = Ia->sgn[Ja_list];
        Iaij = Ia->ij[Ja_list];

        zero_arr(V, jlen);

        for (Ia_ex = 0; Ia_ex < Jacnt; Ia_ex++) {
          kl = *Iaij++;
          I = *Iaridx++;
          tval = *Iasgn++;
          ijkl = INDEX(ij, kl);
          VS = tval * tei[ijkl];
          CprimeI0 = Cprime[I];

#ifdef UBLAS
          C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
          for (J = 0; J < jlen; J++) {
            V[J] += VS * CprimeI0[J];
          }
#endif
        }

        /* scatter */
        for (J = 0; J < jlen; J++) {
          RJ = R[J];
          S[Ia_idx][RJ] += V[J];
        }

      } /* end loop over Ia */
      timer_off("CIWave: s3_mt");

    } /* end loop over j */
  }   /* end loop over i */
}

int form_ilist(struct stringwr *alplist, int Ja_list, int nas, int kl, int *L,
               int *R, double *Sgn) {
  int inum = 0, Ia_idx, Ia_ex, Iacnt, ij;
  int *Iaij;
  struct stringwr *Ia;
  unsigned int *Iaridx;
  signed char *Iasgn;

  /* loop over Ia */
  for (Ia = alplist, Ia_idx = 0; Ia_idx < nas; Ia_idx++, Ia++) {
    /* loop over excitations E^a_{kl} from |A(I_a)> */

    Iacnt = Ia->cnt[Ja_list];
    if (!Iacnt) continue;
    Iaridx = Ia->ridx[Ja_list];
    Iasgn = Ia->sgn[Ja_list];
    Iaij = Ia->ij[Ja_list];
    Ia_ex = 0;
    while (Ia_ex < Iacnt && (ij = *Iaij++) < kl) Ia_ex++;
    if (ij == kl) {
      *R++ = Ia_idx;
      *L++ = Iaridx[Ia_ex];
      *Sgn++ = (double)Iasgn[Ia_ex];
      inum++;
    }
  } /* end loop over Ia */
    /*  if(inum) {
          outfile->Printf("form_ilist: nas = %d\n", nas);
          outfile->Printf("form_ilist: jlen = %d\n", inum);

        }
    */
  return (inum);
}

/*
** S3_BLOCK_VDIAG_ROTF()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/

void s3_block_vdiag_rotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
                         signed char **Sn[2], double **C, double **S,
                         double *tei, int nas, int nbs, int cnas, int Ib_list,
                         int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
                         double **Cprime, double *F, double *V, double *Sgn,
                         int *L, int *R, int norbs, int *orbsym) {
  int Ia_ex;
  int ij, i, j, kl, I, J, RJ;
  double tval, VS, *CprimeI0, *CI0;
  int jlen, Ia_idx, Jacnt, *Iaij;
  int *Iaridx;
  signed char *Iasgn;
  double *Tptr;

  /* loop over i, j */
  for (i = 0; i < norbs; i++) {
    for (j = 0; j <= i; j++) {
      if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
      ij = ioff[i] + j;
      jlen = form_ilist_rotf(Cnt[1], Ridx[1], Sn[1], Ij[1], nbs, ij, L, R, Sgn);

      if (!jlen) continue;

      Tptr = tei + ioff[ij];

      /* gather operation */
      for (I = 0; I < cnas; I++) {
        CprimeI0 = Cprime[I];
        CI0 = C[I];
        for (J = 0; J < jlen; J++) {
          tval = Sgn[J];
          CprimeI0[J] = CI0[L[J]] * tval;
        }
      }

      /* loop over Ia */
      for (Ia_idx = 0; Ia_idx < nas; Ia_idx++) {
        /* loop over excitations E^a_{kl} from |A(I_a)> */
        Jacnt = Cnt[0][Ia_idx];
        Iaridx = Ridx[0][Ia_idx];
        Iasgn = Sn[0][Ia_idx];
        Iaij = Ij[0][Ia_idx];

        zero_arr(V, jlen);

        /* rotf doesn't yet ensure kl's in order */
        for (Ia_ex = 0; Ia_ex < Jacnt; Ia_ex++) {
          kl = *Iaij++;
          I = *Iaridx++;
          tval = *Iasgn++;
          if (kl > ij) continue;
          if (ij == kl) tval *= 0.5;
          VS = Tptr[kl] * tval;
          CprimeI0 = Cprime[I];

#ifdef USE_BLAS
          C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
          for (J = 0; J < jlen; J++) {
            V[J] += VS * CprimeI0[J];
          }
#endif
        }

        /* scatter */
        for (J = 0; J < jlen; J++) {
          RJ = R[J];
          S[Ia_idx][RJ] += V[J];
        }

      } /* end loop over Ia */

    } /* end loop over j */
  }   /* end loop over i */
}

/*
** S3_BLOCK_VROTF()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For non-diagonal blocks of s3
**
*/
void s3_block_vrotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
                    signed char **Sn[2], double **C, double **S, double *tei,
                    int nas, int nbs, int cnas, int Ib_list, int Ja_list,
                    int Jb_list, int Ib_sym, int Jb_sym, double **Cprime,
                    double *F, double *V, double *Sgn, int *L, int *R,
                    int norbs, int *orbsym) {
  int Ia_ex;
  int ij, i, j, kl, ijkl, I, J, RJ;
  double tval, VS, *CprimeI0, *CI0;
  int jlen, Ia_idx, Jacnt, *Iaij;
  int *Iaridx;
  signed char *Iasgn;
  double *Tptr;

  /* loop over i, j */
  for (i = 0; i < norbs; i++) {
    for (j = 0; j <= i; j++) {
      if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
      ij = ioff[i] + j;
      jlen = form_ilist_rotf(Cnt[1], Ridx[1], Sn[1], Ij[1], nbs, ij, L, R, Sgn);

      if (!jlen) continue;

      Tptr = tei + ioff[ij];

      /* gather operation */
      for (I = 0; I < cnas; I++) {
        CprimeI0 = Cprime[I];
        CI0 = C[I];
        for (J = 0; J < jlen; J++) {
          tval = Sgn[J];
          CprimeI0[J] = CI0[L[J]] * tval;
        }
      }

      /* loop over Ia */
      for (Ia_idx = 0; Ia_idx < nas; Ia_idx++) {
        /* loop over excitations E^a_{kl} from |A(I_a)> */
        Jacnt = Cnt[0][Ia_idx];
        Iaridx = Ridx[0][Ia_idx];
        Iasgn = Sn[0][Ia_idx];
        Iaij = Ij[0][Ia_idx];

        zero_arr(V, jlen);

        for (Ia_ex = 0; Ia_ex < Jacnt; Ia_ex++) {
          kl = *Iaij++;
          I = *Iaridx++;
          tval = *Iasgn++;
          ijkl = INDEX(ij, kl);
          VS = tval * tei[ijkl];
          CprimeI0 = Cprime[I];

#ifdef USE_BLAS
          C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
          for (J = 0; J < jlen; J++) {
            V[J] += VS * CprimeI0[J];
          }
#endif
        }

        /* scatter */
        for (J = 0; J < jlen; J++) {
          RJ = R[J];
          S[Ia_idx][RJ] += V[J];
        }

      } /* end loop over Ia */

    } /* end loop over j */
  }   /* end loop over i */
}

int form_ilist_rotf(int *Cnt, int **Ridx, signed char **Sn, int **Ij, int nas,
                    int kl, int *L, int *R, double *Sgn) {
  int inum = 0, Ia_idx, Ia_ex, Iacnt, ij;
  int *Iaij;
  int *Iaridx;
  signed char *Iasgn;

  /* loop over Ia */
  for (Ia_idx = 0; Ia_idx < nas; Ia_idx++) {
    /* loop over excitations E^a_{kl} from |A(I_a)> */

    Iacnt = Cnt[Ia_idx];
    if (!Iacnt) continue;
    Iaridx = Ridx[Ia_idx];
    Iasgn = Sn[Ia_idx];
    Iaij = Ij[Ia_idx];
    Ia_ex = 0;
    for (Ia_ex = 0; Ia_ex < Iacnt; Ia_ex++) {
      ij = *Iaij++;
      if (ij == kl) {
        *R++ = Ia_idx;
        *L++ = Iaridx[Ia_ex];
        *Sgn++ = (double)Iasgn[Ia_ex];
        inum++;
      }
    }
  } /* end loop over Ia */

  return (inum);
}
}}  // namespace psi::detci
