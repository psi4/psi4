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
    \brief Enter brief description of file here
*/

#include <cstdio>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/detci/structs.h"

namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

/*
** S3_BLOCK_BZ1()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** Try routine from Bendazzoli and Evangelisti
**
** Pass in Ialist, Iblist, Jalist, Jblist.  Need to convert these
** to irreps for now; later we can make even more general.
**
*/

//void s3_block_bz1(int Ialist, int Iblist, int Jalist, int Jblist,
//      int nas, int nbs,
//      int cnbs, double *tei, double **C, double **S,
//      double **Cprime, double **Sprime)
//{
//
//   int Iasym, Jasym, Ibsym, Jbsym;
//   int norbs, *orbsym;
//   int i, j, k, l, ij, fullij, fullji, kl, fullkl, fulllk, tmpi, tmpj;
//   int I, J, I1, I2, J1, J2, S2;   /* also try S2 as double */
//   int ilen, jlen;
//   double V, VS;
//   double *Tptr;
//   int *OVptr, *OVptr2;
//   int signmask, nsignmask;
//   double *SprimeI0, *SprimeI1, *SprimeI2, *SprimeI3;
//   double *SprimeI4, *SprimeI5, *SprimeI6, *SprimeI7;
//   double *CprimeI0, *CprimeI1, *CprimeI2, *CprimeI3;
//   double *CprimeI4, *CprimeI5, *CprimeI6, *CprimeI7;
//
//   orbsym = CalcInfo.orbsym + CalcInfo.num_drc_orbs;
//   Iasym = Ialist;
//   Jasym = Jalist;
//   Ibsym = Iblist;
//   Jbsym = Jblist;
//
//   norbs = CalcInfo.num_ci_orbs;
//
//   signmask = 1 << (sizeof(int)*8-1);
//   nsignmask = ~signmask;
//
//   /* loop over k,l */
//   for (k=0; k<norbs; k++) {
//      for (l=0; l<norbs; l++) {
//         if ((orbsym[k] ^ orbsym[l] ^ Jasym ^ Iasym) != 0) continue;
//
//         kl = INDEX(k,l);
//         Tptr = tei + ioff[kl];
//         fullkl = k * norbs + l;
//         fulllk = l * norbs + k;
//
//         ilen = OV[Jalist][fulllk][0];
//
//         if (ilen == 0) continue;
//
//         for (I=0,OVptr=OV[Jalist][fulllk]+1; I<ilen; I++) {
//            tmpi = *OVptr++;
//            I2 = tmpi & nsignmask;
//            S2 = (tmpi & signmask) ? -1 : 1;
//            for (J=0; J<cnbs; J++) {
//               Cprime[I][J] = C[I2][J] * S2;
//               }
//            zero_arr(Sprime[I], nbs);
//            }
//
//         for (i=0; i<norbs; i++) {
//            for (j=0; j<norbs; j++) {
//               if ((orbsym[i] ^ orbsym[j] ^ Jbsym ^ Ibsym) != 0) continue;
//               ij = INDEX(i,j);
//               if (ij > kl) continue;
//               V = Tptr[ij];
//               if (ij==kl) V = V/2.0;
//               fullij = i * norbs + j;
//               fullji = j * norbs + i;
//               jlen = OV[Jblist][fullji][0];
//               OVptr = OV[Jblist][fullji] + 1;
//               OVptr2 = OV[Iblist][fullij] + 1;
//
//
//               for (J=0; J<jlen; J++) {
//                  tmpi = OVptr[J];
//                  tmpj = OVptr2[J];
//                  J1 = tmpi & nsignmask;
//                  J2 = tmpj & nsignmask;
//                  VS = (tmpj & signmask) ? -V : V;
//                  for (I=0; I<ilen%8; I++) {
//                     Sprime[I][J2] += Cprime[I][J1] * VS;
//                     }
//                  }
//
//               for (; I<ilen; I+=8) {
//
//                  SprimeI0=Sprime[I];    SprimeI1=Sprime[I+1];
//                  SprimeI2=Sprime[I+2];  SprimeI3=Sprime[I+3];
//                  SprimeI4=Sprime[I+4];  SprimeI5=Sprime[I+5];
//                  SprimeI6=Sprime[I+6];  SprimeI7=Sprime[I+7];
//                  CprimeI0=Cprime[I];    CprimeI1=Cprime[I+1];
//                  CprimeI2=Cprime[I+2];  CprimeI3=Cprime[I+3];
//                  CprimeI4=Cprime[I+4];  CprimeI5=Cprime[I+5];
//                  CprimeI6=Cprime[I+6];  CprimeI7=Cprime[I+7];
//
//                  for (J=0; J<jlen; J++) {
//                     tmpi = OVptr[J];
//                     tmpj = OVptr2[J];
//                     J1 = tmpi & nsignmask;
//                     J2 = tmpj & nsignmask;
//                     VS = (tmpj & signmask) ? -V : V;
//
//                     SprimeI0[J2] += CprimeI0[J1] * VS;
//                     SprimeI1[J2] += CprimeI1[J1] * VS;
//                     SprimeI2[J2] += CprimeI2[J1] * VS;
//                     SprimeI3[J2] += CprimeI3[J1] * VS;
//                     SprimeI4[J2] += CprimeI4[J1] * VS;
//                     SprimeI5[J2] += CprimeI5[J1] * VS;
//                     SprimeI6[J2] += CprimeI6[J1] * VS;
//                     SprimeI7[J2] += CprimeI7[J1] * VS;
//                     }
//                  }
//
//               } /* end loop over j */
//            } /* end loop over i */
//
//         for (I=0,OVptr=OV[Ialist][fullkl]+1; I<ilen; I++) {
//            tmpi = *OVptr++;
//            I1 = tmpi & nsignmask;
//            for (J=0; J<nbs; J++) {
//               S[I1][J] += Sprime[I][J];
//               }
//            }
//
//         } /* end loop over l */
//      } /* end loop over k */
//}

/*
** S3_BLOCK_BZ2()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** Try routine from Bendazzoli and Evangelisti
**
** This one reverses the direction of vectorization
**
** Pass in Ialist, Iblist, Jalist, Jblist.  Need to convert these
** to irreps for now; later we can make even more general.
**
*/
void s3_block_bz(int Ialist, int Iblist, int Jalist, int Jblist,
      int nas, int nbs,
      int cnas, double *tei, double **C, double **S,
      double **Cprime, double **Sprime, struct calcinfo *CInfo, int ***OV)
{

   int Iasym, Jasym, Ibsym, Jbsym;
   int norbs, *orbsym;
   int i, j, k, l, ij, fullij, fullji, kl, fullkl, fulllk, tmpi, tmpj;
   int I, J, I1, I2, J1, J2, S2;   /* also try S2 as double */
   int ilen, jlen;
   double V, VS;
   double *Tptr;
   int *OVptr, *OVptr2;
   int signmask, nsignmask;
   double *SprimeI, *CprimeI;
   double *SprimeI0,*SprimeI1,*SprimeI2,*SprimeI3;
   double *SprimeI4,*SprimeI5,*SprimeI6,*SprimeI7;
   double *SI0, *SI1, *SI2, *SI3, *SI4, *SI5, *SI6, *SI7;

   orbsym = CInfo->orbsym + CInfo->num_drc_orbs;
   Iasym = Ialist;
   Jasym = Jalist;
   Ibsym = Iblist;
   Jbsym = Jblist;

   norbs = CInfo->num_ci_orbs;

   signmask = 1 << (sizeof(int)*8-1);
   nsignmask = ~signmask;

   /* loop over i,j */
   for (i=0; i<norbs; i++) {
      for (j=0; j<norbs; j++) {
         if ((orbsym[i] ^ orbsym[j] ^ Jbsym ^ Ibsym) != 0) continue;

         ij = INDEX(i,j);
         Tptr = tei + ioff[ij];
         fullij = i * norbs + j;
         fullji = j * norbs + i;

         jlen = OV[Jblist][fullji][0];

         if (jlen == 0) continue;

         for (J=0,OVptr=OV[Jblist][fullji]+1; J<jlen; J++) {
            tmpi = *OVptr++;
            J2 = tmpi & nsignmask;
            S2 = (tmpi & signmask) ? -1 : 1;
            for (I=0; I<cnas%8; I++) {
               Cprime[I][J] = C[I][J2] * S2;
               }
            }

         OVptr = OV[Jblist][fullji]+1;
         for (; I<cnas; I+=8) {
            SprimeI0 = Cprime[I];    SI0 = C[I];
            SprimeI1 = Cprime[I+1];  SI1 = C[I+1];
            SprimeI2 = Cprime[I+2];  SI2 = C[I+2];
            SprimeI3 = Cprime[I+3];  SI3 = C[I+3];
            SprimeI4 = Cprime[I+4];  SI4 = C[I+4];
            SprimeI5 = Cprime[I+5];  SI5 = C[I+5];
            SprimeI6 = Cprime[I+6];  SI6 = C[I+6];
            SprimeI7 = Cprime[I+7];  SI7 = C[I+7];

            for (J=0; J<jlen; J++) {
               tmpi = OVptr[J];
               J2 = tmpi & nsignmask;
               S2 = (tmpi & signmask) ? -1 : 1;

               SprimeI0[J] = SI0[J2] * S2;
               SprimeI1[J] = SI1[J2] * S2;
               SprimeI2[J] = SI2[J2] * S2;
               SprimeI3[J] = SI3[J2] * S2;
               SprimeI4[J] = SI4[J2] * S2;
               SprimeI5[J] = SI5[J2] * S2;
               SprimeI6[J] = SI6[J2] * S2;
               SprimeI7[J] = SI7[J2] * S2;
               }
            }

         zero_mat(Sprime, nas, nbs);

         for (k=0; k<norbs; k++) {
            for (l=0; l<norbs; l++) {
               if ((orbsym[k] ^ orbsym[l] ^ Jasym ^ Iasym) != 0) continue;
               kl = INDEX(k,l);
               if (kl > ij) continue;
               V = Tptr[kl];
               if (ij==kl) V = V/2.0;
               fullkl = k * norbs + l;
               fulllk = l * norbs + k;
               ilen = OV[Jalist][fulllk][0];
               OVptr = OV[Jalist][fulllk] + 1;
               OVptr2 = OV[Ialist][fullkl] + 1;

               for (I=0; I<ilen; I++) {
                  tmpi = *OVptr++;
                  tmpj = *OVptr2++;
                  I1 = tmpi & nsignmask;
                  I2 = tmpj & nsignmask;
                  VS = (tmpj & signmask) ? -V : V;

               #ifdef USE_BLAS
                  C_DAXPY(jlen, VS, Cprime[I1], 1, Sprime[I2], 1);
               #else
                  SprimeI = Sprime[I2];
                  CprimeI = Cprime[I1];

                  for (J=0; J<jlen%8; J++) {
                     SprimeI[J] += CprimeI[J] * VS;
                     }

                  for (; J<jlen; J+=8) {
                     SprimeI[J] += CprimeI[J] * VS;
                     SprimeI[J+1] += CprimeI[J+1] * VS;
                     SprimeI[J+2] += CprimeI[J+2] * VS;
                     SprimeI[J+3] += CprimeI[J+3] * VS;
                     SprimeI[J+4] += CprimeI[J+4] * VS;
                     SprimeI[J+5] += CprimeI[J+5] * VS;
                     SprimeI[J+6] += CprimeI[J+6] * VS;
                     SprimeI[J+7] += CprimeI[J+7] * VS;
                     }
               #endif

                  } /* end loop over I */
               } /* end loop over l */
            } /* end loop over k */

         for (J=0,OVptr=OV[Iblist][fullij]+1; J<jlen; J++) {
            tmpi = *OVptr++;
            J1 = tmpi & nsignmask;
            for (I=0; I<nas%8; I++) {
               S[I][J1] += Sprime[I][J];
               }
            }

         OVptr = OV[Iblist][fullij]+1;

         for (; I<nas; I+=8) {
            SI0 = S[I];  SprimeI0 = Sprime[I];
            SI1 = S[I+1];  SprimeI1 = Sprime[I+1];
            SI2 = S[I+2];  SprimeI2 = Sprime[I+2];
            SI3 = S[I+3];  SprimeI3 = Sprime[I+3];
            SI4 = S[I+4];  SprimeI4 = Sprime[I+4];
            SI5 = S[I+5];  SprimeI5 = Sprime[I+5];
            SI6 = S[I+6];  SprimeI6 = Sprime[I+6];
            SI7 = S[I+7];  SprimeI7 = Sprime[I+7];

            for (J=0; J<jlen; J++) {
               tmpi = OVptr[J];
               J1 = tmpi & nsignmask;
               SI0[J1] += SprimeI0[J];
               SI1[J1] += SprimeI1[J];
               SI2[J1] += SprimeI2[J];
               SI3[J1] += SprimeI3[J];
               SI4[J1] += SprimeI4[J];
               SI5[J1] += SprimeI5[J];
               SI6[J1] += SprimeI6[J];
               SI7[J1] += SprimeI7[J];
               }
            }

         } /* end loop over j */
      } /* end loop over i */
}

}} // namespace psi::detci
