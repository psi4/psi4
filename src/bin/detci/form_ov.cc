/*! \file
**  \ingroup DETCI
**  \brief Form OV arrays of Bendazzoli and Evangelisti, JCP 98, 3141 (1993)
**
** David Sherrill
** University of Georgia
** 8 April 1996
**
*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

/*
** FORM_OV()
** This will only work for Full CI's right now (where Parameters.fci=true)
*/
void form_ov(struct stringwr **alplist)
{

   int i, j, nirreps, norbs;
   int irrep, strnum, strsym, cnt=0;
   int fullij, idx, ovcnt;
   struct stringwr *strlist;
   int signmask,nsignmask;


   /* bitwise sign stuff */

   signmask = 1 << (sizeof(int)*8-1);
   nsignmask = ~signmask;


   /* allocate memory for OV[list][fullij][string] */

   norbs = CalcInfo.num_ci_orbs;
   nirreps = AlphaG->nirreps;
   OV = (int ***) malloc (sizeof(int **) * nirreps);
   for (i=0; i<nirreps; i++) {
      OV[i] = (int **) malloc (sizeof(int *) * norbs * norbs);
      for (j=0; j<norbs*norbs; j++) {
         OV[i][j] = (int *) malloc (sizeof(int) * AlphaG->max_str_per_irrep+1);
         OV[i][j][0] = 0;
         }
      }


   /* now fill up OV by walking through the stringwr lists */

   for (irrep=0; irrep < nirreps; irrep++) {
      strnum = AlphaG->sg[irrep][0].num_strings;
      cnt=0;
      strlist = alplist[irrep];
      while (cnt != strnum) { 
         for (strsym=0; strsym < nirreps; strsym++) {
            for (i=0; i<strlist->cnt[strsym]; i++) {
               fullij = strlist->oij[strsym][i];
               /* idx = cnt + 1; */
               idx = cnt;
               if (strlist->sgn[strsym][i] != 1) idx = idx | signmask;
               ovcnt = OV[irrep][fullij][0];
               ovcnt++;
               OV[irrep][fullij][ovcnt] = idx;
               OV[irrep][fullij][0] = ovcnt;
               }  
            }
         strlist++;
         cnt++;
         }
      }


   /* print out the OV data */

   if (Parameters.print_lvl > 3) {
      for (irrep=0; irrep < nirreps; irrep++) {
         for (fullij=0; fullij<norbs*norbs; fullij++) {
            fprintf(outfile, "OV[irrep=%d][oij=%d]:  ", irrep, fullij);
            for (i=0; i<OV[irrep][fullij][0]; i++) {
               idx = OV[irrep][fullij][i+1];
               fprintf(outfile, "%c", (idx & signmask) ? '-' : '+');
               idx = idx & nsignmask;
               fprintf(outfile, "%2d ", idx);
               }
            fprintf(outfile, "\n");
            }
         }
      }


}

}} // namespace psi::detci

