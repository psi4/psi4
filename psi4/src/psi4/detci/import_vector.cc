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
**  \ingroup DETCI
**  \brief Enter brief description of file here
**
** contains code to import a vector from a previous calculation using
** StringSet, etc., functions in libqt
**
** C. David Sherrill
** August 2003
*/

#include <cstdlib>
#include <cstdio>
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libqt/slaterdset.h"

#include "psi4/physconst.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
  int *listnum);
void stringset_translate_addr(StringSet *sset, int new_nel, int new_ndrc,
  int *pitz2corr, struct olsen_graph *Graph, int *new_list, int *new_idx);


/*
** parse_import_vector
**
** This function takes an imported CI vector using the SlaterDetVector, etc,
** structures and converts it to the new CI space by determining
** an alpha list, alpha relative index, beta list, beta relative index,
** and CI block number for each imported determinant in vec.
**
** C. David Sherrill
** August 2003
*/
void CIWavefunction::parse_import_vector(SlaterDetSet *sdset, int *ialplist, int *ialpidx,
  int *ibetlist, int *ibetidx, int *blknums)
{
  int i,j;
  StringSet *alphastrings, *betastrings;
  SlaterDet *dets;
  int alphastr, betastr;
  int *new_alphastr_list, *new_alphastr_idx;
  int *new_betastr_list, *new_betastr_idx;

  dets = sdset->dets;
  alphastrings = sdset->alphastrings;
  betastrings = sdset->betastrings;

  /* now figure out how the frozen stuff is going to map */
  if (CalcInfo_->num_drc_orbs > alphastrings->ndrc ||
      CalcInfo_->num_drc_orbs > betastrings->ndrc) {
    outfile->Printf( "(parse_import_vector): Can't drop more orbitals now" \
      " than in the imported guess!\n");
    abort();
  }

  if (alphastrings->ndrc != betastrings->ndrc) {
    outfile->Printf( "(parse_import_vector): alpha ndrc != beta ndrc!\n");
    abort();
  }

  new_alphastr_list = init_int_array(alphastrings->size);
  new_alphastr_idx  = init_int_array(alphastrings->size);
  new_betastr_list  = init_int_array(betastrings->size);
  new_betastr_idx   = init_int_array(betastrings->size);

  stringset_translate_addr(alphastrings, CalcInfo_->num_alp_expl,
    CalcInfo_->num_drc_orbs, CalcInfo_->reorder.data(), AlphaG_, new_alphastr_list,
    new_alphastr_idx);
  stringset_translate_addr(betastrings, CalcInfo_->num_bet_expl,
    CalcInfo_->num_drc_orbs, CalcInfo_->reorder.data(), BetaG_, new_betastr_list,
    new_betastr_idx);

  /* loop over all the dets in the imported vector and translate
     each of them to the new determinant number.
  */
  for (i=0; i<sdset->size; i++) {
    alphastr = dets[i].alphastring;
    ialplist[i] = new_alphastr_list[alphastr];
    ialpidx[i]  = new_alphastr_idx[alphastr];
    betastr = dets[i].betastring;
    ibetlist[i] = new_betastr_list[betastr];
    ibetidx[i]  = new_betastr_idx[betastr];

    /* figure out what block we're in */
    j = CIblks_->decode[ialplist[i]][ibetlist[i]];
    if (j == -1) {
      outfile->Printf( "Import vector: can't find CI block!\n");
      outfile->Printf( "Determinant number %d\n", i);
      outfile->Printf( "\nialplist=%d, ialpidx=%d, ibetlist=%d, ibetidx=%d\n",
        ialplist[i], ialpidx[i], ibetlist[i], ibetidx[i]);
      abort();
    }
    else blknums[i] = j;
  } /* end loop over determinants */

  free(new_alphastr_list);  free(new_alphastr_idx);
  free(new_betastr_list);  free(new_betastr_idx);
}


/*
** stringset_translate_addr
**
** This function takes a StringSet and translates the occupations
** for each string, taking into account that some formerly frozen orbitals
** may now be unfrozen, and that the StringSet occupations are stored
** in Pitzer order and we need correlated order, and produces a new
** occupation list which is used to generate a new graph list
** (new_list[s]) and relative index (new_idx[s]).  The original StringSet
** is untouched because it does not carry enough storage space to hold
** both a list number and relative index (and the number of active electrons
** in the occupations array may also have changed).
**
** C. David Sherrill
** August 2003
*/
void stringset_translate_addr(StringSet *sset, int new_nel, int new_ndrc,
  int *pitz2corr, struct olsen_graph *Graph, int *new_list, int *new_idx)
{

  int i, j, l, s;
  int old_nel;
  int *former_drc_occ, num_former_drc, *tmpocc;
  short int *old_occ;

  old_nel = sset->nelec - sset->ndrc;

  former_drc_occ = init_int_array(sset->ndrc);

  for (i=0,num_former_drc=0; i<sset->ndrc; i++) {
    j = (int) sset->drc_occ[i];
    j = pitz2corr[j] - new_ndrc;
    if (j >= 0) former_drc_occ[num_former_drc++] = j;
  }

  if (num_former_drc + old_nel > new_nel) {
    outfile->Printf( "(stringset_translate_addr): num_former_drc + old_nel" \
      " > new_nel!\n");

    abort();
  }

  tmpocc = init_int_array(new_nel);

  /* Loop over all the strings in the imported stringset and translate
     each of them to the new lexical string address.  We won't store
     the translated occs array itself simply because the size might
     have changed.
  */
  for (s=0; s<sset->size; s++) {
    old_occ = sset->strings[s].occ;

    for (i=0,l=0; i<num_former_drc; i++)
      tmpocc[l++] = former_drc_occ[i]; /* these are already translated */

    for (i=0; i<old_nel; i++) {
      j = (int) old_occ[i];
      j = pitz2corr[j] - new_ndrc;
      if (j >= 0) tmpocc[l++] = j;
    }

    if (l != new_nel) {
      outfile->Printf( "(stringset_translate_addr): Imported string has wrong" \
        " number of electrons, %d vs. %d\n", l, new_nel);

      abort();
    }

    new_idx[s] = og_lex_addr(Graph,tmpocc,new_nel,&(new_list[s]));
  }

  free(former_drc_occ);
  free(tmpocc);

}

}} // namespace psi::detci
