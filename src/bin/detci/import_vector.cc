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
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterdset.h>
#include <physconst.h>
#include "structs.h"
#include "ci_tol.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
  int *listnum);
void stringset_translate_addr(StringSet *sset, int new_nel, int new_nfzc,
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
void parse_import_vector(SlaterDetSet *sdset, int *ialplist, int *ialpidx,
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
  if (CalcInfo.num_fzc_orbs > alphastrings->nfzc ||
      CalcInfo.num_fzc_orbs > betastrings->nfzc) {
    fprintf(outfile, "(parse_import_vector): Can't freeze more orbitals now" \
      " than in the imported guess!\n");
    abort();
  }

  if (alphastrings->nfzc != betastrings->nfzc) {
    fprintf(outfile, "(parse_import_vector): alpha nfzc != beta nfzc!\n");
    abort();
  }

  new_alphastr_list = init_int_array(alphastrings->size);
  new_alphastr_idx  = init_int_array(alphastrings->size);
  new_betastr_list  = init_int_array(betastrings->size);
  new_betastr_idx   = init_int_array(betastrings->size);

  stringset_translate_addr(alphastrings, CalcInfo.num_alp_expl, 
    CalcInfo.num_fzc_orbs, CalcInfo.reorder, AlphaG, new_alphastr_list,
    new_alphastr_idx);
  stringset_translate_addr(betastrings, CalcInfo.num_bet_expl, 
    CalcInfo.num_fzc_orbs, CalcInfo.reorder, BetaG, new_betastr_list,
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
    j = CIblks.decode[ialplist[i]][ibetlist[i]];
    if (j == -1) {
      fprintf(outfile, "Import vector: can't find CI block!\n");
      fprintf(outfile, "Determinant number %d\n", i);
      fprintf(outfile, "\nialplist=%d, ialpidx=%d, ibetlist=%d, ibetidx=%d\n", 
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
void stringset_translate_addr(StringSet *sset, int new_nel, int new_nfzc,
  int *pitz2corr, struct olsen_graph *Graph, int *new_list, int *new_idx)
{

  int i, j, l, s;
  int old_nel;
  int *former_fzc_occ, num_former_fzc, *tmpocc;
  short int *old_occ;

  old_nel = sset->nelec - sset->nfzc;

  former_fzc_occ = init_int_array(sset->nfzc);

  for (i=0,num_former_fzc=0; i<sset->nfzc; i++) {
    j = (int) sset->fzc_occ[i];
    j = pitz2corr[j] - new_nfzc;
    if (j >= 0) former_fzc_occ[num_former_fzc++] = j;
  }

  if (num_former_fzc + old_nel > new_nel) {
    fprintf(outfile, "(stringset_translate_addr): num_former_fzc + old_nel" \
      " > new_nel!\n");
    fflush(outfile);
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

    for (i=0,l=0; i<num_former_fzc; i++)
      tmpocc[l++] = former_fzc_occ[i]; /* these are already translated */

    for (i=0; i<old_nel; i++) {
      j = (int) old_occ[i];
      j = pitz2corr[j] - new_nfzc;
      if (j >= 0) tmpocc[l++] = j;
    }

    if (l != new_nel) {
      fprintf(outfile, "(stringset_translate_addr): Imported string has wrong" \
        " number of electrons, %d vs. %d\n", l, new_nel);
      fflush(outfile);
      abort();
    }

    new_idx[s] = og_lex_addr(Graph,tmpocc,new_nel,&(new_list[s]));
  }

  free(former_fzc_occ);
  free(tmpocc);

}

}} // namespace psi::detci

