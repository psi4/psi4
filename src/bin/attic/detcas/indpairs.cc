/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** INDPAIRS.C
** 
** Contains code pertaining to the "independent pairs" of orbitals for which
** the energy is not invariant.  Only these pairs need to be considered in
** the MCSCF orbital rotation procedure.
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

#include <cstdio>
#include <cstdlib>
#include "indpairs.h"

namespace psi { namespace detcas {

IndepPairs::IndepPairs() // Default constructor
{
  npairs = 0;
  nirreps = 0;
  p = NULL;
  q = NULL;
  map_pair_ir = NULL;
  map_pair_rel = NULL;
  ir_npairs = NULL;
  ir_p = NULL;
  ir_q = NULL;
  ir_p_rel = NULL;
  ir_q_rel = NULL;
  ir_map_pair = NULL;
}


// Regular constructor
IndepPairs::IndepPairs(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
  int *fzc, int **fzc_orbs, int *cor, int **cor_orbs,
  int *vir, int **vir_orbs, int *fzv, int **fzv_orbs,
  int *ci2relpitz, int ignore_ras_ras, int ignore_fz) 
{

  set(nirr, num_ras, ras_opi, ras_orbs, fzc, fzc_orbs, cor, cor_orbs,
      vir, vir_orbs, fzv, fzv_orbs, ci2relpitz, ignore_ras_ras, ignore_fz);

}

/*!
** set()
**
** This function sets the independent pairs info
**
** \param nirr           = number of irreps in point group
** \param num_ras        = number of RAS spaces
** \param ras_opi        = number of orbitals per irrep per RAS space
** \param ras_orbs       = ras_orbs[ras][irr][cnt] gives an orbital number
** \param fzc            = number of frozen core orbs per irrep
** \param fzc_orbs       = fzc_orbs[irr][cnt] gives an orbital number
** \param cor            = number of restricted core orbs per irrep
** \param cor_orbs       = cor_orbs[irr][cnt] gives an orbital number
** \param vir            = number of restricted virtual orbitals per irrep
** \param vir_orbs       = vir_orbs[irr][cnt] gives an orbital number
** \param fzv            = number of frozen virtual orbitals per irrep 
** \param fzv_orbs       = fzv_orbs[irr][cnt] gives an orbital number
** \param ci2relpitz     = maps CI orbitals to relative Pitzer indices
** \param ignore_ras_ras = ignore RAS/RAS independent pairs
** \param ignore_fz      = ignore FZC and FZV from all independent pairs
** 
*/

void IndepPairs::set(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
                     int *fzc, int **fzc_orbs,
                     int *cor, int **cor_orbs,
                     int *vir, int **vir_orbs,
                     int *fzv, int **fzv_orbs,
                     int *ci2relpitz, int ignore_ras_ras, int ignore_fz) 
{
  int count, ir_count, *ir_cnt;
  int rasi, rasj, irrep;
  int i, j;

  nirreps = nirr;
  ir_npairs = new int[nirreps];

  // count everything up!

  for (irrep=0; irrep<nirreps; irrep++) ir_npairs[irrep] = 0;

  if (!ignore_fz) {
    // FZC / COR
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += fzc[irrep] * cor[irrep];

    // FZC / RAS
    for (rasj=0; rasj<num_ras; rasj++)
      for (irrep=0; irrep<nirreps; irrep++) 
        ir_npairs[irrep] += ras_opi[rasj][irrep] * fzc[irrep];

    // FZC / VIR
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += fzc[irrep] * vir[irrep];

    // FZC / FZV
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += fzc[irrep] * fzv[irrep];

  }

  // COR / RAS
  for (rasj=0; rasj<num_ras; rasj++)
    for (irrep=0; irrep<nirreps; irrep++) 
      ir_npairs[irrep] += ras_opi[rasj][irrep] * cor[irrep];

  // COR / VIR
  for (irrep=0; irrep<nirreps; irrep++)
    ir_npairs[irrep] += cor[irrep] * vir[irrep];

  // RAS / RAS
  if (!ignore_ras_ras) {
    for (rasj=0; rasj<num_ras; rasj++)
      for (rasi=rasj+1; rasi<num_ras; rasi++)
        for (irrep=0; irrep<nirreps; irrep++)
          ir_npairs[irrep] += ras_opi[rasi][irrep] * ras_opi[rasj][irrep];
  }

  // VIR / RAS
  for (rasj=0; rasj<num_ras; rasj++)
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += ras_opi[rasj][irrep] * vir[irrep];

  if (!ignore_fz) {
    // FZV / COR
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += fzv[irrep] * cor[irrep];
    
    // FZV / RAS
    for (rasj=0; rasj<num_ras; rasj++)
      for (irrep=0; irrep<nirreps; irrep++)
        ir_npairs[irrep] += ras_opi[rasj][irrep] * fzv[irrep];

    // FZV / VIR
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += fzv[irrep] * vir[irrep];
  }

  npairs = 0;
  for (irrep=0; irrep<nirreps; irrep++) 
    npairs += ir_npairs[irrep];

  if (npairs==0)  {
    printf("(IndepPairs): Constructor called but no pairs!!\n");
    fprintf(stderr, "(IndepPairs): Constructor called but no pairs!!\n");
    return;
  }

  // Start making the pairs --- first, allocate the arrays
  p = new int[npairs];
  q = new int[npairs];
  map_pair_ir  = new int[npairs];
  map_pair_rel = new int[npairs]; 

  ir_p = new int*[nirreps];
  ir_q = new int*[nirreps];
  ir_p_rel = new int*[nirreps];
  ir_q_rel = new int*[nirreps];
  ir_map_pair = new int*[nirreps];

  for (irrep=0; irrep<nirreps; irrep++) {
    i = ir_npairs[irrep];
    if (!i) continue;
    ir_p[irrep] = new int[i];
    ir_q[irrep] = new int[i];
    ir_p_rel[irrep] = new int[i];
    ir_q_rel[irrep] = new int[i];
    ir_map_pair[irrep] = new int[i];
  }

  ir_cnt = new int[nirreps];
  for (irrep=0; irrep<nirreps; irrep++) ir_cnt[irrep] = 0; 
  count = 0;

  // Now put everything in the proper arrays

  if (!ignore_fz) {

    // FZC / COR
    set_part(count, ir_cnt, cor, fzc, cor_orbs, fzc_orbs, ci2relpitz);

    // FZC / RAS
    for (rasj=0; rasj<num_ras; rasj++) {
      set_part(count, ir_cnt, ras_opi[rasj], fzc, ras_orbs[rasj], fzc_orbs,
               ci2relpitz); 
    }

    // FZC / VIR
    set_part(count, ir_cnt, vir, fzc, vir_orbs, fzc_orbs, ci2relpitz);

    // FZC / FZV
    set_part(count, ir_cnt, fzv, fzc, fzv_orbs, fzc_orbs, ci2relpitz);

  }

  // COR / RAS
  for (rasj=0; rasj<num_ras; rasj++) {
    set_part(count, ir_cnt, ras_opi[rasj], cor, ras_orbs[rasj], cor_orbs,
             ci2relpitz); 
  }

  // COR / VIR
  set_part(count, ir_cnt, vir, cor, vir_orbs, cor_orbs, ci2relpitz); 

  // COR / FZV
  if (!ignore_fz)
    set_part(count, ir_cnt, fzv, cor, fzv_orbs, cor_orbs, ci2relpitz);
 
  // RAS / RAS
  if (!ignore_ras_ras) {
    for (rasj=0; rasj<num_ras; rasj++) {
      for (rasi=rasj+1; rasi<num_ras; rasi++) {
        set_part(count, ir_cnt, ras_opi[rasi], ras_opi[rasj],
          ras_orbs[rasi], ras_orbs[rasj], ci2relpitz);
      }
    }
  }

  // RAS / VIR
  for (rasj=0; rasj<num_ras; rasj++) {
    set_part(count, ir_cnt, vir, ras_opi[rasj], vir_orbs, ras_orbs[rasj], 
             ci2relpitz); 
  }

  // RAS / FZV
  if (!ignore_fz) {
    for (rasj=0; rasj<num_ras; rasj++) {
      set_part(count, ir_cnt, fzv, ras_opi[rasj], fzv_orbs, ras_orbs[rasj], 
               ci2relpitz); 
    }
  }

  // VIR / FZV
  if (!ignore_fz)
    set_part(count, ir_cnt, fzv, vir, fzv_orbs, vir_orbs, ci2relpitz);

  // check things
  if (count != npairs) {
    printf("(IndepPairs::set): mismatch in counted pairs!\n");
    fprintf(stderr, "(IndepPairs::set): mismatch in counted pairs!\n");
  }


}


/*!
** set_part()
**
** This function does the dirty work in setting up the independent pair
** arrays.
**
** \param count      = number of independent pairs so far
** \param ir_cnt     = count per irrep
** \param num_orbs_i = number of orbitals per irrep, pair i
** \param num_orbs_j = number of orbitals per irrep, pair j
** \param orbs_i     = orbital number for [irrep][cnt], pair i
** \param orbs_j     = orbital number for [irrep][cnt], pair j
** \param ci2relpitz = map CI orbital to relative Pitzer numbering
*/
void IndepPairs::set_part(int &count, int *ir_cnt, 
  int *num_orbs_i, int *num_orbs_j, int **orbs_i, int **orbs_j, 
  int *ci2relpitz)
{
  int i, j, irrep;
  int ir_count;

  for (irrep=0; irrep<nirreps; irrep++) {
    ir_count = ir_cnt[irrep];
    for (i=0; i<num_orbs_i[irrep]; i++) {
      for (j=0; j<num_orbs_j[irrep]; j++) {
        p[count] = orbs_i[irrep][i];
        q[count] = orbs_j[irrep][j];
        map_pair_ir[count] = irrep;
        map_pair_rel[count] = ir_count;
        ir_p[irrep][ir_count] = p[count];
        ir_q[irrep][ir_count] = q[count];
        ir_p_rel[irrep][ir_count] = ci2relpitz[p[count]];
        ir_q_rel[irrep][ir_count] = ci2relpitz[q[count]];
        ir_map_pair[irrep][ir_count] = count;
        count++;
        ir_count++;
      }
    }
  ir_cnt[irrep] = ir_count;
  }

}


IndepPairs::~IndepPairs() // Destructor
{
  int h;

  if (npairs) {
    delete [] p;
    delete [] q;
    delete [] map_pair_ir;
    delete [] map_pair_rel;
  }

  if (ir_npairs != NULL) {
    for (h=0; h<nirreps; h++) {
      if (ir_npairs[h]) {
        delete [] ir_p[h];
        delete [] ir_q[h];
        delete [] ir_p_rel[h];
        delete [] ir_q_rel[h];
        delete [] ir_map_pair[h];
      }  
    }
  }

  if (npairs) {
    delete [] ir_npairs;
    delete [] ir_p;
    delete [] ir_q;
    delete [] ir_p_rel;
    delete [] ir_q_rel;
    delete [] ir_map_pair;
  }

}


void IndepPairs::print(FILE *outfile)
{
  int h;

  fprintf(outfile, "\nList of all independent pairs:\n");
  print_selected(npairs, p, q, outfile);
  fprintf(outfile, "\nLists of independent pairs by irrep:\n");
  
  for (h=0; h<nirreps; h++) {
    if (!ir_npairs[h]) continue;
    fprintf(outfile, "\n\t Irrep %d:\n", h);
    print_selected(ir_npairs[h], ir_p[h], ir_q[h], 
                   ir_p_rel[h], ir_q_rel[h], outfile);
  }

}


void IndepPairs::print_selected(int num, int *parr, int *qarr, FILE *outfile)
{
  int ii;

  fprintf(outfile, "\n  %4d Independent Pairs\n", num);
  fprintf(outfile, "\t p\t q\n");
  fprintf(outfile,   "    -------------------\n");

  for (ii=0; ii<num; ii++) {
    fprintf(outfile,"\t %2d\t",parr[ii]);
    fprintf(outfile,"%2d\n",qarr[ii]);
  }
  fflush(outfile);

}

void IndepPairs::print_selected(int num, int *parr, int *qarr, 
                                int *prel, int *qrel, FILE *outfile)
{
  int ii;

  fprintf(outfile, "\n\t %4d Independent Pairs\n", num);
  fprintf(outfile, "\t p\t q\t P\t Q\n");
  fprintf(outfile,   "    ----------------------------------\n");

  for (ii=0; ii<num; ii++) {
    fprintf(outfile,"\t %2d\t%2d\t%2d\t%2d\n",
            parr[ii],qarr[ii],prel[ii],qrel[ii]);
  }
  fflush(outfile);

}


void IndepPairs::print_vec(double *arr, const char *label, FILE *outfile) 
{
  int pair;

  fprintf(outfile, "%s\n", label);
  for (pair=0; pair<npairs; pair++) {
    fprintf(outfile, "Pair (%2d,%2d) = %12.7lf\n",
            p[pair], q[pair], arr[pair]);
  }
  fprintf(outfile, "\n");
  fflush(outfile);

}

/*
** get_irrep_vec()
**
** This function gets the piece of a vector belonging to a given irrep.
**  Assumes that the ordering of the input array is consistent with the
**  ordering of independent pairs in the total list of pairs, and has
**  the same length as the total number of independent pairs.
*/
double * IndepPairs::get_irrep_vec(int irrep, double *arr)
{

  int pair, target_len;
  double *newarr;

  target_len = ir_npairs[irrep];
  if (target_len==0)  return (NULL);

  newarr = new double[target_len];
  for (pair=0; pair<target_len; pair++) {
    newarr[pair] = arr[ir_map_pair[irrep][pair]];
  }

  return(newarr);

}


/*
** put_irrep_vec()
**
** This function scatters a vector for an irrep to the overall vector
**  using the mapping arrays in the independent pairs class.  This is
**  basically the reverse of get_irrep_vec().
*/
void IndepPairs::put_irrep_vec(int irrep, double *ir_vec, double *tot_vec)
{
  int pair;

  for (pair=0; pair<ir_npairs[irrep]; pair++) {
    tot_vec[ir_map_pair[irrep][pair]] = ir_vec[pair];
  }

} 


}} // end namespace psi::detcas

