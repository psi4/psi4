/*! \file shell_pairs.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<libciomr/libciomr.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"

namespace psi { namespace cints {

/*-------------------------------
  Explicit function delcarations
 -------------------------------*/
static double s_ovlp(double, double, double, double, double);
static long int memory_to_store_shell_pairs();

/*!----------------------------------------------------
  shell_pair includes pointers to dynamically
  allocated data. To minimize overhead and to be able
  to predict how much memory is needed I create a
  locally visible stack that's used to distribute
  memory chunks
  ----------------------------------------------------*/
static double *stack;

void init_shell_pairs()
{
  struct shell_pair **pairs, *sp;
  struct coordinates P, PA, PB, AB;
  struct coordinates *atom_i, *atom_j;
  int i, j;
  int si, sj, np_i, np_j;
  long int memd;
  double a1, a2, ab2, gam, inorm, jnorm;
  double *curr_stack_ptr;

  /*---
    estimate storage requirement for the dynamically
    allocated parts of shell_pair structure 
   ---*/
  memd = memory_to_store_shell_pairs();
  if (memd > UserOptions.memory)
    throw std::domain_error("Not enough memory to store shell pair data. Need more memory.");
  UserOptions.memory -= memd;
  stack = (double *) malloc(sizeof(double)*memd);
  curr_stack_ptr = stack;

  /*--- allocate room for the shell pairs ---*/
  pairs = (struct shell_pair **)malloc(sizeof(struct shell_pair *)*BasisSet.num_shells);
  for(i=0; i<BasisSet.num_shells; i++)
    pairs[i] = (struct shell_pair *)malloc(sizeof(struct shell_pair)*BasisSet.num_shells);
  /*  UserOptions.memory -= BasisSet.num_shells*BasisSet.num_shells*sizeof(struct shell_pair)/sizeof(double);*/

  /*--- loop over all shell pairs si, sj and create primitive pairs pairs ---*/
  for (si=0; si<BasisSet.num_shells; si++){
    atom_i = &(Molecule.centers[BasisSet.shells[si].center-1]);
    for (sj=0; sj<BasisSet.num_shells; sj++){
      atom_j = &(Molecule.centers[BasisSet.shells[sj].center-1]);
      AB.x = atom_i->x - atom_j->x;
      AB.y = atom_i->y - atom_j->y;
      AB.z = atom_i->z - atom_j->z;
      ab2 = AB.x*AB.x+AB.y*AB.y+AB.z*AB.z;

      sp = &(pairs[si][sj]);
      sp->i = si;
      sp->j = sj;
      sp->AB[0] = AB.x;
      sp->AB[1] = AB.y;
      sp->AB[2] = AB.z;
      np_i = BasisSet.shells[si].n_prims;
      np_j = BasisSet.shells[sj].n_prims;

      sp->a1 = curr_stack_ptr; curr_stack_ptr += np_i;
      sp->a2 = curr_stack_ptr; curr_stack_ptr += np_j;

      sp->gamma = (double **) malloc(np_i*sizeof(double *));
      for(i=0;i<np_i;i++) {
	sp->gamma[i] = curr_stack_ptr; curr_stack_ptr += np_j;
      }
      
      sp->inorm = curr_stack_ptr; curr_stack_ptr += np_i;
      sp->jnorm = curr_stack_ptr; curr_stack_ptr += np_j;
      
      sp->Sovlp = (double **) malloc(np_i*sizeof(double *));
      for(i=0;i<np_i;i++) {
	sp->Sovlp[i] = curr_stack_ptr; curr_stack_ptr += np_j;
      }

      sp->P  = (double ***) malloc(np_i*sizeof(double **));
      sp->PA = (double ***) malloc(np_i*sizeof(double **));
      sp->PB = (double ***) malloc(np_i*sizeof(double **));
      for(i=0;i<np_i;i++) {
	sp->P[i]  = (double **) malloc(np_j*sizeof(double *));
	sp->PA[i]  = (double **) malloc(np_j*sizeof(double *));
	sp->PB[i]  = (double **) malloc(np_j*sizeof(double *));
	for(j=0;j<np_j;j++) {
	  sp->P[i][j]  = curr_stack_ptr; curr_stack_ptr += 3;
	  sp->PA[i][j]  = curr_stack_ptr; curr_stack_ptr += 3;
	  sp->PB[i][j]  = curr_stack_ptr; curr_stack_ptr += 3;
	}
      }

      /*--- loop over primitives here ---*/
      for (i = 0; i < np_i; i++){
	a1 = BasisSet.cgtos[BasisSet.shells[si].fprim+i-1].exp;
	inorm = BasisSet.cgtos[BasisSet.shells[si].fprim+i-1].ccoeff[BasisSet.shells[si].am-1];
	sp->a1[i] = a1;
	sp->inorm[i] = inorm;
        for (j = 0; j < np_j; j++){
          a2 = BasisSet.cgtos[BasisSet.shells[sj].fprim+j-1].exp;
          gam = a1 + a2;
	  jnorm = BasisSet.cgtos[BasisSet.shells[sj].fprim+j-1].ccoeff[BasisSet.shells[sj].am-1];

	  P.x = (atom_i->x*a1 + atom_j->x*a2)/gam;
	  P.y = (atom_i->y*a1 + atom_j->y*a2)/gam;
	  P.z = (atom_i->z*a1 + atom_j->z*a2)/gam;
          PA.x = P.x - atom_i->x;
          PA.y = P.y - atom_i->y;
          PA.z = P.z - atom_i->z;
          PB.x = P.x - atom_j->x;
          PB.y = P.y - atom_j->y;
          PB.z = P.z - atom_j->z;

          /*--- copy init data into pairs array in prep for eri evaluation ---*/
          sp->a2[j] = a2;
          sp->gamma[i][j] = a1+a2;
          sp->jnorm[j] = jnorm;
          sp->P[i][j][0] = P.x;
	  sp->P[i][j][1] = P.y;
	  sp->P[i][j][2] = P.z;
	  sp->PA[i][j][0] = PA.x;
	  sp->PA[i][j][1] = PA.y;
	  sp->PA[i][j][2] = PA.z;
	  sp->PB[i][j][0] = PB.x;
	  sp->PB[i][j][1] = PB.y;
	  sp->PB[i][j][2] = PB.z;
          sp->Sovlp[i][j] = s_ovlp(a1,a2,inorm,jnorm,ab2);

	  sp->dmat = NULL;
	  sp->dmato = NULL;
	  sp->dmata = NULL;
	  sp->dmatb = NULL;
	}
      }
    }
  }

  BasisSet.shell_pairs = pairs;
  return;
}



void init_unique_shell_pairs(void)
{
  struct unique_shell_pair **pairs;
  int i,j,count;
  int us,usi,usj,si,sj;
  int nbfi, nbfj;
  int bf_i,bf_j,max_bf_j;
  int shell,first_bf,last_bf,num_bf,bf;
  int so,irrep,irr_i,irr_j,max_irr_j;
  int so_i, so_j;
  int *count_irr;
  int ***bf_so;

  count_irr = init_int_array(Symmetry.nirreps);

  bf_so = (int ***) malloc(Symmetry.num_unique_shells*sizeof(int **));
  for(us=0;us<Symmetry.num_unique_shells;us++) {
    shell = Symmetry.us2s[us];
    /*--- fill in bf_so ---*/
    first_bf = BasisSet.shells[shell].fbf-1;
    num_bf = BasisSet.puream ? 2*BasisSet.shells[shell].am-1 : ioff[BasisSet.shells[shell].am];
    last_bf = first_bf + num_bf;
    bf_so[us] = init_int_matrix(num_bf,Symmetry.nirreps);
    for(bf=first_bf;bf<last_bf;bf++) {
      so = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	bf_so[us][bf-first_bf][irrep] = -1;
	for(i=0;i<Symmetry.sopi[irrep];i++,so++)
	  if (Symmetry.usotao[so][bf] != 0.0) {
	    bf_so[us][bf-first_bf][irrep] = so;
	    so += Symmetry.sopi[irrep] - i;
	    break;
	  }
      }
    }
  }

  pairs = (struct unique_shell_pair **)malloc(sizeof(struct unique_shell_pair *)*Symmetry.num_unique_shells);
  for(usi=0; usi<Symmetry.num_unique_shells; usi++){
    si = Symmetry.us2s[usi];
    nbfi = BasisSet.puream ? 2*BasisSet.shells[si].am-1 : ioff[BasisSet.shells[si].am];
    pairs[usi] = (struct unique_shell_pair *)malloc(sizeof(struct unique_shell_pair)*Symmetry.num_unique_shells);
    for(usj=0; usj<Symmetry.num_unique_shells; usj++){
      sj = Symmetry.us2s[usj];
      pairs[usi][usj].SOpair_npi = init_int_array(Symmetry.nirreps);
      nbfj = BasisSet.puream ? 2*BasisSet.shells[sj].am-1 : ioff[BasisSet.shells[sj].am];

      for(bf_i=0;bf_i<nbfi;bf_i++) {
	max_bf_j = (usi == usj) ? bf_i+1 : nbfj;
	for(bf_j=0;bf_j<max_bf_j;bf_j++)
	  for(irr_i=0;irr_i<Symmetry.nirreps;irr_i++) {
	    if (bf_so[usi][bf_i][irr_i] == -1) continue;
	    max_irr_j = (usi == usj && bf_i == bf_j) ? irr_i+1 : Symmetry.nirreps;
	    for(irr_j=0;irr_j<max_irr_j;irr_j++)
	      if (bf_so[usj][bf_j][irr_j] != -1)
		pairs[usi][usj].SOpair_npi[Symmetry.dp_table[irr_i][irr_j]]++;
	  }
      }

      pairs[usi][usj].SOpair_so_i = (int **) malloc(Symmetry.nirreps*sizeof(int *));
      pairs[usi][usj].SOpair_so_j = (int **) malloc(Symmetry.nirreps*sizeof(int *));
      pairs[usi][usj].SOpair_bf_i = (int **) malloc(Symmetry.nirreps*sizeof(int *));
      pairs[usi][usj].SOpair_bf_j = (int **) malloc(Symmetry.nirreps*sizeof(int *));
      for(irrep=0;irrep<Symmetry.nirreps;irrep++)
	if (i = pairs[usi][usj].SOpair_npi[irrep]) {
	  pairs[usi][usj].SOpair_so_i[irrep] = (int *) malloc(i*sizeof(int));
	  pairs[usi][usj].SOpair_so_j[irrep] = (int *) malloc(i*sizeof(int));
	  pairs[usi][usj].SOpair_bf_i[irrep] = (int *) malloc(i*sizeof(int));
	  pairs[usi][usj].SOpair_bf_j[irrep] = (int *) malloc(i*sizeof(int));
	}

      memset(count_irr,0,Symmetry.nirreps*sizeof(int));
      for(bf_i=0;bf_i<nbfi;bf_i++) {
	max_bf_j = (usi == usj) ? bf_i+1 : nbfj;
	for(bf_j=0;bf_j<max_bf_j;bf_j++)
	  for(irr_i=0;irr_i<Symmetry.nirreps;irr_i++) {
	    if (bf_so[usi][bf_i][irr_i] == -1) continue;
	    max_irr_j = (usi == usj && bf_i == bf_j) ? irr_i+1 : Symmetry.nirreps;
	    for(irr_j=0;irr_j<max_irr_j;irr_j++)
	      if (bf_so[usj][bf_j][irr_j] != -1) {
		irrep = Symmetry.dp_table[irr_i][irr_j];
		count = (count_irr[irrep]++);
		if (usi == usj && (so_i = bf_so[usi][bf_i][irr_i]) < (so_j = bf_so[usj][bf_j][irr_j])) {
		  pairs[usi][usj].SOpair_so_i[irrep][count] = so_j;
		  pairs[usi][usj].SOpair_so_j[irrep][count] = so_i;
		  pairs[usi][usj].SOpair_bf_i[irrep][count] = bf_j;
		  pairs[usi][usj].SOpair_bf_j[irrep][count] = bf_i;
		}
		else {
		  pairs[usi][usj].SOpair_so_i[irrep][count] = bf_so[usi][bf_i][irr_i];
		  pairs[usi][usj].SOpair_so_j[irrep][count] = bf_so[usj][bf_j][irr_j];
		  pairs[usi][usj].SOpair_bf_i[irrep][count] = bf_i;
		  pairs[usi][usj].SOpair_bf_j[irrep][count] = bf_j;
		}
	      }
	  }
      }
      
    }
  }

  free(count_irr);
  for(us=0;us<Symmetry.num_unique_shells;us++) {
    shell = Symmetry.us2s[us];
    num_bf = BasisSet.puream ? 2*BasisSet.shells[shell].am-1 : ioff[BasisSet.shells[shell].am];
    free_int_matrix(bf_so[us]);
  }
  free(bf_so);

  Symmetry.us_pairs = pairs;
  return;
}


void dealloc_pairs(void)
{
  int i, j, k, l, si, sj;
  struct shell_pair *sp;
  int np_i;

  free(stack);
  for(si=0;si<BasisSet.num_shells;si++)
    for(sj=0;sj<BasisSet.num_shells;sj++) {
      np_i = BasisSet.shells[si].n_prims;
      sp = &(BasisSet.shell_pairs[si][sj]);
      if (sp->dmato != NULL)
	free_block(sp->dmat);
      if (sp->dmat != NULL)
	free_block(sp->dmato);
      if (sp->dmato != NULL)
	free_block(sp->dmata);
      if (sp->dmat != NULL)
	free_block(sp->dmatb);

    free(sp->gamma);
    free(sp->Sovlp);

    if (sp->P != NULL) {
      for(i=0;i<np_i;i++)
	    free(sp->P[i]);
      free(sp->P);
    }
    if (sp->PA != NULL) {
      for(i=0;i<np_i;i++)
	    free(sp->PA[i]);
      free(sp->PA);
    }
    if (sp->PB != NULL) {
      for(i=0;i<np_i;i++)
	    free(sp->PB[i]);
      free(sp->PB);
    }
  }

  for(si=0;si<BasisSet.num_shells;si++)
    free(BasisSet.shell_pairs[si]);
  free(BasisSet.shell_pairs);

  if (Symmetry.nirreps > 1 && UserOptions.symm_ints)
    for(i=0;i<Symmetry.num_unique_shells;i++) {
      for(j=0;j<Symmetry.num_unique_shells;j++) {
	for(k=0;k<Symmetry.nirreps;k++)
	  if (l = Symmetry.us_pairs[i][j].SOpair_npi[k]) {
	    free(Symmetry.us_pairs[i][j].SOpair_so_i[k]);
	    free(Symmetry.us_pairs[i][j].SOpair_so_j[k]);
	    free(Symmetry.us_pairs[i][j].SOpair_bf_i[k]);
	    free(Symmetry.us_pairs[i][j].SOpair_bf_j[k]);
	  }
	free(Symmetry.us_pairs[i][j].SOpair_npi);
	free(Symmetry.us_pairs[i][j].SOpair_so_i);
	free(Symmetry.us_pairs[i][j].SOpair_so_j);
	free(Symmetry.us_pairs[i][j].SOpair_bf_i);
	free(Symmetry.us_pairs[i][j].SOpair_bf_j);
      }
      free(Symmetry.us_pairs[i]);
    }
  return;
}

/*!-----------------------------------------------------------------
  This computes overlap of 2 s-functions with exponents a1 and a2,
  normalization coefficients of norm1 and norm2, and square of
  distance between the centers of ab2
 -----------------------------------------------------------------*/
double s_ovlp(double a1, double a2, double norm1, double norm2, double ab2)
{
  double t, zeta, gamma;

  gamma = a1 + a2;
  t = M_PI/gamma;
  zeta = -a1*a2/gamma;

  return t*sqrt(t)*exp(zeta*ab2)*norm1*norm2;
}

/*!----------------------------------------------------
  This evaluates how much memory (in double words) is
  needed to store basic (i.e. without dmat or lagr)
  shell pair data
  ---------------------------------------------------*/
long int memory_to_store_shell_pairs()
{
  int i, j, np_i, np_j;
  long int memd = 0;

  for(i=0;i<BasisSet.num_shells;i++) {
    np_i = BasisSet.shells[i].n_prims;
    for(j=0;j<BasisSet.num_shells;j++) {
      np_j = BasisSet.shells[j].n_prims;
      memd += (2*(np_i + np_j) + 11*np_i*np_j);
    }
  }
    
  return memd;
}
}}
