/*! \file basisset.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>

#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"shell_pairs.h"
#include"small_fns.h"

namespace psi { namespace cints {

/*--------------------------------
  Explicit functions declarations
 --------------------------------*/
static void get_primitives(void);
static void get_shell_info(void);

void init_basisset()
{
  BasisSet.num_shells = chkpt_rd_nshell();
  BasisSet.num_prims = chkpt_rd_nprim();
  BasisSet.num_ao = chkpt_rd_nao();
  BasisSet.am2shell = chkpt_rd_am2canon_shell_order();
  BasisSet.shells_per_am = chkpt_rd_shells_per_am();
  BasisSet.max_am = chkpt_rd_max_am()+1;
  BasisSet.puream = chkpt_rd_puream();
/* BasisSet.cgtos = */ get_primitives();
/* BasisSet.shells = */ get_shell_info();
/* BasisSet.shell_pairs = */ init_shell_pairs();

  /*-----------------------------------------------
    Namespaces are not well defined in CINTS,
    because of this here're some overlapping inits
   -----------------------------------------------*/
  if (BasisSet.puream && UserOptions.symm_ints)
    Symmetry.usotao = chkpt_rd_usotbf();
  else
    Symmetry.usotao = chkpt_rd_usotao();
  if (Symmetry.nirreps > 1 && UserOptions.symm_ints)
/* Symmetry.us_pairs = */ init_unique_shell_pairs();

  return;
}


void cleanup_basisset()
{
  int i;
  
  dealloc_pairs();
  for(i=0;i<BasisSet.num_shells;i++)
    free(BasisSet.shells[i].trans_vec);
  free(BasisSet.shells);
  free(BasisSet.cgtos);

  return;
}

void get_shell_info()
{
   int i, j, l, g, count, stab_index;
   int *shell_center;		/* atomic center of each shell */
   int *shell_type;		/* angular mom. of shell */
   int *shell_num_prims;	/* number of primitives per shell */
   int *prim_pointers;		/* first primitive in shell */
   int *shell_fbf;              /* first basisfn in shell */
   int *shell_fao;              /* first AO in shell */
   int **shell_trans_table;     /* shell transformation table */

   /*--- retrieve location of shells (which atom it's centered on) ---*/
   shell_center = chkpt_rd_snuc();

   /*--- retrieve angular momentum of each shell (1=s, 2=p, 3=d, etc  ) ---*/
   shell_type = chkpt_rd_stype();

   /*--- retrieve number of primitives per shell ---*/
   shell_num_prims = chkpt_rd_snumg();

   /*--- retrieve pointer to first primitive in shell ---*/
   prim_pointers = chkpt_rd_sprim();

   /*--- retrieve pointer to first basisfn in shell ---*/
   shell_fbf = chkpt_rd_sloc_new();

   /*--- retrieve pointer to first AO in shell ---*/
   shell_fao = chkpt_rd_sloc();
   
   /*--- retrieve shell tranformation table ---*/
   shell_trans_table = chkpt_rd_shell_transm();

   /*--- retrieve maximum am ---*/
   if (BasisSet.max_am > CINTS_MAX_AM)
     throw std::domain_error("Angular momentum limit of CINTS exceeded, reconfigure and recompile");
   
   BasisSet.shells = (struct shell_def *) malloc(sizeof(struct shell_def)*
						 BasisSet.num_shells);
   BasisSet.max_num_prims = 0;
   for (i=0; i<BasisSet.num_shells; i++){
      BasisSet.shells[i].center = shell_center[i];
      BasisSet.shells[i].am = shell_type[i];
      BasisSet.shells[i].n_prims = shell_num_prims[i];
      if (shell_num_prims[i] > BasisSet.max_num_prims)
        BasisSet.max_num_prims = shell_num_prims[i];
      BasisSet.shells[i].fprim = prim_pointers[i];
      BasisSet.shells[i].trans_vec = init_int_array(Symmetry.nirreps);
      for(j=0; j<Symmetry.nirreps; ++j)
    	  BasisSet.shells[i].trans_vec[j] = shell_trans_table[i][j];
      BasisSet.shells[i].fbf = shell_fbf[i];
      BasisSet.shells[i].fao = shell_fao[i];
      /*--- compute index of the stabilizer for the shell ---*/
      count = 1;
      for(g=1;g<Symmetry.nirreps;g++)
	if (i == BasisSet.shells[i].trans_vec[g]-1)
	  count++;
      stab_index = Symmetry.nirreps/count;
      if (Symmetry.max_stab_index < stab_index)
	Symmetry.max_stab_index = stab_index;
   }

   free(shell_center);
   free(shell_type);
   free(shell_num_prims);
   free(prim_pointers);
   free(shell_fbf);
   free(shell_fao);
   free_int_matrix(shell_trans_table);
   
   return;
}


void get_primitives(void)
{
   int i, j;
   double *exponents;      /* primitive gaussian exponents */
   double **ccoeffs;       /* primitive gaussian cont. coeff. for each ang. momentum*/

   /*--- read in exponents of primitive gaussians ---*/
   exponents = chkpt_rd_exps();

   /*--- read in coefficients of primitive gaussians ---*/
   ccoeffs = chkpt_rd_contr_full();

   /*--- allocate prims structure ---*/
   BasisSet.cgtos = (struct gaussian_function *)malloc(sizeof(struct gaussian_function)*BasisSet.num_prims);

   /*--- fill prims structure ---*/
   for (i=0; i<BasisSet.num_prims; i++){
     BasisSet.cgtos[i].exp = exponents[i];
     for(j=0;j<CINTS_MAX_AM;j++) 
       BasisSet.cgtos[i].ccoeff[j] = ccoeffs[i][j];
   }

   free(exponents);
   free_block(ccoeffs);

   return;
}


}}

