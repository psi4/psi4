/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <math.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

/*-----------------------------------------------------------------------------------------------------------------
  This function builds SALCs of cartesian displacements
 -----------------------------------------------------------------------------------------------------------------*/

void build_cartdisp_salcs()
{
  int num_cd = 3*num_atoms;
  int i, j, u, xyz, irrep, G, cd;
  /* salc_symblk[i][j] is the pointer to j-th SALC in irrep i. Each SALC is an array of num_cd coefficients */
  double*** salc_symblk = (double***) malloc(sizeof(double**)*nirreps);
  for(i=0; i<nirreps; i++)
    salc_symblk[i] = (double**) malloc(sizeof(double*)*num_cd);
  
  /* work array */
  double* salc = init_array(num_cd);
  
  cdsalc_pi = init_int_array(nirreps);
  cdsalc2cd = block_matrix(num_cd,num_cd);
  
  /* compute SALCs */
  for(u=0; u<num_uniques; u++) {
    int atom = u2a[u];
    /* project each displacement */
    for(xyz=0; xyz<3; xyz++) {
      int cd = 3*atom + xyz;
      /* on each irrep */
      for(irrep=0; irrep<nirreps; irrep++) {
        memset((void*)salc,0,sizeof(double)*num_cd);
        
        /* this is the order of the atom stabilizer (subgroup which preserves atom's position unchanged) */
        int stab_order = 0;
        /* apply the projector: */
        for(G=0; G<nirreps; G++) {
          int Gatom = atom_orbit[atom][G];
          if (Gatom == atom)
            ++stab_order;
          int Gcd = 3*Gatom + xyz;
          double coeff = ao_type_transmat[1][G][xyz] * irr_char[irrep][G];
          salc[Gcd] += coeff;
        }
        /* normalize this SALC -- this is why stab_order was needed */
        for(cd=0; cd<num_cd; cd++)
          salc[cd] /= sqrt((double)nirreps*stab_order);
        
        /* if result is non-zero then add this salc to salc_symblk and increment salc counter */
        for(i=0; i<num_cd; i++) {
          if (fabs(salc[i])>1e-10 ) {
            salc_symblk[irrep][cdsalc_pi[irrep]] = init_array(num_cd);
            memcpy(salc_symblk[irrep][cdsalc_pi[irrep]],salc,sizeof(double)*num_cd);
            ++cdsalc_pi[irrep];
            break;
          }
        }
      }
    }
  }
  
  /* copy salc_symblk to cdsalc2cd */
  {
    int c = 0;
    for(irrep=0; irrep<nirreps; irrep++) {
      int num_per_irrep = cdsalc_pi[irrep];
      for(i=0; i<num_per_irrep; i++,c++) {
        for(j=0; j<num_cd; j++) {
          cdsalc2cd[j][c] = salc_symblk[irrep][i][j];
        }
        free(salc_symblk[irrep][i]);
      }
      free(salc_symblk[irrep]);
    }
    free(salc_symblk);
  }
  
  
  if (print_lvl >= PRINTUSOTAO) {
    fprintf(outfile,"    -Cartesian displacement SALCs per irrep:\n");
    fprintf(outfile,"    Irrep  #SALCs\n");
    fprintf(outfile,"    -----  ------\n");
    for(irrep=0;irrep<nirreps;irrep++) {
      fprintf(outfile,"    %3d    %4d\n",irrep,cdsalc_pi[irrep]);
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"    -Cartesian displacement SALCs:\n");
    print_mat(cdsalc2cd,num_cd,num_cd,outfile);
    fprintf(outfile,"\n");
  }
  
}

}} // namespace psi::input
