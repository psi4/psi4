/*! \file
    \ingroup OPTKING
    \brief IRREP.CC :  symmetrizes set of delocalized internal coordinates
            also places irreps of coordinates in global array 'irr'

  &simples  -- address of an object of class internals
 **di_coord -- pointer to delocalized internal corrdinate matrix 
 (eigenvectors of G)
  returns:  -- pointer to 'symm_coordinate', a vector containing pointers to
               properly symmetrized delocalized internal coordinate vectors
significant modifications by J. Kenny June '00
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"

namespace psi { namespace optking {

static double **irrep_reduce(double **coord_mat, const simples_class &simples, int coords, int order);

double **irrep(const simples_class &simples, double **di_coord) {

  int i, j, k, a, ops, coord_num, irrep, offset,index,
  order,                  /* order of point group */
  num_spanned_big,        /* number of irreps spanned with significant characters*/    
  num_spanned,            /* number of irreps spanned */
  *irr_tmp,               /* if num_spanned > num_zero holds irrep info temporarily */
  *spanned_arr;           /* numbers of irreps spanned by coordinates */

  double norm, dot, 
  *tmp_evect,             /* effect of operation on coordinate */
  *evect_proj,            /* projection of coordinate onto irrep */
  **symm_coord,           /* holds pointers to properly symmetrizes coordinates */
  **coord_tmp, **coord_tmp2,
  **irreps_spanned;       /* irrep reduction coefficients, (coordinates,irreps) */
  char buffer[MAX_LINELENGTH], *err;

  /* allocate some memory */
  // symm_coord = (double **) malloc(num_nonzero * sizeof(double *));
  symm_coord = init_matrix(num_nonzero,simples.get_num());
  spanned_arr = init_int_array(num_nonzero);

  /* determine order of point group */
  order = 0;
  for(irrep=0;irrep<syminfo.nirreps;irrep++)
    order += ops_in_class[irrep]; 

  /*** find irrep reduction coefficients for delocalized coordinates ***/

  irreps_spanned = irrep_reduce(di_coord, simples, num_nonzero, order);

  /* count spanned irreps */
  num_spanned_big = 0;
  for(coord_num=0;coord_num<num_nonzero;++coord_num) {
    num_spanned = 0;
    for(irrep=0;irrep<syminfo.nirreps;++irrep) {
      if( irreps_spanned[coord_num][irrep] > EVAL_TOL )
        ++num_spanned;     /* slight contamination */
      if( irreps_spanned[coord_num][irrep] > SPANNED_IRREP_TOL )
        ++num_spanned_big; /* characters large enough to be projected */
    }
    spanned_arr[coord_num] = num_spanned;
  }

  /* count number of irreps spanned in total (reuse 'num_spanned') */ 
  num_spanned = 0;
  for(coord_num=0;coord_num<num_nonzero;coord_num++) {
    num_spanned += spanned_arr[coord_num];
  }

  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\n\nIrrep info for unsymmetrized coordinates (coordinate, irrep):\n");
    print_mat(irreps_spanned,num_nonzero,syminfo.nirreps,outfile);
    fprintf(outfile,"\n\nForming symmetrized coordinates.\n");
  }

  // if each coordinate spans only one irrep, no projection, form symm_coord and quit
  if(num_spanned == num_nonzero) {
    if(optinfo.print_delocalize == 1) {
      fprintf(outfile,"\nG Eigenvectors naturally spin-adapted.\n");
    }
    for(coord_num=0;coord_num<num_nonzero;coord_num++) { 
      for(irrep=0;irrep<syminfo.nirreps;++irrep) 
        if(irreps_spanned[coord_num][irrep] > SPANNED_IRREP_TOL)
          irr[coord_num] = irrep; 
      for (j=0; j<simples.get_num(); ++j)
        symm_coord[coord_num][j] = di_coord[coord_num][j];
    }
    free_matrix(irreps_spanned);
    free_int_array(spanned_arr);
    return symm_coord;
  }
  else if(num_spanned > num_nonzero) {
    // use projection operator for coordinates that span more than one irrep
    evect_proj = init_array(simples.get_num());
    tmp_evect = init_array(simples.get_num());
    irr_tmp = init_int_array(num_spanned_big);
    coord_tmp = init_matrix(num_spanned_big,simples.get_num()); 

    index = 0;     /* this variable keeps track of where to put next projection */
    for(coord_num=0;coord_num<num_nonzero;coord_num++) {

      /* again if only one irrep spanned, don't bother with projection operator */
      if(spanned_arr[coord_num] == 1) {
        for (i=0;i<simples.get_num();++i){
          coord_tmp[index][i] = di_coord[coord_num][i]; }
        ++index;
      }      

      /* projection operator time (where modifications are needed for higher symmetry) */

      else if(spanned_arr[coord_num] > 1) {

        /*loop over irreps*/
        for (irrep=0;irrep<syminfo.nirreps;irrep++) {
          zero_array(evect_proj,simples.get_num());

          /* if irrep spanned, project */
          if(irreps_spanned[coord_num][irrep] > SPANNED_IRREP_TOL) {
            /* loop over ops */
            for (ops = 0;ops<syminfo.nirreps;++ops) {
              zero_array(tmp_evect, simples.get_num());

              /* find effect of op on vector (Rvec) */
              for (j=0;j<simples.get_num();++j) {
                a = simples.id_to_index(syminfo.ict_ops[j][ops]);
                tmp_evect[a] += syminfo.ict_ops_sign[j][ops] * di_coord[coord_num][j];
              }  

              /* add char*Rvec to projection */
              for(j=0;j<simples.get_num();++j) {
                evect_proj[j] += syminfo.ct[irrep][ops]*tmp_evect[j];	       
              }
            }

            /* let's normalize now */
            norm=0;
            for(j=0;j<simples.get_num();++j) {
              norm += evect_proj[j] * evect_proj[j];
            }  
            for(j=0;j<simples.get_num();++j){
              coord_tmp[index][j] = evect_proj[j] / sqrt(norm);
            }
            ++index;
          }
        }
      }
    }

    free_array(tmp_evect);
    free_array(evect_proj);
    free_int_array(spanned_arr); 

    if (optinfo.print_delocalize == 1) {
      fprintf(outfile,"\n\nSymmetrized coordinates (each row is a coordinate):\n");    
      print_mat(coord_tmp,num_spanned_big,simples.get_num(),outfile);
    }

    /* check characters of projected coordinates to make sure they span only one irrep each */
    free_matrix(irreps_spanned); 
//    irreps_spanned = init_matrix(num_spanned_big,simples.get_num());

    irreps_spanned = irrep_reduce(coord_tmp, simples, num_spanned_big, order);

    /* count irreps spanned for error checking and put irrep info in 'irr_tmp' */
    num_spanned = 0;
    for (coord_num=0;coord_num<num_spanned_big;coord_num++) { 
      for (irrep=0;irrep<syminfo.nirreps;++irrep) {
        if ( irreps_spanned[coord_num][irrep] > SPANNED_IRREP_TOL ) {
          ++num_spanned;
          irr_tmp[coord_num] = irrep;
        }
      }
    }

    if (optinfo.print_delocalize == 1) {
      fprintf(outfile,"\nIrrep info for symmetrized coordinates: ");
      print_mat(irreps_spanned,num_spanned_big,syminfo.nirreps,outfile);
      fprintf(outfile,"\nOrthonormalizing symmetrized coordinates.\n");
    }

    // check that number of irreps spanned by projections is as expected
    if (num_spanned != num_spanned_big) {
      punt("projected coordinate(s) span incorrect number of irreps");
    }

    // orthogonalize the vectors
    offset=0;     /* offset + index tells where next coordinate goes in matrix */
    coord_tmp2 = init_matrix(num_spanned_big, simples.get_num());

    /* orthogonalize vectors for each symmetry block independently */
    for (irrep=0;irrep<syminfo.nirreps;irrep++) {
      index=0;

      /* include only if proper symmetry */
      for (coord_num=0;coord_num<num_spanned_big;++coord_num) {
        if (irr_tmp[coord_num] == irrep ) {
          if (index == 0) {
            /* since we only use pointer coord_tmp needs to stick around */
            coord_tmp2[offset + index] = coord_tmp[coord_num]; 
          }

          else if(index != 0 && irr_tmp[coord_num] == irrep) {
            tmp_evect = coord_tmp[coord_num]; 
            for (i=1;i<(index+1);++i) {
              dot=0.0;
              for (j=0;j<simples.get_num();++j) {
                dot += coord_tmp[coord_num][j] * coord_tmp2[offset+index-i][j];
              }			
              for (j=0;j<simples.get_num();++j) {
                tmp_evect[j] -= dot * coord_tmp2[offset+index-i][j];
              }
            }

            // find norm of evect_tmp, keep if big enough
            norm=0.0;
            for (i=0;i<simples.get_num();++i) {
              norm += tmp_evect[i] * tmp_evect[i];
            }
            norm = sqrt(norm);
            if (norm > 1.0E-2) {
              for (i=0;i<simples.get_num();++i) {
                coord_tmp2[offset+index][i]=tmp_evect[i]/norm;
              }
              ++index;
            }
          }  
          if(index==0) ++index; 
        }   
      }   
      offset += index;
    }

    free_int_array(irr_tmp); 

    /* check characters and orthogonality of orthogonalized coordinates */

    irreps_spanned = irrep_reduce(coord_tmp2, simples, num_spanned_big, order);

    zero_int_array(irr,num_nonzero);

    /* put final irrep info in global array 'irr' and check proper number of irreps spanned */
    num_spanned = 0;
    for (coord_num=0;coord_num<num_spanned_big;coord_num++) { 
      for (irrep=0;irrep<syminfo.nirreps;++irrep) {      
        if ( irreps_spanned[coord_num][irrep] > SPANNED_IRREP_TOL) {
          irr[coord_num] = irrep;
          ++num_spanned;
        }
      }
    }
    if (optinfo.print_delocalize == 1) {

      fprintf(outfile,"\nOrthonormalized coordinate matrix (there may be some rows of all zero's):");
      print_mat(coord_tmp2,num_spanned_big,simples.get_num(),outfile);
      fprintf(outfile,"\nIrrep info for orthonormalized coordinates: ");
      print_mat(irreps_spanned,num_spanned_big,syminfo.nirreps,outfile);
      fprintf(outfile,"\nChecking for orthogonality and proper number of coordinates.\n");

    }
    if (num_spanned != num_nonzero) {
      punt("orthogonalized coordinates span incorrect number of irreps");
    }

    /* check symmetry blocks, make sure orthogonal */
    for(irrep=0;irrep<syminfo.nirreps;++irrep) {
      for(i=0;i<num_spanned_big;++i) {
        for(j=(i+1);j<num_spanned_big;++j) {
          if (irr[i] == irr[j] == irrep) {
            dot = 0.0;
            for (k=0;k<simples.get_num();++k) {
              dot += coord_tmp2[i][k]*coord_tmp2[j][k];
            }
            if (dot > 1.0E-8 ) {
              fprintf(outfile,"\nwarning -- dot product of coord %i and coord %i: %lf",i,j,dot);
            }
          }
        }
      }
    }

    /* check that coord_tmp2 has 3n-6 nonzero rows (otherwise we don't have a proper basis for
       non-redundant subspace of primitives) */
    if (offset != num_nonzero) {
      punt("orthogonalization yields incorrect number of coordinates");
    }

    /* if everything is ok we can copy pointers and return */
    for(i=0;i<num_nonzero;++i) {
      for (j=0; j<simples.get_num(); ++j)
        symm_coord[i][j] = coord_tmp2[i][j];
    }

    fprintf(outfile,"\nSuccessful formation of symmetry adapted coordinates.\n");

    return symm_coord;
  }

  punt("num_spanned < num_nonzero, it should not happen");
  return NULL;
}


/*** IRREP_REDUCE finds reduction coefficients of delocalized coordinates
parameters: coord_mat  -- matrix with coordinates in rows 
&simples   -- address of object of class internals
num_coords -- number of coordinates
order      -- order of irrep
returns:    characters -- matrix of characters, (coordinate,irrep) ***/

double **irrep_reduce(double **coord_mat, const simples_class &simples, int num_coords, int order) {

  int coord, operations, irreducible, i, j, a;

  double red_coef,         /* reduction coefficient for irrep */
  *chars,                  /* character of current vector */
  *tmp_vec,                /* holds effect of operation on vector */
  **coef_mat;              /* matrix of reduction coefficients */

  chars = init_array(syminfo.nirreps);
  tmp_vec = init_array(simples.get_num());
  coef_mat = init_matrix(num_coords,syminfo.nirreps);

  /* loop over number of coordinate vectors */
  for(coord=0;coord<num_coords;++coord) { 

    zero_array(chars,syminfo.nirreps);

    /* determine the character of the vector under different operations */
    for (operations = 0;operations<syminfo.nirreps;++operations) {
      zero_array(tmp_vec, simples.get_num());
      for (i=0;i<simples.get_num();++i) {
        a = simples.id_to_index(syminfo.ict_ops[i][operations]);
        tmp_vec[a] += syminfo.ict_ops_sign[i][operations] * coord_mat[coord][i];
      }
      for (i=0;i<simples.get_num();++i) {
        chars[operations] += tmp_vec[i]*coord_mat[coord][i];
      }
    }

    /* determine the irrep of the vector */
    for (irreducible=0;irreducible<syminfo.nirreps;++irreducible) {
      red_coef=0.0;
      for (j=0;j<syminfo.nirreps;++j) {
        red_coef += chars[j] * syminfo.ct[irreducible][j] * ops_in_class[j];
      }
      red_coef = red_coef/ (float)order;
      coef_mat[coord][irreducible] = red_coef;
    }
  }

  free_array(chars); 
  free_array(tmp_vec);   
  return coef_mat;
}

}} /* namespace psi::optking */

