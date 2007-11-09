#include <psifiles.h>
#include "includes.h"

extern "C" {

double *** block_mat_alloc(int n_so_typs, int num_ir, int* num_so)
   {
      double ***array;
      int i,blk,nn;

      if ((array = (double ***) malloc(sizeof(double **)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_mat_alloc\n");
          exit(PSI_RETURN_FAILURE);
          }

      for (i=0,blk=0; i < num_ir ; i++) {
         if (nn=num_so[i]) {
            array[blk] = (double **) init_matrix(nn,nn);
            blk++;
            }
         }

      return(array);
      }

void block_mat_dealloc(double*** array, int num_ir, int* num_so)
   {
      int i;
      int blk=0;

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            free_matrix(array[blk],num_so[i]);
            blk++;
            }
         }

      free(array);
      }

double ** block_arr_alloc(int n_so_typs, int num_ir, int* num_so)
   {
      double **array;
      int i;
      int j=0;

      if ((array = (double **) malloc(sizeof(double *)*n_so_typs))==NULL) {
          fprintf(stderr,"trouble in block_arr_alloc\n");
          exit(PSI_RETURN_FAILURE);
          }

      for (i=0; i < num_ir ; i++) {
         if (num_so[i]) {
            int nget = num_so[i]*(num_so[i]+1)/2;
            if (j>=n_so_typs) {
                fprintf(stderr,"block_arr_alloc: too many rows\n");
                exit(PSI_RETURN_FAILURE);
              }
            array[j] = (double *) init_array(nget);
            j++;
            }
         }
      return(array);
      }

void block_arr_dealloc(double** array, int n_so_typs)
   {
      int i;

      for (i=0; i < n_so_typs ; i++) {
         free(array[i]);
         }

      free(array);
      }

} /* extern "C" */
