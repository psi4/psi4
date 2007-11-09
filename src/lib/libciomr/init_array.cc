/*
 \file init_array.cc
 \ingroup (CIOMR)
*/

#include <psifiles.h>
#include "includes.h"

extern "C" {

/*!
** init_array(): This function initializes an array of doubles of
** length 'size' and returns a pointer to the first element
** \ingroup (CIOMR)
*/
double * init_array(unsigned long int size)
   {
    double *array;

    if ((array = (double *) malloc(size*(unsigned long int)sizeof(double)))
        == NULL) {
       fprintf(stderr,"init_array:  trouble allocating memory \n");
       fprintf(stderr,"size = %ld\n",size);
       exit(PSI_RETURN_FAILURE);
       }
    bzero(array,size*(unsigned long int)sizeof(double));
    return(array);
   }

} /* extern "C" */
