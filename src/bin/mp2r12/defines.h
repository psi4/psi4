/*! \file
    \ingroup MP2R12
    \brief Enter brief description of file here 
*/
/*--- Useful and not-so-useful macro definitions ---*/
#ifndef _psi_src_bin_mp2r12_h_
#define _psi_src_bin_mp2r12_h_


#define USE_BLAS 1

#define IOFF 32641
#define MAXIOFF3 255
#define MAX_NUM_IRREPS 8
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define MIN(i,j) (((i) >= (j)) ? (j) : (i))

#endif /* Header guard */
