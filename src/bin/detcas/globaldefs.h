/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** GLOBALDEFS.H
**
** Global defines for DETCAS
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
*/

#ifndef _psi_src_bin_detcas_globaldefs_h
#define _psi_src_bin_detcas_globaldefs_h

namespace psi { namespace detcas {

#define MAX_RAS_SPACES 4
#define IOFF_MAX       50604
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define MAX_COMMENT 10

}} // end namespace psi::detcas

#endif // header guard

