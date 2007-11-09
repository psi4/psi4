/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.1  2000/02/04 22:52:30  evaleev
 * Initial revision
 *
/* Revision 1.2  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: gprgid.cc 3592 2007-09-28 13:01:33Z evaleev $";

extern "C" char *gprgid()
{
   char *prgid = "SCF";

   return(prgid);
   }
