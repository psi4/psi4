/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
/*---------------------------------
  Definitions of general constants
 ---------------------------------*/

#define CHECKPOINTFILE 30               /*Checkpoint file number*/
#define NUM_ELEMENTS 108                /*Number of known elements*/
#define MAXATOM 300 
#define MAXIOFF 4000
#define MAXSTRING 132
#define MAXCONTRACTION 1500              /*Maximum number of primitives*/
#define MAXBASISCOLUMNS 2               /*Number of columns in basis_set[][]. Used to be equal to MAXCONTRACTION*/
#define ZERO (1.0E-11)                  /*Tolerance in the symmetry determination*/

/*------------------------
  Conditional compilation
 ------------------------*/
#define DEBUG 0


/*----------------
  Printing levels
 ----------------*/

#define DEBUGPRINT 5                    /*Dumps a lot of debugging information*/
#define PRINTSYMINFO 3                  /*Prints symmetry orbits, etc.*/
#define PRINTUSOTAO 4                   /*Prints unitary AO to SO matrix*/


/*-----------------------------------------------------
  Constants used in the geometry manipulation routines
 -----------------------------------------------------*/

#define ZERO_BOND_DISTANCE 1.0E-2	/*Cutoff on bond lengths - lesser values are forbidden*/
#define ZERO_BOND_ANGLE 1.0E-2		/*The same for bond angles*/
#define LINEAR_CUTOFF 1.0E-6		/*Cutoff used in the constructing of a linear fragment in the beginning of Z-matrix*/
#define ZERO_MOMENT_INERTIA 1.0E-10     /*Tolerance for degenerate rotational constants*/

/*-----------------
  Inline functions
 -----------------*/

#define sign(x) (((x) == 0) ? 0.0 : (((x) > 0) ? 1.0 : -1.0))  /*Returns 1.0 if int > 0, -1.0 if int < 0, and 0.0 otherwise*/
#define parity(m) ((m)%2 ? -1 : 1)  /*Returns (-1)^m */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/*--- dmalloc stuff to include into each source file ---*/
#ifdef DMALLOC
#include <dmalloc.h>
#endif






