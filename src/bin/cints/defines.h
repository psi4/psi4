#ifndef _psi_src_bin_cints_defines_h
#define _psi_src_bin_cints_defines_h

/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
/*---------------------
  Customizable values
 ---------------------*/
#define CINTS_MAX_AM 10            /* Max ang. mom. + 1 */
#define MAX_NUM_AO 4000            /* Maximum number of AOs */
#define CUTOFF 15                  /* Default cutoff on the integrals */
#define ROT_INV_TOLER 1E-4         /* Tolerance on the "rotational variance" */
#define MAX_NUM_DOUBLES 2500000    /* Default number of double words to use in MP2 */

/*----------------------------------
  Flags for conditional compilation
 ----------------------------------*/
#define PRINT 0          /* Print two-electron integrals out? */
#define PRINT_DERIV1 0   /* Print derivative integrals? Do this *only* for 
			    no-symmetry, non-puream cases */
#define USE_MM 1         /* Use matrix multiplies rather than just nested loops
                            in cart->puream transformation and normalization (should always be
			    set to 1 except for testing purposes)*/
#define SPARSE_C2P 1     /* Use sparcity of cartesian to spherical harmonics transformation */
#define USE_BLAS 0       /* Use routines from vendor-provided BLAS library for
			    4-index transformations such as cart->puream
			    (only if SPARSE_C2P is set 0, otherwise sparse matrix
			    multiplies will be used), AO->MO, etc.;
			    Use only if you have libblas.a or its other analog available */
#define SCF_ONLY 0       /* If you want to be able to compute integrals needed in
			    SCF only if WFN=SCF - set this to 1 */

/*------------------------------------------------------------
  Predefined cutoffs used while computing Fm(T) and integrals
 ------------------------------------------------------------*/
#define SOFT_ZERO 1E-5             /* Definition of a soft floating-point "zero" */
#define ZERO 1E-15                 /* Definition of a hard floating-point "zero" */
#define EPS 1.0e-17                /* Another definition of floating-point "zero"
				      used in computing auxiliary function */
#undef USE_TAYLOR_FM              /* Use Taylor interpolation formula to compute Fm(T) */
#define TAYLOR_ORDER 6             /* Order of Taylor interpolation used to compute Fm(T) */

/*----------------------------------
  Thresholds for printing out stuff
 ----------------------------------*/
#define PRINT_DEBUG 5    /* Print everything */
#define PRINT_INTRO 1    /* Print level to print intro overhead */
#define PRINT_OPTIONS 1  /* Print level to print options out */
#define PRINT_BASIS 3    /* Print level to print basis set information */
#define PRINT_GEOMETRY 2 /* Print level to print cartesian geometry */
#define PRINT_CCOEFF 4   /* Print level to print coupling coefficients */
#define PRINT_OPDM 4     /* Print level to print onepdms */
#define PRINT_OEI 3      /* Print level to print one-electron integrals */
#define PRINT_OEDERIV 2      /* Print level to print oe and nuclear contrinution to gradient */
#define PRINT_TEDERIV 2      /* Print level to print te contribution to gradient  */
#define PRINT_MOINFO_CORR 1  /* Print level to print orbital information for correlated calculations */

/*-------------------------
  Default sizes for arrays
 -------------------------*/
#define MAXFACT 100
#define MAXNIRREPS 8
#define IOFFMAX MAX_NUM_AO

/*----------------
  Macro functions
 ----------------*/
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define INDEX(a,b) (((a)>(b)) ? ioff[(a)] + (b) : ioff[(b)] + (a))
#define SWAP(a,b) {dum = (a); (a) = (b); (b) = dum;}

/*------
  Misc.
 ------*/
#define NUM_TE_TYPES 4   /* Number of types of integrals to compute
			    0 - ERIs
			    1 - r12
			    2 - [r12,T1]
			    3 - [r12,T2]
			    */

#define MkPT2_USE_IWL 1 /* How the MkPT2 routine dumps integrals to disk */
#define MkPT2_TEST 0    /* Whether to test the MkPT2 integrals by printing them on writing to IWL and reading them back in */
/*----------
  DFT cutoffs
  ----------*/

#define WEIGHT_CUTOFF 1E-10  /* if the weighting function is lower
				 than this, do not contract into the 
				 Fock Matrix
			      */

#define DEN_CUTOFF 1E-10     /* If density at a point is less than
				this, do not contract into the Fock
				matrix 
			     */

#endif
