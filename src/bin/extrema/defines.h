/*! \defgroup EXTREMA extrema: An optimizer (experimental code) */

/*###########################################################################*/
/*! 
   \file
   \ingroup EXTREMA
   \brief #define parameters 
*/
/*###########################################################################*/

#ifndef _psi_bin_extrema_defines_h_
#define _psi_bin_extrema_defines_h_

#define NORMAL_PRINT 1 /*!< normal print level */
#define RIDICULOUS_PRINT 3 /*!< level at which printing becomes ridiculous */
#define BT_CONV 1.0e-10 /*!< convergence level for internals->cartesians
			  iterative back transformation */
#define POS_NEG_TORS 1.0e-6 /*!< tolerance for pos/neg torsion pairs */
#define BT_LOOP 1000 /*!< max iterations for internals->cartesians 
		     iterative back transformation */
#define MAX_LINELENGTH 133 /*!< max length of lines in file11 */
#define EQUIV_GRAD 1.0e-6 /*!< tolerance for equivalent gradients */
#define ALMOST_ONE 1.0e-8 /*!< tolerance for numbers near 1.0 in 
			    cart_to_internals */
#define BOND_LIM 0.1 /*!< bond limit in angstroms */
#define ANGLE_LIM 5.0 /*!< angle limit in degrees */
#define DELOC_EV_TOL 1.0e-4 /*!< nonzero eigenvalue if G matrix tolerance */
#define BOND_TYPE 0 /*!< simple bond type id */
#define ANGLE_TYPE 1 /*!< simple angle type id */
#define TORS_TYPE 2 /*!< simple tors type id */
#define IRREP_TOL 0.05 /*!< smallest valid irrep coefficient */
#define ALMOST_ZERO 1.0e-14 /*!< used when computing norms */
#define ANOTHER_ZERO 1.0e-8 /*!< another zero */
#define NEAR_180 179.0 /*!< valence angle near 180? */
#define NOT_180 (180.0 - 1e-5) /*!< greater than this is 180 */

#define CART_TYPE 1 /*!< cartesian coordinate id */
#define ZMAT_TYPE 2 /*!< zmat coordinate id */ 
#define DELOC_TYPE 3 /*!< deloc coordinate id */

#endif // header guard
