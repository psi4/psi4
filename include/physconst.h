/*
** PHYSCONST.H : Includes physical constants
** D. Sherrill, April 1993
** Edward F. Valeev, March 1998
*/

/* revised 5/7/93 to agree w/ 1986 CODATA
 * recommended values J. Phys. Chem. Ref. Data 17, 1795 (1988)
 * all values revised
 *
 * revised 4/28/02 to agree with NIST online database
 * updated au2amu, hartree2J, hartree2wavenumbers
 *
 * revised 10/28/02 by TDC to include vacuum permittivity (_e0)
 * Avagadro's number (_na).
 * 
 * Added electron rest mass from NIST database 6/27/03.
 * -TDC
 */

#ifndef _psi_include_physconst_h_
#define _psi_include_physconst_h_

#define _pi    3.14159265358979323846264338327950288
#define _twopi 6.2831853071795862320E0
#define _h 6.6260755E-34  /* J s    */
#define _c 2.99792458E8  /* m s^-1 */
#define _kb 1.380658E-23  /* J K^-1 */
#define _psi3_R 8.314510        /* J K^-1 mol^-1 */
#define _bohr2angstroms 0.529177249
#define _bohr2m 0.529177249E-10
#define _bohr2cm 0.529177249E-8
#define _amu2g 1.6605402E-24
#define _amu2kg 1.6605402E-27
#define _au2amu 5.485799110E-4 
#define _hartree2J 4.35974381E-18
#define _hartree2aJ 4.35974381
#define _cal2J 4.184
#define _dipmom_au2si    8.47835791E-30
#define _dipmom_au2debye 2.54175
#define _dipmom_debye2si 3.33564E-30
#define _c_au 137.0359895
#define _hartree2ev 27.211396
#define _hartree2wavenumbers 219474.6313710
#define _e0 8.854187816E-12  /* vacuum permittivity */
#define _na 6.022136736E23 /* Avagadro's number */
#define _me 9.10938188E-31 /* electron rest mass */

/* For Cray X1 compilers */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* header guard */
