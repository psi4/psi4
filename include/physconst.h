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

/*
 * Make sure that you comment any new additions to this, so they are inlined into the manual
 * as well as the phyconst.py python file.  Note the format of the existing comment markers.
 */
#define _pi    3.14159265358979323846264338327950288 
#define _twopi 6.2831853071795862320E0
#define _h 6.62606896E-34                   /*- The Planck constant (Js) -*/
#define _c 2.99792458E8                     /*- Speed of light (ms$^{-1}$) -*/
#define _kb 1.3806504E-23                   /*- The Boltzmann constant (JK$^{-1}$) -*/
#define _R 8.314472                         /*- Universal gas constant (JK$^{-1}$mol$^{-1}$) -*/
#define _bohr2angstroms 0.52917720859       /*- Bohr to Angstroms conversion factor -*/
#define _bohr2m 0.52917720859E-10           /*- Bohr to meters conversion factor -*/
#define _bohr2cm 0.52917720859E-8           /*- Bohr to centimeters conversion factor -*/
#define _amu2g 1.660538782E-24              /*- Atomic mass units to grams conversion factor -*/
#define _amu2kg 1.660538782E-27             /*- Atomic mass units to kg conversion factor -*/
#define _au2amu 5.485799097E-4              /*- Atomic units (m$@@e$) to atomic mass units conversion factor -*/
#define _hartree2J 4.359744E-18             /*- Hartree to joule conversion factor -*/
#define _hartree2aJ 4.359744                /*- Hartree to attojoule (10$^{-18}$J) conversion factor -*/
#define _cal2J 4.184                        /*- Calorie to joule conversion factor -*/
#define _dipmom_au2si    8.47835281E-30     /*- Atomic units to SI units (Cm) conversion factor for dipoles -*/
#define _dipmom_au2debye 2.54174623         /*- Atomic units to Debye conversion factor for dipoles -*/
#define _dipmom_debye2si 3.335640952E-30    /*- Debye to SI units (Cm) conversion factor for dipoles -*/
#define _c_au 137.035999679                 /*- Speed of light in atomic units -*/
#define _hartree2ev 27.21138                /*- Hartree to eV conversion factor -*/
#define _hartree2wavenumbers 219474.6       /*- Hartree to cm$^{-1}$ conversion factor -*/
#define _hartree2kcalmol 627.5095           /*- Hartree to kcal mol$^{-1}$ conversion factor -*/
#define _hartree2MHz 6.579684E9             /*- Hartree to MHz conversion factor -*/
#define _kcalmol2wavenumbers 349.7551       /*- kcal mol$^{-1}$ to cm$^{-1}$ conversion factor -*/
#define _e0 8.854187817E-12                 /*- Vacuum permittivity (Fm$^{-1}$)-*/
#define _na 6.02214179E23                   /*- Avagadro's number -*/
#define _me 9.10938215E-31                  /*- Electron rest mass (in kg)-*/

/* For Cray X1 compilers */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* header guard */
