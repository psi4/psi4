#ifndef libscf_solver_integrator_def_H
#define libscf_solver_integrator_def_H
/*
 *  integrator_.h
 *  Declaration of class Integrator for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/15/10.
 *
 */
using namespace psi;

namespace psi { namespace scf {
/*! \ingroup SCF */
//! Lebedev Unit Sphere
struct SphericalQuadrature {
    int n;
    double *x;
    double *y;
    double *z;
    double *w;
};

/*! \ingroup SCF */
//! Radial Quadrature
struct RadialQuadrature {
    int n;
    double *r;
    double *w;
};
/// Treutler radial mapping radii (See Treutler 1995, Table 1, pp. 348)
const double treutlerRadii_[] = 
{
    1.0, /* Ghost Atom */
    0.8, /* H */
    0.9, /* He */
    1.8, /* Li */
    1.4, /* Be */
    1.3, /* B */
    1.1, /* C */
    0.9, /* N */
    0.9, /* O */
    0.9, /* F */
    0.9, /* Ne */
    1.4, /* Na */
    1.3, /* Mg */
    1.3, /* Al */
    1.2, /* Si */
    1.1, /* P */
    1.0, /* S */
    1.0, /* Cl */
    1.0, /* Ar */
    1.5, /* K */
    1.4, /* Ca */
    1.3, /* Sc */
    1.2, /* Ti */
    1.2, /* V */
    1.2, /* Cr */
    1.2, /* Mn */
    1.2, /* Fe */
    1.2, /* Co */
    1.1, /* Ni */
    1.1, /* Cu */
    1.1, /* Zn */
    1.1, /* Ga */
    1.0, /* Ge */
    0.9, /* As */
    0.9, /* Se */
    0.9, /* Br */
    0.9/* Kr */
};
/*
 * Adjusted van der Waals radii (Angstrom) from atomic ROHF/TZV computations
 * temporarily acting as adjusted Bragg-Slater radii
 * S. Grimme, J. Comput. Chem. 27, 1787-1799 (2006)
 * TODO: Replace with official Bragg-Slater Radii
 */
const double braggSlaterRadii_[] =
{
    1.000, /* Default element or ghost */
    1.001, /* H  (untouched) */
    0.5*1.012, /* He */
    0.5*0.825, /* Li */
    0.5*1.408, /* Be */
    0.5*1.485, /* B  */
    0.5*1.452, /* C  */
    0.5*1.397, /* N  */
    0.5*1.342, /* O  */
    0.5*1.287, /* F  */
    0.5*1.243, /* Ne */
    0.5*1.144, /* Na */
    0.5*1.364, /* Mg */
    0.5*1.639, /* Al */
    0.5*1.716, /* Si */
    0.5*1.705, /* P  */
    0.5*1.683, /* S  */
    0.5*1.639, /* Cl */
    0.5*1.595, /* Ar */
    0.5*1.485, /* K  */
    0.5*1.474, /* Ca */
    0.5*1.562, /* Sc */
    0.5*1.562, /* Ti */
    0.5*1.562, /* V  */
    0.5*1.562, /* Cr */
    0.5*1.562, /* Mn */
    0.5*1.562, /* Fe */
    0.5*1.562, /* Co */
    0.5*1.562, /* Ni */
    0.5*1.562, /* Cu */
    0.5*1.562, /* Zn */
    0.5*1.650, /* Ga */
    0.5*1.727, /* Ge */
    0.5*1.760, /* As */
    0.5*1.771, /* Se */
    0.5*1.749, /* Br */
    0.5*1.727, /* Kr */
    0.5*1.628, /* Rb */
    0.5*1.606, /* Sr */
    0.5*1.639, /* Y  */
    0.5*1.639, /* Zr */
    0.5*1.639, /* Nb */
    0.5*1.639, /* Mo */
    0.5*1.639, /* Tc */
    0.5*1.639, /* Ru */
    0.5*1.639, /* Rh */
    0.5*1.639, /* Pd */
    0.5*1.639, /* Ag */
    0.5*1.639, /* Cd */
    0.5*1.672, /* In */
    0.5*1.804, /* Sn */
    0.5*1.881, /* Sb */
    0.5*1.892, /* Te */
    0.5*1.892, /* I  */
    0.5*1.881  /* Xe */
};
/*! \ingroup SCF */
/*
*Declaration of SG1_radii for use with SG1 grid
*Gill, M.W., Johnson, B.G., Pople, J.A., Chem. Phys. Lett.,
*209, July 1993, pp. 506 
*See Table 1, pp. 508
* TODO: Add more radii (Ar just wont cut it)
*/
const double SG1Radii_[] =
{
    2.0000, /* Ghost */
    1.0000, /* H */
    0.5882, /* He */
    3.0769, /* Li */
    2.0513, /* Be */
    1.5385, /* B */
    1.2308, /* C */
    1.0256, /* N */
    0.8791, /* O */
    0.7692, /* F */
    0.6838, /* Ne */
    4.0909, /* Na */
    3.1579, /* Mg */
    2.5714, /* Al */
    2.1687, /* Si */
    1.8750, /* P */
    1.6514, /* S */
    1.4754, /* Cl */
    1.3333 /* Ar */
};
/*! \ingroup SCF */
/*
*Declaration of SG1 spherical grid orders used
*Gill, M.W., Johnson, B.G., Pople, J.A., Chem. Phys. Lett.,
*209, July 1993, pp. 506 
*See pp. 509, right column
*/
const int SG1GridOrders_[] =
{
    6,
    38,
    86,
    194,
    86
};
/*! \ingroup SCF */
/*
* Declaration of SG1 sphere cutoffs
* \alpha in \alpha R
* Gill, M.W., Johnson, B.G., Pople, J.A., Chem. Phys. Lett.,
* 209, July 1993, pp. 506 
* See Table 4 pp. 509
* TODO: Add higher atoms
*/
const double SG1Alpha_[] =
{
    0.25, 0.5, 1.0, 4.5,/* H-He */
    0.1667, 0.50, 0.9, 3.5,/* Li-Ne */
    0.10, 0.4, 0.8, 2.5,/* Na-Ar */
};
/*! \ingroup SCF */
/// Enumeration defining available nuclear weights
enum nuclear_scheme {naive_n, becke_n, treutler_n};
/*! \ingroup SCF */
/// Enumeration defining available spherical weights
enum spherical_scheme {lebedev_s};
/*! \ingroup SCF */
/// Enumeration defining available radial weights
enum radial_scheme {becke_r, treutler_r, em_r};
/*! \ingroup SCF */
/// Enumeration defining special grids
enum special_grid {SG1, NONE};
/*! \ingroup SCF */
//! Declaration of 1/ln(2)
const double INVLN2 = 1.0/log(2.0);
}}
#endif

