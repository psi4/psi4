#ifndef libmints_integrator_def_H
#define libmints_integrator_def_H
/*
 *  integrator.h
 *  Declaration of class Integrator for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/15/10.
 *
 */
using namespace psi;

namespace psi {
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
/*! \ingroup SCF */
/// Enumeration defining available pruning modes
enum PruningScheme {constant_p, SG0_p, SG1_p, automatic_p};
/*! \ingroup SCF */
/// Enumeration defining available Voronoi closeness functions 
enum CoordinateScheme {elliptical_c, projection_c};
/*! \ingroup SCF */
/// Enumeration defining available point grouping functions 
enum GroupingScheme {boxes_g, voronoi_g};
/*! \ingroup SCF */
/// Enumeration defining available fuzzy voronoi functions 
enum VoronoiScheme {stratmann_v, becke_v};
/*! \ingroup SCF */
/// Enumeration defining available nuclear weights
enum NuclearScheme {naive_n, becke_n, treutler_n};
/*! \ingroup SCF */
/// Enumeration defining available spherical weights
enum SphericalScheme {lebedev_s, em_s};
/*! \ingroup SCF */
/// Enumeration defining available radial weights
enum RadialScheme {becke_r, treutler_r, em_r, multi_exp_r, mura_r};
/*! \ingroup SCF */
/// Enumeration defining special grids
enum SpecialScheme {SG1, SG0, NONE};
/*! \ingroup SCF */
//! Declaration of 1/ln(2)
const double INVLN2 = 1.0/log(2.0);

/// Treutler radial mapping radii (See Treutler 1995, Table 1, pp. 348)
const int maxTreutlerIndex_ = 36;
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
const int maxBraggSlaterIndex_ = 54;
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
const int maxSG1Index_ = 18;
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
const int maxSG0Index_ = 18;
const double SG0Radii_[] =
{
    2.0000, /* Ghost *///Ghost not defined
    1.3000, /* H */
    1.3000, /* He */// Noble gases not defined
    1.9500, /* Li */
    2.2000, /* Be */
    1.4500, /* B */
    1.2000, /* C */
    1.1000, /* N */
    1.1000, /* O */
    1.2000, /* F */
    1.2000, /* Ne */// Noble gases not defined
    2.3000, /* Na */
    2.2000, /* Mg */
    2.1000, /* Al */
    1.3000, /* Si */
    1.3000, /* P */
    1.1000, /* S */
    1.4500, /* Cl */
    1.4500 /* Ar */ // Noble gases not defined
};
const int SG0NRad_[] = 
{
    23, /* Ghost *///Ghost not defined
    23, /* H */
    23, /* He */// Noble gases not defined
    23, /* Li */
    23, /* Be */
    23, /* B */
    23, /* C */
    23, /* N */
    23, /* O */
    23, /* F */
    23, /* Ne */// Noble gases not defined
    26, /* Na */
    26, /* Mg */
    26, /* Al */
    26, /* Si */
    26, /* P */
    26, /* S */
    26, /* Cl */
    26 /* Ar */ // Noble gases not defined
};
const int SG0Orders_[] = 
{
      6,  6,  6,  6,  6,  6, 18, 18, 18, 26, 38, 74,110,146,146,146,146,146,146, 86, 50, 38, 18,  0,  0,  0, /* Ghost *///Ghost not defined
      6,  6,  6,  6,  6,  6, 18, 18, 18, 26, 38, 74,110,146,146,146,146,146,146, 86, 50, 38, 18,  0,  0,  0, /* H */
      6,  6,  6,  6,  6,  6, 18, 18, 18, 26, 38, 74,110,146,146,146,146,146,146, 86, 50, 38, 18,  0,  0,  0, /* He */// Noble gases not defined
      6,  6,  6,  6,  6,  6, 18, 18, 18, 26, 38, 74,110,146,146,146,146,146,146, 86, 50, 38, 18,  0,  0,  0, /* Li */
      6,  6,  6,  6, 18, 18, 26, 38, 38, 74, 86,110,110,146,146,146,146,146, 50, 38, 18,  6,  6,  0,  0,  0, /* Be */
      6,  6,  6,  6, 26, 26, 26, 26, 38, 38, 38, 86, 86, 86,146,146,146,146,146,146, 38,  6,  6,  0,  0,  0, /* B */
      6,  6,  6,  6,  6,  6, 18, 18, 26, 38, 38, 50, 50, 86,110,146,170,170,146,146, 86, 38, 18,  0,  0,  0, /* C */
      6,  6,  6,  6,  6,  6, 18, 18, 18, 26, 38, 38, 74, 74,110,170,170,146,146,146, 86, 50, 50,  0,  0,  0, /* N */
      6,  6,  6,  6,  6, 18, 26, 26, 38, 50, 50, 50, 50, 86,110,110,110,110,110, 86, 50, 38,  6,  0,  0,  0, /* O */
      6,  6,  6,  6, 38, 38, 50, 50, 50, 50, 74, 74,110,110,146,146,110,110, 86, 86, 86, 50,  6,  0,  0,  0, /* F */
      6,  6,  6,  6, 38, 38, 50, 50, 50, 50, 74, 74,110,110,146,146,110,110, 86, 86, 86, 50,  6,  0,  0,  0, /* Ne */// Noble gases not defined
      6,  6,  6,  6,  6,  6, 18, 18, 26, 26, 26, 38, 50, 50,110,110,110,110,110,110,110,110, 74, 74,  6,  6, /* Na */
      6,  6,  6,  6,  6, 18, 18, 26, 26, 38, 38, 50, 50, 74,110,110,146,146,146,146,110, 86, 38, 38, 18,  6, /* Mg */
      6,  6,  6,  6,  6,  6, 18, 18, 26, 38, 38, 50, 50, 74, 86,146,146,170,170,110,110, 86, 74, 26, 18,  6, /* Al */
      6,  6,  6,  6,  6, 18, 18, 18, 18, 38, 38, 38, 38, 50, 50, 50, 74,110,110,146,170,170,170, 86, 50,  6, /* Si */
      6,  6,  6,  6,  6, 18, 18, 18, 18, 38, 38, 38, 38, 50, 50, 50, 74,110,110,146,170,170,170, 86, 50,  6, /* P */
      6,  6,  6,  6, 18, 26, 26, 26, 26, 26, 26, 26, 26, 38, 38, 50, 74, 74,110,170,170,170,146,110, 50,  6, /* S */
      6,  6,  6,  6, 18, 18, 18, 18, 18, 18, 18, 26, 26, 38, 38, 50, 74,110,110,170,170,170,146,110, 86,  6, /* Cl */
      6,  6,  6,  6, 18, 18, 18, 18, 18, 18, 18, 26, 26, 38, 38, 50, 74,110,110,170,170,170,146,110, 86,  6 /* Ar */ // Noble gases not defined
};

}
#endif

