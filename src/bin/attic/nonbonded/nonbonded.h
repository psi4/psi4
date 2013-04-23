/*
 * nonbonded.h
 *
 * Contains empirical parameters relevant for nonbonded interactions
 *
 * C. David Sherrill
 * January 2008
 */

#ifndef _psi_include_masses_h
#define _psi_include_masses_h

#define LAST_VDW_RADII_GRIMME_INDEX (54)

/*
 * van der Waals radii (Angstrom) from atomic ROHF/TZV computations
 * used in the revised DFT-D method, B97-D
 * S. Grimme, J. Comput. Chem. 27, 1787-1799 (2006)
 */
double vdw_radii_grimme[] =
{
  2.000, /* Default element or ghost */
  1.001, /* H  */
  1.012, /* He */
  0.825, /* Li */
  1.408, /* Be */
  1.485, /* B  */
  1.452, /* C  */
  1.397, /* N  */
  1.342, /* O  */
  1.287, /* F  */
  1.243, /* Ne */
  1.144, /* Na */
  1.364, /* Mg */
  1.639, /* Al */
  1.716, /* Si */
  1.705, /* P  */
  1.683, /* S  */
  1.639, /* Cl */
  1.595, /* Ar */
  1.485, /* K  */
  1.474, /* Ca */
  1.562, /* Sc */
  1.562, /* Ti */
  1.562, /* V  */
  1.562, /* Cr */
  1.562, /* Mn */
  1.562, /* Fe */
  1.562, /* Co */
  1.562, /* Ni */
  1.562, /* Cu */
  1.562, /* Zn */
  1.650, /* Ga */
  1.727, /* Ge */
  1.760, /* As */
  1.771, /* Se */
  1.749, /* Br */
  1.727, /* Kr */
  1.628, /* Rb */
  1.606, /* Sr */
  1.639, /* Y  */
  1.639, /* Zr */
  1.639, /* Nb */
  1.639, /* Mo */
  1.639, /* Tc */
  1.639, /* Ru */
  1.639, /* Rh */
  1.639, /* Pd */
  1.639, /* Ag */
  1.639, /* Cd */
  1.672, /* In */
  1.804, /* Sn */
  1.881, /* Sb */
  1.892, /* Te */
  1.892, /* I  */
  1.881  /* Xe */
};


#define LAST_VDW_C6_GRIMME_INDEX (54)
/*
 * van der Waals C_6 parameters (in J nm^6 mol^-1) 
 * obtained for the revised DFT-D method, B97-D
 * S. Grimme, J. Comput. Chem. 27, 1787-1799 (2006)
 */
double vdw_C6_grimme[] =
{
  0.00, /* Default element or ghost */
  0.14, /* H  */
  0.08, /* He */
  1.61, /* Li */
  1.61, /* Be */
  3.13, /* B  */
  1.75, /* C  */
  1.23, /* N  */
  0.70, /* O  */
  0.75, /* F  */
  0.63, /* Ne */
  5.71, /* Na */
  5.71, /* Mg */
 10.79, /* Al */
  9.23, /* Si */
  7.84, /* P  */
  5.57, /* S  */
  5.07, /* Cl */
  4.61, /* Ar */
 10.80, /* K  */
 10.80, /* Ca */
 10.80, /* Sc */
 10.80, /* Ti */
 10.80, /* V  */
 10.80, /* Cr */
 10.80, /* Mn */
 10.80, /* Fe */
 10.80, /* Co */
 10.80, /* Ni */
 10.80, /* Cu */
 10.80, /* Zn */
 16.99, /* Ga */
 17.10, /* Ge */
 16.37, /* As */
 12.64, /* Se */
 12.47, /* Br */
 12.01, /* Kr */
 24.67, /* Rb */
 24.67, /* Sr */
 24.67, /* Y  */
 24.67, /* Zr */
 24.67, /* Nb */
 24.67, /* Mo */
 24.67, /* Tc */
 24.67, /* Ru */
 24.67, /* Rh */
 24.67, /* Pd */
 24.67, /* Ag */
 24.67, /* Cd */
 37.32, /* In */
 38.71, /* Sn */
 38.44, /* Sb */
 31.74, /* Te */
 31.50, /* I  */
 29.99  /* Xe */
};


#endif /* header guard */

