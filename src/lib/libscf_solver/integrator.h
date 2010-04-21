#ifndef libscf_solver_integrator_H
#define libscf_solver_integrator_H
/*
 *  integrator.h
 *  Definition of class Integrator for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/15/10.
 *
 */

#include <libpsio/psio.hpp>
#include <libmints/properties.h>
#include <libmints/basispoints.h>
#include <libmints/basisset.h>
#include <libmints/molecule.h>
#include <libmints/gridblock.h> 
#include <libmints/vector3.h>
#include "hf.h"
#include <stdlib.h>
#include <vector>
#include <string>
using namespace psi;

namespace psi { namespace scf {

/*! \ingroup SCF */
//! Integration Point/Weight container struct
// DEPRECATED 
typedef struct IntegrationPoint {
    Vector3 point;
    double weight;
};

/*! \ingroup SCF */
//! Lebedev Unit Sphere
typedef struct SphericalQuadrature {
    int n;
    double *x;
    double *y;
    double *z;
    double *w;
};

/*! \ingroup SCF */
//! Radial Quadrature
typedef struct RadialQuadrature {
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

/*! \ingroup SCF */
//! Integration Points/Weights container class
/*! Generates integration grids and computes points and weights
    for common molecular integration schemes
    References: Becke JCP 1988, Treutler+Alrichs JCP 1995
*/
class Integrator {
private:
    /// Grid Block
    SharedGridBlock grid_block_;
    /// Current Radial Quadrature
    RadialQuadrature rad_;
    /// Vector of Spherical Quadratures
    std::vector<SphericalQuadrature> sphere_;
    /// Current Atom
    int nuclearIndex_;
    /// Current Sphere
    int sphereIndex_;
    /// Current Spherical Point
    int sphericalIndex_;
    /// Current Radial Point
    int radialIndex_;
    /// Have any points been requested?
    bool started_;
protected:
    /// Returns the array of number of spherical grid points
    int *nspherical_;
    /// Returns the number of radial grid points per atom
    int nradial_;
    /// Returns the number of nuclei
    int nnuclei_;
    /// Number of spherical grids for pruned grids
    int nspheres_;
    /// Maximum number of grid points in a block
    int block_size_; 
    /// Selected nuclear weight scheme
    nuclear_scheme nuclear_scheme_;
    /// Selected spherical weight scheme
    spherical_scheme spherical_scheme_;
    /// Selected radial weight scheme
    radial_scheme radial_scheme_;
    /// Special grid scheme
    special_grid special_grid_;
    /// Returns the molecule this Integrator is associated with
    shared_ptr<Molecule> molecule_;
    /// Returns true if all points have been accessed, false otherwise
    bool done_;
    /// The nuclear weight at the current point, no heteroatomic size adjustment (see Becke 1988)
    double getNuclearWeightNaive(Vector3 v, int nuc);
    /// The nuclear weight at the current point, Becke size adjustment (see Becke 1988, appendix)
    double getNuclearWeightBecke(Vector3 v, int nuc);
    /// The nuclear weight at the current point, Treutler size adjustment (see Treutler 1995)
    double getNuclearWeightTreutler(Vector3 v, int nuc);
    /// The radial quadrature, with the Treutler weights (see Treutler 1995)
    RadialQuadrature getRadialQuadratureTreutler(int n, double xi, double alpha = 0.6);
    /// The radial quadrature, with the Becke weights (see Becke 1988)
    RadialQuadrature getRadialQuadratureBecke(int n, double xi);
    /// The radial quadrature, with the Euler-Maclaurin weights (see Gill 1993)
    RadialQuadrature getRadialQuadratureEulerMaclaurin(int n, double xi);
    /// Compute the Lebedev Sphere for the desired number of spherical points
    SphericalQuadrature getSphericalQuadrature(int deg);
    /// Compute the Lebedev reccurence points (used by getSphericalQuadrature)
    int getLebedevReccurencePoints(int type, int start, double a, double b, double v, SphericalQuadrature l);
    /// The adjusted Bragg-Slater radius of nucleus with charge n (See Becke 1988, etc)
    double getBraggSlaterRadius(int n);
    /// The optimized Treutler radius of nucleus with charge n (See Treutler 1995)
    double getTreutlerRadius(int n);
    /// The SG1 radius for SG1 grids (See Gill et. al., 1993)
    double getSG1Radius(int n);
    /// The SG1 cutoff alpha (See Gill et. al., 1993)
    double getSG1Alpha(int n, int index);
    /// Does this grid have pruning applied?
    const bool hasPrunedGrids() const {return nspheres_>1 ;}

public:
    /** Factory Constructor, sets up Integrator associated with Molecule m,
	* with grid preferences/defaults in Options opt
	* @param m Molecule to be integrated over  
	* @param opt Options object specifying grid type/make
	* @ return SharedIntegtor type
	*/
    static Integrator * createIntegrator(shared_ptr<Molecule> m, Options & opt)
    {
        return new Integrator(m,opt);
    }
    /** Constructor, sets up Integrator associated with Molecule m,
	* with grid preferences/defaults in Options opt
	* @param m Molecule to be integrated over  
	* @param opt Options object specifying grid type/make
	*/
    Integrator(shared_ptr<Molecule> m, Options & opt);
    /// Default Constructor, does nothing
    Integrator() {}
    /// Destructor
    ~Integrator();
    /// Checks that the sum of all coded Lebedev Grids is 1.0
    void checkLebedev(FILE* out);
    /// Checks that all Lebedev grids and radial integrators
    /// Are operational
    void checkSphericalIntegrators(FILE* out);
    /// Checks that all full molecular integrators are working
    void checkMolecularIntegrators(FILE* out);
    /// Prints available LebedevGrids
    void printAvailableLebedevGrids(FILE* out);
    /// Restarts the integration traverse
    void reset();
    /** Is the integration traverse completed
	* @return true if all points have been asked for by getNextPoints, false otherwise
	*/
    bool isDone() {return done_; }
    /** The next integration grid block
    * @ return the GridBlock object containing up to the next block_size_
        points
	*/
    SharedGridBlock getNextBlock();
    /** The next integration point and weight
	* @ return The IntegrationPoint struct corresponding to the next point and weight
    * DEPRECATED	
    */
    IntegrationPoint getNextPoint();
    /** The molecule associated with this Integrator
	* @return The molecule associated with this Integrator
	*/
    const shared_ptr<Molecule> getMolecule() const {return molecule_;}
    /** Number of radial points per atom
	* @return The number of radial points per nucleus
	*/
    const int getNRadial() const {return nradial_; }
    /** Number of spherical points per nucleus in the given sphere
	* @return The number of spherical points per nucleus
	* @param ind The index of the desired sphere
	*/
    const int getNSpherical(int ind = 0) const {return nspherical_[ind];}
    /** Number of spheres
	* @return The number of spheres in the pruned grid scheme
	* (either 1 or 5 at the moment)
	*/
    const int getNSpheres() const {return nspheres_;}
    /** Selected nuclear scheme
	* @return Nuclear grid scheme, such as Treutler
	*/
    const nuclear_scheme getNuclearScheme() const {return nuclear_scheme_; }
    /** Selected radial scheme
	* @return Radial grid scheme, such as Treutler
	*/
    const radial_scheme getRadialScheme() const {return radial_scheme_; }
    /** Selected spherical scheme
	* @return Spherical grid scheme, such as Lebedev
	*/
    const spherical_scheme getSphericalScheme() const {return spherical_scheme_; }
    /** Selected special grid scheme
	* @return Special grid scheme, such as SG1
	*/
    const special_grid getSpecialGridScheme() const {return special_grid_; }
    /** String describing the grid setup
	* @return descriptive string 
	*/
    string getString();
};

typedef shared_ptr<Integrator> SharedIntegrator;

}}
#endif
