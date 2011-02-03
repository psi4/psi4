#ifndef libscf_solver_integrator_H
#define libscf_solver_integrator_H
/*
 *  integrator.h
 *  Declaration of class Integrator for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/15/10.
 *
 */

#include "integrator_defines.h"
#include <libmints/gridblock.h>
#include <stdlib.h>
#include <vector>
#include <string>
using namespace psi;
using namespace boost;

namespace psi {

class Properties;
class BasisPoints;
class BasisSet;
class Molecule;
class Vector3;
class Options;

namespace scf {
/*! \ingroup SCF */
//! Integration Points/Weights container class
/*! Generates integration grids and computes points and weights
    for common molecular integration schemes
    References: Becke JCP 1988, Treutler+Alrichs JCP 1995
*/
class Integrator {
public:
    /** Factory Constructor, sets up Integrator associated with Molecule m,
    * with grid preferences/defaults in Options opt
    * @param m Molecule to be integrated over
    * @param opt Options object specifying grid type/make
    * @ return SharedIntegtor type
    */
    static boost::shared_ptr<Integrator> createIntegrator(shared_ptr<Molecule> m, Options & opt)
    {
        return boost::shared_ptr<Integrator>(new Integrator(m,opt));
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
    * @return true if all points have been asked for by getNextBlock, false otherwise
    */
    bool isDone() {return done_; }
    /** The next integration grid block
    * @return the GridBlock object containing up to the next block_size_
        points
    */
    SharedGridBlock getNextBlock();
    /** The integration grid block corresponding to block number
    * THREAD SAFE, allows a single Integrator to be used by many threads
    * @param block_number index of block number
    * @param block the PREALLOCATED GridBlock object containing up to the next block_size_
        points
    */
    void getBlock(unsigned long int block_number, SharedGridBlock block);
    /** The number of blocks in the integrator
    * @return the number of blocks in the Integrator
    **/
    int getNumberOfBlocks();
    /** The molecule associated with this Integrator
    * @return The molecule associated with this Integrator
    */
    shared_ptr<Molecule> getMolecule() const {return molecule_;}
    /** Number of radial points per atom
    * @return The number of radial points per nucleus
    */
    int getNRadial() const {return nradial_; }
    /** Number of spherical points per nucleus in the given sphere
    * @return The number of spherical points per nucleus
    * @param ind The index of the desired sphere
    */
    int getNSpherical(int ind = 0) const {return nspherical_[ind];}
    /** Number of spheres
    * @return The number of spheres in the pruned grid scheme
    * (either 1 or 5 at the moment)
    */
    int getNSpheres() const {return nspheres_;}
    /** Number of points per block
    * @return The number of points per block
    */
    int getBlockSize() const {return block_size_;}
    /** Selected nuclear scheme
    * @return Nuclear grid scheme, such as Treutler
    */
    nuclear_scheme getNuclearScheme() const {return nuclear_scheme_; }
    /** Selected radial scheme
    * @return Radial grid scheme, such as Treutler
    */
    radial_scheme getRadialScheme() const {return radial_scheme_; }
    /** Selected spherical scheme
    * @return Spherical grid scheme, such as Lebedev
    */
    spherical_scheme getSphericalScheme() const {return spherical_scheme_; }
    /** Selected special grid scheme
    * @return Special grid scheme, such as SG1
    */
    special_grid getSpecialGridScheme() const {return special_grid_; }
    /** String describing the grid setup
    * @return descriptive string
    */
    std::string getString();
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
    bool hasPrunedGrids() const {return nspheres_>1 ;}
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
    /// Chi values (Becke or Treutler)
    double** chi_values_;
    /// Inverse distance map between atoms
    double** inv_distance_map_;
};

typedef shared_ptr<Integrator> SharedIntegrator;

}}
#endif
