#ifndef libmints_cubature_H
#define libmints_cubature_H

#include <psi4-dec.h>
#include <psiconfig.h>
#include <libmints/vector3.h>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class GridBlock; // Deprecated
class Vector3;
class BasisExtents;
class BlockOPoints;

class SphericalGrid {

protected:
    /// Master map of order to number of points 
    static std::map<int, int> lebedev_order_to_points_;
    /// Helper function to generate lebedev reccurence points
    static int lebedevReccurence(int, int, double, double, double, SphericalGrid* leb);
    /// Helper to build Lebedev grid
    static boost::shared_ptr<SphericalGrid> buildLebedevGrid(int order);
    /// Static method to build and store all Lebedev grids (idempotent)
    static void buildAllLebedevGrids();
    /// Static map of all Lebedev grids to prevent regeneration 
    static std::map<int, boost::shared_ptr<SphericalGrid> > lebedev_grids_;

    /// Scheme
    std::string scheme_;
    /// Number of points
    int npoints_;
    /// Spherical harmonic order
    int order_;
    /// x coordinates on unit sphere 
    double* x_;
    /// y coordinates on unit sphere 
    double* y_;
    /// z coordinates on unit sphere 
    double* z_;
    /// weights  on unit sphere (normalized to 4\pi)
    double* w_;

    /// Lookup method for Lebedev grids
    static boost::shared_ptr<SphericalGrid> lookupLebedevGrid(int order);

public:
    SphericalGrid();
    virtual ~SphericalGrid();   
    static boost::shared_ptr<SphericalGrid> buildGrid(const std::string& scheme, int order);
    static std::map<int, int> lebedevOrdersToPoints();

    /// Number of grid points
    int npoints() const { return npoints_; } 
    /// Order of quadrature 
    int order() const { return order_; } 
    /// Scheme used to build the quadrature
    std::string scheme() const { return scheme_; }
    /// Print a trace of this spherical quadrature
    virtual void print(FILE* out = outfile); 
    
    /// The x points. You do not own this 
    double* x() const { return x_; }
    /// The y points. You do not own this 
    double* y() const { return y_; }
    /// The z points. You do not own this 
    double* z() const { return z_; }
    /// The weights, normalized to 4\pi on Omega. You do not own this 
    double* w() const { return w_; }

};

class RadialGrid {

protected:
    /// Scheme
    std::string scheme_;
    /// Number of points
    int npoints_;
    /// r coordinates
    double* r_;
    /// weights without r^2 (1D)
    double* w_x_;
    /// weights with r^2 (3D)
    double* w_r_;

    static boost::shared_ptr<RadialGrid> buildTreutlerGrid(int n, double xi, double alpha = 0.6);
    static boost::shared_ptr<RadialGrid> buildBeckeGrid(int n, double xi);
    static boost::shared_ptr<RadialGrid> buildMuraGrid(int n, double xi);
    static boost::shared_ptr<RadialGrid> buildEMGrid(int n, double xi);
    static boost::shared_ptr<RadialGrid> buildMultiExpGrid(int n, double xi);
public:
    RadialGrid();
    virtual ~RadialGrid();   

    /// Build the radial grid
    static boost::shared_ptr<RadialGrid> buildGrid(const std::string& scheme, int n, double xi);

    /// Number of grid points
    int npoints() const { return npoints_; } 
    /// Scheme used to build the quadrature
    std::string scheme() const { return scheme_; }
    /// Print a trace of this radial quadrature
    virtual void print(FILE* out = outfile); 
    
    /// The r points. You do not own this 
    double* r() const { return r_; }
    /// The weights, without the r^2 element. You do not own this
    double* w_x() const { return w_x_; }
    /// The weights, with the r^2 element. You do not own this
    double* w_r() const { return w_r_; }
         
};

class AtomicGrid {

protected:
    /// The radial quadrature
    boost::shared_ptr<RadialGrid> radial_;
    /// The vector of spherical grids (unit spheres)  
    std::vector<boost::shared_ptr<SphericalGrid> > sphericals_;
    /// The rotation matrix to bring this grid to standard orientation
    SharedMatrix rotation_;
    /// The center of this atom
    Vector3 center_;

    /// Total points for this atom
    int npoints_;
    /// Full x points. 
    double* x_;
    /// Full y points. 
    double* y_;
    /// Full z points. 
    double* z_;
    /// Full weights
    double* w_;

public:
    AtomicGrid();
    virtual ~AtomicGrid();   

    void buildGrid(const Vector3& center,
              SharedMatrix rotation,
              boost::shared_ptr<RadialGrid> radial,          
              std::vector<boost::shared_ptr<SphericalGrid> > spheres); 
    
    /// Number of grid points
    int npoints() const { return npoints_; } 
    /// Print a trace of this atomic quadrature
    virtual void print(FILE* out = outfile); 
    
    /// The x points. You do not own this 
    double* x() const { return x_; }
    /// The y points. You do not own this 
    double* y() const { return y_; }
    /// The z points. You do not own this 
    double* z() const { return z_; }
    /// The weights, normalized to 1 on R3. You do not own this 
    double* w() const { return w_; }

    /// Get the radial grid
    boost::shared_ptr<RadialGrid> radial() const { return radial_; }
    /// Get the spherical grid vector
    std::vector<boost::shared_ptr<SphericalGrid> > spherical() const { return sphericals_; } 
};

class MolecularGrid {

protected:
    /// Bragg-Slater Radii
    static std::vector<double> BSRadii_;
    /// Helper to get BS Radii
    static void getBSRadii(); 

    /// Scheme
    std::string scheme_;
    /// The vector of atomic quadratures
    std::vector<boost::shared_ptr<AtomicGrid> > atoms_;
    /// The molecule this grid is built on
    boost::shared_ptr<Molecule> molecule_;
    
    /// Total points for this molecule 
    int npoints_;
    /// Maximum number of points in a block
    int max_points_;
    /// Maximum number of functions in a block
    int max_functions_;
    /// Full x points. 
    double* x_;
    /// Full y points. 
    double* y_;
    /// Full z points. 
    double* z_;
    /// Full weights
    double* w_;
 
    /// Vector of blocks 
    std::vector<boost::shared_ptr<BlockOPoints> > blocks_;
   
    /// Points to basis extents, built internally
    boost::shared_ptr<BasisExtents> extents_;
    /// BasisSet from extents_
    boost::shared_ptr<BasisSet> primary_;
 
    /// Helpers for weights
    void applyNaiveWeights();
    void applyBeckeWeights();
    void applyTreutlerWeights();
    void applyStratmannWeights();

    /// Helper for the first three types
    void applyStandardWeights(SharedMatrix chi);

    /// Sieve and block
    void sieve();    
    void remove_zero_points();
    void remove_distant_points(double Rcut);
    void block(int max_points, int min_points);

public:
    MolecularGrid(boost::shared_ptr<Molecule> molecule);
    virtual ~MolecularGrid();   

    /// Helper function to produce standard grid orientation for molecule mol
    SharedMatrix standard_orientation(boost::shared_ptr<Molecule> mol);
    /// Build the grid 
    void buildGrid(std::vector<boost::shared_ptr<AtomicGrid> >& atoms, const std::string& nuclear_scheme, 
        boost::shared_ptr<BasisExtents> extents, int max_points, int min_points);

    /// Legacy method to GridBlock object. You do not own this.
    boost::shared_ptr<GridBlock> fullGrid();

    /// Number of grid points
    int npoints() const { return npoints_; } 
    /// Maximum number of grid points in a block
    int max_points() const { return max_points_; }
    /// Maximum number of funtions in a block 
    int max_functions() const { return max_functions_; }
    /// Print a trace of this molecular quadrature
    virtual void print(FILE* out = outfile, int print = 2); 
 
    /// The x points. You do not own this 
    double* x() const { return x_; }
    /// The y points. You do not own this 
    double* y() const { return y_; }
    /// The z points. You do not own this 
    double* z() const { return z_; }
    /// The weights, normalized to 1 on R3. You do not own this 
    double* w() const { return w_; }

    /// Pointer to basis extents
    boost::shared_ptr<BasisExtents> extents() const { return extents_; }
    /// Set of spatially sieved blocks of points, generated by sieve() internally 
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks() const { return blocks_; }

};

class PseudospectralGrid : public MolecularGrid {

protected:
    /// The primary basis
    boost::shared_ptr<BasisSet> primary_;
    /// The dealias basis 
    boost::shared_ptr<BasisSet> dealias_;
    /// The Options object
    Options& options_;  
    /// The filename used to optionally build the grid
    std::string filename_; 

    /// Master builder methods
    void buildGridFromOptions();
    void buildGridFromFile(); 

public:

    /// Constructor to use for autogeneration
    PseudospectralGrid(boost::shared_ptr<Molecule> molecule,
                       boost::shared_ptr<BasisSet> primary, 
                       boost::shared_ptr<BasisSet> dealias,
                       Options& options); 
    /// Construtor to use for semiautomatic generation with grid file
    PseudospectralGrid(boost::shared_ptr<Molecule> molecule,
                       boost::shared_ptr<BasisSet> primary, 
                       boost::shared_ptr<BasisSet> dealias,
                       const std::string& filename, 
                       Options& options); 
    virtual ~PseudospectralGrid();

};

class DFTGrid : public MolecularGrid {

protected:
    /// The primary basis 
    boost::shared_ptr<BasisSet> primary_; 
    /// The Options object
    Options& options_;  
    /// Master builder methods
    void buildGridFromOptions();

public:
    DFTGrid(boost::shared_ptr<Molecule> molecule,
            boost::shared_ptr<BasisSet> primary,
            Options& options);
    virtual ~DFTGrid();
};

class BlockOPoints {

protected:
    /// number of points in this block
    int npoints_;
    /// Pointer to x (does not own)
    double* x_;
    /// Pointer to y (does not own)
    double* y_;
    /// Pointer to z (does not own)
    double* z_;
    /// Pointer to w (does not own)
    double* w_;
    /// Relevant shells, local -> global 
    std::vector<int> shells_local_to_global_;
    /// Relevant functions, local -> global 
    std::vector<int> functions_local_to_global_;
    /// Reference to the extents object
    boost::shared_ptr<BasisExtents> extents_;

    /// Center of this BlockOPoints
    Vector3 xc_;
    /// Bounding radius of the BlockOPoints
    double R_;

    /// Populate significant functions given information in extents
    void populate();
    /// Compute bounding sphere
    void bound();

public:
    BlockOPoints(int npoints, double* x, double* y, double* z, double* w, 
        boost::shared_ptr<BasisExtents> extents);     
    virtual ~BlockOPoints();

    /// Refresh populations (if extents_->delta() changes)
    void refresh() { populate(); }

    /// Number of grid points
    int npoints() const { return npoints_; } 
    /// Print a trace of this BlockOPoints
    void print(FILE* out = outfile, int print = 2); 
    
    /// The x points. You do not own this 
    double* x() const { return x_; }
    /// The y points. You do not own this 
    double* y() const { return y_; }
    /// The z points. You do not own this 
    double* z() const { return z_; }
    /// The weights. You do not own this 
    double* w() const { return w_; }

    /// Relevant shells, local -> global 
    const std::vector<int>& shells_local_to_global() const { return shells_local_to_global_; }
    /// Relevant functions, local -> global 
    const std::vector<int>& functions_local_to_global() const { return functions_local_to_global_; }
};

class BasisExtents {

protected: 
    /// Basis this corresponds to
    boost::shared_ptr<BasisSet> primary_;
    /// Cutoff value for basis values 
    double delta_;
    /// Significant extent of shells
    boost::shared_ptr<Vector> shell_extents_;
    /// Maximum extent
    double maxR_;   
 
    /// Recompute and shell_extents_
    void computeExtents();
public:
    BasisExtents(boost::shared_ptr<BasisSet> primary, double delta);
    virtual ~BasisExtents(); 

    /// Print a trace of these extents
    void print(FILE* out = outfile);
    /// Reset delta and recompute extents
    void set_delta(double delta) { delta_ = delta; computeExtents(); }   

    /// The cutoff value 
    double delta() const { return delta_; } 
    /// The basis set this BasisExtents is built on
    boost::shared_ptr<BasisSet> basis() const { return primary_; }
    /// WCS significant extent of each shell
    boost::shared_ptr<Vector> shell_extents() const { return shell_extents_; }
    /// Maximum spatial extent over all atoms
    double maxR() const { return maxR_; }
};

}
#endif
