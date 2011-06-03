#ifndef three_index_cubature_H
#define three_index_cubature_H

#include <psi4-dec.h>
#include <psiconfig.h>
#include <libmints/vector3.h>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class GridBlock;
class Vector3;

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
    boost::shared_ptr<Matrix> rotation_;
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
              boost::shared_ptr<Matrix> rotation,
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
    /// Full x points. 
    double* x_;
    /// Full y points. 
    double* y_;
    /// Full z points. 
    double* z_;
    /// Full weights
    double* w_;
    
    /// Helpers for weights
    void applyNaiveWeights();
    void applyBeckeWeights();
    void applyTreutlerWeights();
    void applyStratmannWeights();

    // Helper for the first three types
    void applyStandardWeights(boost::shared_ptr<Matrix> chi);

public:
    MolecularGrid(boost::shared_ptr<Molecule> molecule);
    virtual ~MolecularGrid();   

    /// Build the grid 
    void buildGrid(std::vector<boost::shared_ptr<AtomicGrid> > atoms,
                   const std::string& scheme);

    /// Legacy method to GridBlock object. You do not own this.
    boost::shared_ptr<GridBlock> fullGrid();

    /// Number of grid points
    int npoints() const { return npoints_; } 
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


}
#endif
