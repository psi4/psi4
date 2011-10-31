#ifndef _psi_src_lib_libmints_extern_potential_h_
#define _psi_src_lib_libmints_extern_potential_h_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

namespace psi {

class Matrix;
class BasisSet;

/*! \ingroup MINTS
 *  \class ExternalPotential 
 *  Stores external potential field, computes external potential matrix
 *  Like standard potential integrals, this is negative definite (electrons are the test charge)
 */
class ExternalPotential {
protected:
    /// Name of potential
    std::string name_;
    /// <Z,x,y,z> array of charges
    std::vector<boost::tuple<double,double,double,double> > charges_;
    /// <mx,my,mz,x,y,z> array of dipoles
    std::vector<boost::tuple<double,double,double,double,double,double> > dipoles_;
    /// <qxx,qxy,qxz,qyy,qyz,qzz,x,y,z> array of quadrupoles
    std::vector<boost::tuple<double,double,double,double,double,double,double,double,double> > quadrupoles_;

public:
    /// Constructur, does nothing
    ExternalPotential();
    /// Destructor, does nothing
    ~ExternalPotential();

    /// Set name 
    void setName(const std::string & name) { name_ = name; }    

    /// Add a charge Z at (x,y,z)
    void addCharge(double Z,double x, double y, double z);
    /// Add a dipole (mx,my,mz) at (x,y,z)
    void addDipole(double mx, double my, double mz, double x, double y, double z);
    /// Add a quadrupole (qxx,qxy,qxz,qyy,qyz,qzz) at (x,y,z)
    void addQuadrupole(double qxx, double qxy, double qxz, double qyy, double qyz,
         double qzz, double x, double y, double z);

    /// Reset the field to zero (eliminates all entries) 
    void clear();
    
    /// Translate the origin by (dx,dy,dz)
    void translate(double dx, double dy, double dz);   
    /// Rotate about the origin by the rotation matrix R
    void rotate(SharedMatrix R); 

    /// Compute the external potential matrix in the given basis set
    /// C1 for now!
    SharedMatrix computePotentialMatrix(boost::shared_ptr<BasisSet> basis);
    /// Compute the external potential at a single point
    double computePotentialPoint(double x, double y, double z);

    /// Print a trace of the external potential
    void print(FILE* out = outfile) const;

    /// Python print helper
    void py_print() const { print(outfile); }

};

}

#endif
