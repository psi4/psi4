#ifndef _psi_src_lib_libmints_extern_potential_h_
#define _psi_src_lib_libmints_extern_potential_h_

#include <vector>
#include "typedefs.h"
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

    /// Debug flag
    int debug_;
    /// Print flag
    int print_;

    /// Name of potential
    std::string name_;
    /// <Z,x,y,z> array of charges
    std::vector<boost::tuple<double,double,double,double> > charges_;
    /// Auxiliary basis sets (with accompanying molecules and coefs) of diffuse charges
    std::vector<std::pair<boost::shared_ptr<BasisSet>, SharedVector> > bases_;

public:
    /// Constructur, does nothing
    ExternalPotential();
    /// Destructor, does nothing
    ~ExternalPotential();

    /// Set name 
    void setName(const std::string & name) { name_ = name; }    

    /// Add a charge Z at (x,y,z)
    void addCharge(double Z,double x, double y, double z);
    /// Add a basis of S auxiliary functions with DF coefficients
    void addBasis(boost::shared_ptr<BasisSet> basis, SharedVector coefs);

    /// Reset the field to zero (eliminates all entries) 
    void clear();
    
    /// Compute the external potential matrix in the given basis set
    SharedMatrix computePotentialMatrix(boost::shared_ptr<BasisSet> basis);
    /// Compute the contribution to the nuclear repulsion energy for the given molecule
    double computeNuclearEnergy(boost::shared_ptr<Molecule> mol);
    
    /// Print a trace of the external potential
    void print(FILE* out = outfile) const;

    /// Python print helper
    void py_print() const { print(outfile); }

    /// Print flag 
    void set_print(int print) { print_ = print; } 
    /// Debug flag          
    void set_debug(int debug) { debug_ = debug; } 

};

}

#endif
