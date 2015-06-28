#ifndef PROP_H
#define PROP_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>

namespace psi {

class CubicScalarGrid;

class Properties {

protected:
    
    // => Task specification <= //

    /// Global options object
    Options& options_;

    // => Key Member Data <= //

    /// Orbital Basis Set 
    boost::shared_ptr<BasisSet> basisset_;
    /// AO-basis (C1) OPDM for alpha electrons
    boost::shared_ptr<Matrix> Da_;
    /// AO-basis (C1) OPDM for beta electrons
    boost::shared_ptr<Matrix> Db_;
    /// AO-basis (C1) SCF orbital coefficients for alpha electrons
    boost::shared_ptr<Matrix> Ca_;
    /// AO-basis (C1) SCF orbital coefficients for beta electrons
    boost::shared_ptr<Matrix> Cb_;
     
    // => Computers <= //

    /// Grid-based property computer
    boost::shared_ptr<CubicScalarGrid> grid_;

    // => Helper Functions <= //

    /// Common setup across all constructors
    void common_init();

public:
    // => Constructors <= //
    
    /// Construct a Properties object from a Wavefunction (possibly with symmetry in wfn)
    Properties(boost::shared_ptr<Wavefunction> wfn);
    /// Common Destructor
    virtual ~Properties();

    // => High-Level Property Computers <= //

    /// Compute all relevant properties from options object specifications
    void compute_properties();    

    // => Low-Level Property Computers (Do not use unless you are an advanced client code) <= //
    
    /// Obligatory title info 
    void print_header();
    /// Compute a density grid task (key.cube)
    void compute_density(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ESP grid task (Dt.cube and ESP.cube)
    void compute_esp(boost::shared_ptr<Matrix> Dt);
    /// Compute an orbital task (key_N.cube, for 0-based indices of C)
    void compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::string& key);
    /// Compute a basis function task (key_N.cube, for 0-based indices of basisset_)
    void compute_basis_functions(const std::vector<int>& indices, const std::string& key);
    /// Compute a LOL grid task (key.cube)
    void compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ELF grid task (key.cube)
    void compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key);
};

}

#endif
