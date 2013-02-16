#ifndef INFSAPT_H
#define INFSAPT_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>

namespace psi {

namespace dftsapt {

/**
* Class INFSAPT provides the trial implementation of the new
* infinite-order explictly-orthogonalized SAPT0 model
* 
* This class is designed to handle N-body decompositions
* for RHF-reference monomers
**/ 
class INFSAPT {

protected:

    // => Parameters <= //

    // Print flag
    int print_;
    // Debug flag
    int debug_;
    // Bench flag
    int bench_;
    
    // Memory in doubles
    unsigned long int memory_;

    // Schwarz cutoff
    double schwarz_;

    // => State Data <= //

    // Energies table
    std::map<std::string, double> energies_;

    // Dimer primary basis set
    boost::shared_ptr<BasisSet> primary_;
    // Dimer -RI or -MP2FIT auxiliary basis set 
    boost::shared_ptr<BasisSet> mp2fit_;

    // Cluster Hartree-Fock
    boost::shared_ptr<Wavefunction> cluster_;
    // Monomer Hartree-Focks
    std::vector<boost::shared_ptr<Wavefunction> > monomers_;

    // => SCF => PT2 Crossover data <= //

    // Orthogonalized monomer occupied orbitals
    std::vector<boost::shared_ptr<Matrix> > Caocc_;
    // Orthogonalized monomer occupied eigenvalues
    std::vector<boost::shared_ptr<Vector> > eps_aocc_;  

    // Joint orthonormal virtual orbitals
    boost::shared_ptr<Matrix> Cavir_; 
    // Joint orthonormal virtual eigenvalues
    boost::shared_ptr<Vector> eps_avir_; 

    // Print author/sizing/spec info
    virtual void print_header() const;
    // Obligatory
    virtual void print_trailer();

    // Handle the SCF-like terms
    virtual void scf_terms();
    // Handle the PT2-like terms
    virtual void pt2_terms(); 

    // => Helper Methods for PT2 <= //

    // BUild the (ia|Q) disk tensor
    void build_iaQ(
        boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<BasisSet> auxiliary,
        boost::shared_ptr<Matrix> Caocc, 
        boost::shared_ptr<Matrix> Cavir);

    // Protected constructor (use factory below)
    INFSAPT();
    // Common initialization
    void common_init();
public:
    // Destructor, frees memory
    virtual ~INFSAPT();

    // Factory constructor, call this with a converged cluster RHF and a vector of monomer RHFs
    static boost::shared_ptr<INFSAPT> build(boost::shared_ptr<Wavefunction> d,
                                            std::vector<boost::shared_ptr<Wavefunction> > m);

    // Compute the INF-SAPT analysis
    virtual double compute_energy();


};

}} // End namespace

#endif

