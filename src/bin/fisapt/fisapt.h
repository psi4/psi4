#ifndef FISAPT_H
#define FISAPT_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>
#include <boost/tuple/tuple.hpp>

namespace psi {

class JK;

namespace fisapt {

class FISAPT { 

protected:

    /// Global options object 
    Options& options_;
    /// Memory in doubles
    size_t doubles_;
    /// Reference wavefunction 
    boost::shared_ptr<Wavefunction> reference_;

    /// Orbital Basis Set (full molecule) 
    boost::shared_ptr<BasisSet> primary_;

    /// Global JK object
    boost::shared_ptr<JK> jk_;

    /// Map of scalars 
    std::map<std::string, double> scalars_;
    /// Map of matrices
    std::map<std::string, boost::shared_ptr<Vector> > vectors_;
    /// Map of matrices
    std::map<std::string, boost::shared_ptr<Matrix> > matrices_;

    // => FISAPT 0th-Order Wavefunction <= //

    /// Common initialization (bases, orbitals, eigenvalues, etc)
    void common_init();
    /// Print header, bases, sizes, etc
    void print_header();
    /// Localize the active occupied orbitals via IBO2
    void localize();
    /// Partition the nuclei and electrons
    void partition();
    /// Build the overlap integrals S
    void overlap();
    /// Build the kinetic integrals T
    void kinetic();
    /// Build the nuclear potentials V and interaction energies
    void nuclear();
    /// Build the J/K potentials for C, D, and E
    void coulomb();
    /// Solve the relaxed SCF equations for A0 and B0
    void scf();
    /// Freeze the core orbitals 
    void freeze_core();
    /// Produce unified matrices for A', B', and C'
    void unify();
    /// Plot some analysis files
    void plot();

    // => F-SAPT0 <= //

    /// Localize 
    void flocalize();
    /// Electrostatics
    void felst();
    /// Exchange
    void fexch();
    /// Induction
    void find();
    /// Dispersion
    void fdisp();
    /// Output
    void fdrop();

    // => SAPT0 <= //

    /// Delta HF
    void dHF();
    /// Electrostatics
    void elst();
    /// Exchange
    void exch();
    /// Induction
    void ind();
    /// Dispersion
    void disp();
    /// Print SAPT results
    void print_trailer();
     
    // Build the ExchInd20 potential in the monomer A ov space
    boost::shared_ptr<Matrix> build_exch_ind_pot(std::map<std::string, boost::shared_ptr<Matrix> >& vars);
    // Build the Ind20 potential in the monomer A ov space
    boost::shared_ptr<Matrix> build_ind_pot(std::map<std::string, boost::shared_ptr<Matrix> >& vars);

    /// Helper to drop a matrix to filepath/A->name().dat
    void drop(boost::shared_ptr<Matrix> A, const std::string& filepath);
    /// Helper to drop a vector to filepath/A->name().dat
    void drop(boost::shared_ptr<Vector> A, const std::string& filepath);
    /// Helper to extract columns from a matrix
    static boost::shared_ptr<Matrix> extract_columns(
        const std::vector<int>& cols,
        boost::shared_ptr<Matrix> A);

public:
    /// Initialize an FISAPT object with an SCF reference
    FISAPT(boost::shared_ptr<Wavefunction> scf);  
    virtual ~FISAPT();

    /// Gogo!
    void compute_energy(); 

};

class FISAPTSCF {

protected:


    /// Global options object 
    Options& options_;
    
    /// Global JK object
    boost::shared_ptr<JK> jk_;

    /// Map of scalars 
    std::map<std::string, double> scalars_;
    /// Map of matrices
    std::map<std::string, boost::shared_ptr<Vector> > vectors_;
    /// Map of matrices
    std::map<std::string, boost::shared_ptr<Matrix> > matrices_;

    /// Print orbitals
    void print_orbitals(
        const std::string& header,
        int start,
        boost::shared_ptr<Vector> eps
        );
        
   
public:

    FISAPTSCF(
        boost::shared_ptr<JK> jk,    // JK object
        double enuc,                 // Nuclear repulsion energy
        boost::shared_ptr<Matrix> S, // Overlap integrals
        boost::shared_ptr<Matrix> X, // Restricted orthogonalization matrix [nbf x nmo]
        boost::shared_ptr<Matrix> T, // Kinetic integrals
        boost::shared_ptr<Matrix> V, // Potential integrals
        boost::shared_ptr<Matrix> W, // External embedding potential
        boost::shared_ptr<Matrix> C  // Guess for occupied orbitals [nbf x nocc]
        );
    virtual ~FISAPTSCF();

    void compute_energy();
 
    std::map<std::string, double>& scalars()                      { return scalars_; }
    std::map<std::string, boost::shared_ptr<Vector> >& vectors()  { return vectors_; }
    std::map<std::string, boost::shared_ptr<Matrix> >& matrices() { return matrices_; }

};

class CPHF_FISAPT {

friend class FISAPT;

protected:

    // => Global Data <= //

    // Convergence tolerance
    double delta_;
    // Maximum allowed iterations
    int maxiter_;
    // JK Object 
    boost::shared_ptr<JK> jk_;
    
    // => Monomer A Problem <= //

    // Perturbation applied to A
    boost::shared_ptr<Matrix> w_A_;
    // Response of A
    boost::shared_ptr<Matrix> x_A_;
    // Active occ orbital coefficients of A
    boost::shared_ptr<Matrix> Cocc_A_;
    // Active vir orbital coefficients of A
    boost::shared_ptr<Matrix> Cvir_A_;
    // Active occ orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_occ_A_;
    // Active vir orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_vir_A_;

    // => Monomer B Problem <= //

    // Perturbation applied to B
    boost::shared_ptr<Matrix> w_B_;
    // Response of B
    boost::shared_ptr<Matrix> x_B_;
    // Active occ orbital coefficients of B
    boost::shared_ptr<Matrix> Cocc_B_;
    // Active vir orbital coefficients of B
    boost::shared_ptr<Matrix> Cvir_B_;
    // Active occ orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_occ_B_;
    // Active vir orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_vir_B_;

    // Form the s = Ab product for the provided vectors b (may or may not need more iterations)
    std::map<std::string, boost::shared_ptr<Matrix> > product(std::map<std::string, boost::shared_ptr<Matrix> > b);
    // Apply the denominator from r into z
    void preconditioner(boost::shared_ptr<Matrix> r,
                        boost::shared_ptr<Matrix> z,
                        boost::shared_ptr<Vector> o,
                        boost::shared_ptr<Vector> v);

public:
    CPHF_FISAPT();
    virtual ~CPHF_FISAPT();

    void compute_cphf();
};

} // Namespace fisapt 

} // Namespace psi

#endif
