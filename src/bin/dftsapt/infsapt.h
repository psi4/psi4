#ifndef INFSAPT_H
#define INFSAPT_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>

namespace psi {

class JK;

namespace dftsapt {

class DIISN;

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

    // => Startup Data <= //
    
    // Cluster Hartree-Fock
    boost::shared_ptr<Wavefunction> cluster_;
    // Monomer Hartree-Focks
    std::vector<boost::shared_ptr<Wavefunction> > monomers_;

    // => State Data <= //

    // Energies table
    std::map<std::string, double> energies_;

    // Dimer primary basis set
    boost::shared_ptr<BasisSet> primary_;
    // Dimer -RI or -MP2FIT auxiliary basis set 
    boost::shared_ptr<BasisSet> mp2fit_;

    // => SCF => PT2 Crossover data <= //

    // Orthogonalized monomer occupied orbitals
    std::vector<boost::shared_ptr<Matrix> > Caocc_;
    // Orthogonalized monomer occupied eigenvalues
    std::vector<boost::shared_ptr<Vector> > eps_aocc_;  

    // Joint orthonormal virtual orbitals
    std::vector<boost::shared_ptr<Matrix> > Cavir_; 
    // Joint orthonormal virtual eigenvalues
    std::vector<boost::shared_ptr<Vector> > eps_avir_; 
    
    // => Main Runtime Functions <= //

    // Print sizing/spec info
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
        std::vector<boost::shared_ptr<Matrix> >& Caocc, 
        std::vector<boost::shared_ptr<Matrix> >& Cavir);

    // Protected constructor (use factory below)
    INFSAPT();
    void common_init();
public:
    virtual ~INFSAPT();

    // Factory constructor, call this with a converged cluster RHF and a vector of monomer RHFs
    static boost::shared_ptr<INFSAPT> build(boost::shared_ptr<Wavefunction> d,
                                            std::vector<boost::shared_ptr<Wavefunction> > m);

    // Compute the INF-SAPT analysis
    virtual double compute_energy();


};

/**
 * Class PB provides a utility for performing
 * Pauli Blockade RHF in a C1 dimer-centered basis 
 **/
class PB {

friend class INFSAPT;

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

    // => SCF Parameters <= //

    // Penalty parameter
    double lambda_;

    // Canonical orthogonalization cutoff
    double S_cutoff_;

    // Maximum number of PB iterations to attempt
    int maxiter_;
    // Convergence criteria (energy)
    double E_convergence_;
    // Convergence criteria (denominator)
    double D_convergence_;

    // Do DIIS?
    bool diis_;
    // Iteration to start saving Fock matrices
    int diis_start_;
    // Min number of vectors to extrapolate with
    int diis_min_;
    // Max number of vectors to extrapolate with
    int diis_max_;
    // DIIS manager
    boost::shared_ptr<DIISN> diis_manager_;

    // => Startup Data <= //

    // Monomer Hartree-Focks
    std::vector<boost::shared_ptr<Wavefunction> > monomers_;
    
    // => Skeleton Data <= //

    // Cluster overlap matrix
    boost::shared_ptr<Matrix> S_;
    // Cluster orthoognalization matrix
    boost::shared_ptr<Matrix> X_;
    // Cluster kinetic matrix
    boost::shared_ptr<Matrix> T_;
    // Monomer potential matrices 
    std::vector<boost::shared_ptr<Matrix> > V_;
    // JK utility object (holds monomer J,K,D)
    boost::shared_ptr<JK> jk_;
    // Fock matrices (with penalty)
    std::vector<boost::shared_ptr<Matrix> > F_;
    // Commutators
    std::vector<boost::shared_ptr<Matrix> > Q_;

    // => State Data <= //

    // Occupied C matrices (guess from 
    std::vector<boost::shared_ptr<Matrix> > Cocc_;
    // Joint virtual C matrices
    std::vector<boost::shared_ptr<Matrix> > Cvir_;
    // Occupied eigenvalues
    std::vector<boost::shared_ptr<Vector> > eps_occ_;
    // Joint virtual eigenvalues
    std::vector<boost::shared_ptr<Vector> > eps_vir_;

    // Initial RHF energies (recomputed)
    boost::shared_ptr<Vector> E_HF_0_;
    // Last iteration's RHF energies (recomputed)
    boost::shared_ptr<Vector> E_HF_L_;
    // Current RHF energies (recomputed)
    boost::shared_ptr<Vector> E_HF_A_;

    // Current RHF energy (total)
    double E_;
    // Unconstrained RHF energy (total)
    double E_0_;

    // Original <a|b> overlaps, for reference
    boost::shared_ptr<Matrix> S_0_;

    // Vector with starting position of each monomer's occupied space (last element is total nocc)
    std::vector<int> occ_starts_;

    // => Convergence Data <= //
    
    // Sum of E_HF_A - E_HF_L
    double dE_;
    // Max RMS of commutator
    double dQ_;

    // => Helper Methods <= //

    // Compute the symmetric orthogonalization guess, placing the orthonormalized orbitals in Cocc_
    void symmetric_orthogonalize();
    // Build the Fock operators from T/V/J/K and the projectors
    void build_fock();
    // Build the commutators
    void build_commutators();
    // Check to see if we are converged
    bool check_convergence();
    // Perform DIIS extrapolation
    bool diis();
    // Diagonalize the Fock matrices
    void diagonalize();
    // Recompute <a|F_A|a> without the penalty (orbitals do not change, at present)
    void purify_eigenvalues(); 
    // Print a trace of the energies/eigenvalues for each monomer
    void print_summary();
    // Analyze the remaining orthogonality error
    void print_errors();
    
    // Compute the current value of <a|b> for all monomers
    boost::shared_ptr<Matrix> compute_Soo(); 
    // Compute the current HF energy for all monomers
    boost::shared_ptr<Vector> HF_energy();

    // Protected constructor (use factory below)
    PB();
    void common_init();
public:
    virtual ~PB();    

    // Factory constructor, call this with a vector of converged monomer RHFs (non-orthogonal)
    static boost::shared_ptr<PB> build(std::vector<boost::shared_ptr<Wavefunction> > m);
    
    // => Accessors <= //

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
    void set_bench(int bench) { bench_ = bench; }
    void set_memory(unsigned long int memory) { memory_ = memory; }

    boost::shared_ptr<JK> jk() const { return jk_; }
    const std::vector<boost::shared_ptr<Matrix> >& V() const { return V_; }

    boost::shared_ptr<Vector> E_HF_0() const { return E_HF_0_; }
    boost::shared_ptr<Vector> E_HF_A() const { return E_HF_A_; }

    const std::vector<boost::shared_ptr<Matrix> >& Cocc() const { return Cocc_; }
    const std::vector<boost::shared_ptr<Matrix> >& Cvir() const { return Cvir_; }
    const std::vector<boost::shared_ptr<Vector> >& eps_occ() const { return eps_occ_; }
    const std::vector<boost::shared_ptr<Vector> >& eps_vir() const { return eps_vir_; }
    
    // Print sizing/spec info
    virtual void print_header() const;
    
    // => Computers <= // 

    // Initialize the JK object and S/X/T/V/D/J/K/F/Q matrices for the non-orthogonal orbitals
    void initialize();

    // Compute the PB SCF
    void compute();

};

/**
 * Class DIISN provides DIIS capabilities with 
 * arbitrary number of state (x) and residual (f)
 * quantities
 *
 * x and f are both based on Matrix objects, which,
 * at least for now, must all be C1 and the same
 * size from quantity to quantity
 **/
class DIISN {

protected:

    // Name of this DIIS object, for file mirror
    std::string name_;
    // Minimum number of vectors to extrapolate with
    int min_vecs_;
    // Maximum number of vectors to extrapolate with
    int max_vecs_;

    // Number of bodies
    int N_;
    // Number of elements in each x
    size_t x_numel_;
    // Number of elements in each f
    size_t f_numel_;

    // File pointer to DIIS file
    FILE* fh_;
    
    // Current <f|f> metric
    boost::shared_ptr<Matrix> metric_;
    // Current iteration
    int iter_;

    // Namespaced scratch file name for this DIIS object
    std::string filename() const;

public:
    /**
    * Master DIISN constructor
    *
    * @name - unique name for this DIIS object (used in filename)
    * @min_vecs - min vectors to start extrapolating with
    * @max_vecs - max vectors to extrapolate with (worst residual removal)
    * @x/f - used to determine N, x_numel, f_numel
    *
    * x - state variables (extrapolated quantities)
    * f - residual variables (skeleton quantities)
    **/
    DIISN(const std::string& name,
          int min_vecs,
          int max_vecs,
          const std::vector<boost::shared_ptr<Matrix> >& x,
          const std::vector<boost::shared_ptr<Matrix> >& f);

    // Destructor, closes/deletes DIIS file
    virtual ~DIISN();

    // Add entry and possibly extrapolate into x
    bool diis(std::vector<boost::shared_ptr<Matrix> >& x,
              std::vector<boost::shared_ptr<Matrix> >& f); 
     
    
};

}} // End namespace

#endif

