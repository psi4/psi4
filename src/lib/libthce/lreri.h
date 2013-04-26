#ifndef LRERI_H
#define LRERI_H

#include <libmints/typedefs.h>

namespace psi {

class Tensor;
class CoreTensor;
class DiskTensor;

class LRERI {
    
protected:

    // => Utility <= //

    /// Print flag
    int print_;
    /// Debug flag
    int debug_;
    /// Bench flag
    int bench_;
    /// Memory in doubles
    unsigned long int memory_;

    // => Basis Set <= //

    /// Primary orbital basis (nso)
    boost::shared_ptr<BasisSet> primary_;
    
    // => Orbital Spaces <= //
    
    /// Occupation matrix coefficients (nso x nmo)
    boost::shared_ptr<Matrix> C_; 
    /// Orbital spaces, each defined by a keyword and the index range in <start, end+1>. 
    std::map<std::string, std::pair<int, int> > spaces_; 
    /// Orbital spaces order buffer, to keep the printing nice.
    std::vector<std::string> spaces_order_;

    // => Utility Routines <= //
    
    /// Set defaults
    void common_init();
    /// Print info
    virtual void print_header() = 0;

    /// Inverse fitting metric
    boost::shared_ptr<Matrix> Jm12(boost::shared_ptr<BasisSet> auxiliary, double condition);

public:
    // => Constructors <= //

    LRERI(boost::shared_ptr<BasisSet> primary);
    virtual ~LRERI();

    // => Defaults <= //

    /// O: Load the usual orbital spaces for an RHF or UHF wavefunction
    virtual void load_wavefunction(boost::shared_ptr<Wavefunction> ref);
    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => Orbital Space Control <= //

    /// R: Set the overall C matrix (calls clear before starting)
    void set_C(boost::shared_ptr<Matrix> C);
    /// R: Add an orbital subspace to the queue. The subspace will be referred to by key, and ranges from start to end-1.  
    void add_space(const std::string& key, int start, int end); 
    /// Clear the C matrix and orbital spaces list
    virtual void clear();

    // => Computers <= // 

    /// R: Compute the desired ERI factorization
    virtual void compute() = 0; 
    
    // => Setters <= //

    /// Set the print flag
    void set_print(int print) { print_ = print; }
    /// Set the debug flag
    void set_debug(int debug) { debug_ = debug; }
    /// Set the bench flag
    void set_bench(int bench) { bench_ = bench; }
    /// Set the allowed memory in doubles
    void set_memory(unsigned long int memory) { memory_ = memory; }

};

/**
 * DFERI
 **/
class DFERI : public LRERI {

protected:

    // => DF Tech <= //

    /// Auxiliary orbital-pair basis
    boost::shared_ptr<BasisSet> auxiliary_;
    /// Relative condition number in J^-1/2
    double J_cutoff_;
    /// Schwarz sieve tolerance
    double schwarz_cutoff_;

    // => Targets <= //

    /// Three-center integrals, by name, sorted e.g. (ov|Q), DiskTensor
    std::map<std::string, boost::shared_ptr<Tensor> > ints_;
    /// Requested pair spaces
    std::map<std::string, std::pair<std::string, std::string> > pair_spaces_;
    /// Order of pair spaces, to keep printing nice
    std::vector<std::string> pair_spaces_order_;

    // => Utility Routines <= //
    
    /// Set defaults
    void common_init();
    /// Print info
    virtual void print_header();

    void allocate();
    void transform();
    void fit();

public:
    // => Constructors <= //

    DFERI(boost::shared_ptr<BasisSet> primary,
          boost::shared_ptr<BasisSet> auxiliary);
    virtual ~DFERI();

    static boost::shared_ptr<DFERI> build(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary, Options& options);
    static boost::shared_ptr<DFERI> build(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary, Options& options, boost::shared_ptr<Wavefunction> ref);

    // => Defaults <= //
    
    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => Pair Space Control <= //

    // R: add an orbital pair space to the list of tasks
    void add_pair_space(const std::string& name, const std::string& space1, const std::string& space2);
    // Clear this DFERI object of tasks
    virtual void clear();

    // => Computers <= // 

    /// R: Compute the requested DF 3-index integrals
    virtual void compute();
    /// Handle to computed disk tensors, by name in add_pair above
    std::map<std::string, boost::shared_ptr<Tensor> >& ints() { return ints_; }

    // => Setters <= //

    /// Set the relative eigenvalue cutoff in J^{-1/2}
    void set_J_cutoff(double J_cutoff) { J_cutoff_ = J_cutoff; }
    /// Set the schwarz sieve cutoff
    void set_schwarz_cutoff(double schwarz_cutoff) { schwarz_cutoff_ = schwarz_cutoff; }
    
};

/**
 * LSTHCERI
 **/
class LSTHCERI : public LRERI {


protected:

    // => LS-LSTHC Tech <= //

    /// Raw AO-basis X matrix (nso x nP, for now)
    boost::shared_ptr<Matrix> X_;
    /// Auxiliary orbital-pair basis
    boost::shared_ptr<BasisSet> auxiliary_;
    /// Relative condition number in J^-1/2
    double J_cutoff_;
    /// Relative condition number in S^-1
    double S_cutoff_;
    /// Schwarz sieve tolerance
    double schwarz_cutoff_;
    /// Balance the X matrices?
    bool balance_;

    // => Targets <= //

    /// ERI factors, CoreTensor, swapped out
    std::map<std::string, std::vector<boost::shared_ptr<Tensor> > > ints_;
    /// METH factors, CoreTensor, swapped out
    std::map<std::string, std::vector<boost::shared_ptr<Tensor> > > meths_;
    /// Requested ERI spaces
    std::map<std::string, std::vector<std::string> > eri_spaces_;
    /// Order of ERI spaces, to keep printing nice
    std::vector<std::string> eri_spaces_order_;

    // => Utility Routines <= //
    
    /// Set defaults
    void common_init();
    /// Print info
    virtual void print_header();

    /// Build all requred X matrices (np x nP, core)
    std::map<std::string, boost::shared_ptr<Tensor> > build_X(bool meth = false);
    /// Build all required E matrices (nA x nP, disk) [deleted]
    std::map<std::string, boost::shared_ptr<Tensor> > build_E(std::map<std::string, boost::shared_ptr<Tensor> >& Xs); 
    /// Build all required inverse S matrices (nP x nP, core, swapped)
    std::map<std::string, boost::shared_ptr<Tensor> > build_S(std::map<std::string, boost::shared_ptr<Tensor> >& Xs, bool meth = false);
    /// Build all requred L matrices (nP x nA, core, swapped)
    std::map<std::string, boost::shared_ptr<Tensor> > build_L(std::map<std::string, boost::shared_ptr<Tensor> >& Es, 
                                                              std::map<std::string, boost::shared_ptr<Tensor> >& Ss); 
    /// Build all required Z matrices (nP x nP, core, swapped)
    std::map<std::string, boost::shared_ptr<Tensor> > build_Z(std::map<std::string, boost::shared_ptr<Tensor> >& Ls); 
    /// Pack up the integrals 
    void pack(std::map<std::string, boost::shared_ptr<Tensor> >& Xs,
              std::map<std::string, boost::shared_ptr<Tensor> >& Zs,
              std::map<std::string, boost::shared_ptr<Tensor> >& Ls,
              std::map<std::string, boost::shared_ptr<Tensor> >& Ss);
    /// Pack up the meth intermediates
    void pack_meth(std::map<std::string, boost::shared_ptr<Tensor> >& Xs,
                   std::map<std::string, boost::shared_ptr<Tensor> >& Ss);

public:
    // => Constructors <= //

    LSTHCERI(boost::shared_ptr<BasisSet> primary,
           boost::shared_ptr<BasisSet> auxiliary,
           boost::shared_ptr<Matrix> X);
    virtual ~LSTHCERI();

    static boost::shared_ptr<LSTHCERI> build(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> X, Options& options);
    static boost::shared_ptr<LSTHCERI> build(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> X, Options& options, boost::shared_ptr<Wavefunction> ref);

    // => Defaults <= //
    
    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => ERI Space Control <= //

    // R: add an eri space to the list of tasks
    void add_eri_space(const std::string& name, const std::string& space1, const std::string& space2, const std::string& space3, const std::string& space4);
    // Clear this LSTHCERI object of tasks
    virtual void clear();

    // => Computers <= // 

    /// R: Compute the requested LS-LSTHC factors
    virtual void compute();
    /// LS-LSTHC factors [X1,X2,Z,X3,X4,L12,L34,Sinv12,Sinv34]
    std::map<std::string, std::vector<boost::shared_ptr<Tensor> > >& ints() { return ints_; };

    /// O: Compute the METH X and Sinv matrices
    virtual void compute_meth();
    /// METH helper factors [X1,X2,Sinv]
    std::map<std::string, std::vector<boost::shared_ptr<Tensor> > >& meths() { return meths_; }

    // => Setters <= //

    /// Set the relative eigenvalue cutoff in J^{-1/2}
    void set_J_cutoff(double J_cutoff) { J_cutoff_ = J_cutoff; }
    /// Set the relative eigenvalue cutoff in S^{-1}
    void set_S_cutoff(double S_cutoff) { S_cutoff_ = S_cutoff; }
    /// Set the schwarz sieve cutoff
    void set_schwarz_cutoff(double schwarz_cutoff) { schwarz_cutoff_ = schwarz_cutoff; }
    /// Set to balance or not?
    void set_balance(bool balance) { balance_ = balance; } 
};


} // End namespace

#endif

