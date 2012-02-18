#ifndef SOLVER_H 
#define SOLVER_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Matrix;
class Vector;
class Hamiltonian;
class RHamiltonian;
class UHamiltonian;

class Solver {

// => BASE CLASSES <= //

protected:
    
    /// Print flag, defaults to 1
    int print_;	
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Name of solver (set by subclasses)
    std::string name_;
    /// Memory available, in doubles, defaults to 0 => Unlimited storage 
    unsigned long int memory_;	

    /// Convergence criteria, defaults to 1.0E-6
    double criteria_;
    /// Maximum number of iterations, defaults to 100
    int maxiter_;
    /// Converged or not?
    bool converged_;
    /// Convergence measure at this iteration 
    double convergence_;
    /// Current iteration count
    int iteration_;
    /// Preconditioner type
    std::string precondition_;

    /// Common initialization
    void common_init();
    
public:
    // => Constructors < = //

    /// Default Constructor
    Solver();
    /// Destructor
    virtual ~Solver();

    // => Knobs <= // 

    /// Set precondition type (specific to solver type) 
    void set_precondition(const std::string& precondition) { precondition_ = precondition; }
    /// Set maximum vector storage space (defaults to 0 MB => Unlimited storage)
    void set_memory(unsigned long int memory) { memory_ = memory; }
    /// Set maximum number of iterations (defaults to 100)
    void set_maxiter(int maxiter) { maxiter_ = maxiter; }
    /// Set convergence criteria (defaults to 1.0E-6)
    void set_convergence(double criteria) { criteria_ = criteria; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }

    // => Accessors <= //
    
    /// What iteration is it?
    int iteration() const { return iteration_; }
    /// Did the solver converge?
    bool converged() const { return converged_; }
    /**
    * Print header information regarding Solver
    * type on output file
    */
    virtual void print_header() const = 0;
    /// Estimate of maximum memory usage (in doubles)
    virtual unsigned long int memory_estimate() = 0;

    // => Computers <= //

    /** 
     * Warm it all up, everything you got.
     * Cmon you apes, you wanna live forever?!
     */
    virtual void initialize() = 0;
    /**
     * Perform iterative solution,
     * based on current state
     */
    virtual void solve() = 0;
    /**
     * Method to clear off overhead memory 
     * without destroying the object
     */
    virtual void finalize() = 0;
    

};

class RSolver : public Solver {

protected:
    /// Reference to underlying RHamiltonian
    boost::shared_ptr<RHamiltonian> H_;

public:
    // => Constructors < = //

    RSolver(boost::shared_ptr<RHamiltonian> H);
    /// Destructor
    virtual ~RSolver();   

    // => Accessors <= //
    
    /**
    * Knob to swap out a Hamiltonian object
    * @param H new RHamiltonian object
    */
    void set_Hamiltonian(boost::shared_ptr<RHamiltonian> H) { H_ = H; } 

    /**
    * Pointer to the Hamiltonian object
    * @return current RHamiltonian object
    */
    boost::shared_ptr<RHamiltonian> H() const { return H_; }
};

class USolver : public Solver {

protected:
    /// Reference to underlying UHamiltonian
    boost::shared_ptr<UHamiltonian> H_;

public:
    // => Constructors < = //

    /// Reference to underlying RHamiltonian
    USolver(boost::shared_ptr<UHamiltonian> H);
    /// Destructor
    virtual ~USolver();   

    // => Accessors <= //
    
    /**
    * Knob to swap out a Hamiltonian object
    * @param H new UHamiltonian object
    */
    void set_Hamiltonian(boost::shared_ptr<UHamiltonian> H) { H_ = H; } 

    /**
    * Pointer to the Hamiltonian object
    * @return current UHamiltonian object
    */
    boost::shared_ptr<UHamiltonian> H() const { return H_; }
};

// => APPLIED CLASSES <= //

class CGRSolver : public RSolver {

protected:

    /// Force vectors
    std::vector<boost::shared_ptr<Vector> > b_;    
    /// Solution vectors
    std::vector<boost::shared_ptr<Vector> > x_;    
    /// Product vectors
    std::vector<boost::shared_ptr<Vector> > Ap_;    
    /// M^-1 x
    std::vector<boost::shared_ptr<Vector> > z_;    
    /// Residual vectors
    std::vector<boost::shared_ptr<Vector> > r_;    
    /// Conjugate directions
    std::vector<boost::shared_ptr<Vector> > p_;    
    /// Alpha values
    std::vector<double> alpha_;
    /// Beta values
    std::vector<double> beta_;
    /// Residual norm
    std::vector<double> r_nrm2_;
    /// z'r 
    std::vector<double> z_r_;
    /// Which vectors are converged?
    std::vector<bool> r_converged_; 
    /// Number of converged vectors
    int nconverged_;
    
    /// Diagonal M, for guess and Jacobi preconditioning
    boost::shared_ptr<Vector> diag_;
    /// A subspace matrix, for preconditioning 
    SharedMatrix A_;
    /// A subspace indices
    std::vector<std::vector<int> > A_inds_;
    /// Shifts (to solve (A-mI)
    std::vector<std::vector<double> > shifts_;     
    /// Number of guess vectors to use for subspace preconditioner
    int nguess_;    

    void setup();
    void guess();
    void residual();
    void products_x();
    void products_p();
    void alpha();
    void update_x();
    void update_r();
    void check_convergence();
    void update_z();
    void beta();
    void update_p();

public:

    CGRSolver(boost::shared_ptr<RHamiltonian> H);
    virtual ~CGRSolver();
    
    /// Static constructor, uses Options object
    static boost::shared_ptr<CGRSolver> build_solver(Options& options,
        boost::shared_ptr<RHamiltonian> H);

    std::vector<boost::shared_ptr<Vector> >& x() { return x_; }  
    std::vector<boost::shared_ptr<Vector> >& b() { return b_; }  

    void print_header() const;
    unsigned long int memory_estimate();
    void initialize();
    void solve();
    void finalize();

    void set_shifts(const std::vector<std::vector<double> >& shifts) { shifts_ = shifts; }
    void set_A(SharedMatrix A, const std::vector<std::vector<int> > inds) { A_ = A; A_inds_ = inds; }
    void set_nguess(int nguess) { nguess_ = nguess; }
};

class DLRSolver : public RSolver {

protected:
    
    // => Control parameters <= //

    /// Number of desired roots
    int nroot_;
    /// Required norm for subspace expansion
    double norm_;
    /// Maximum allowed subspace size 
    int max_subspace_; 
    /// Minimum allowed subspace size (after collapse)
    int min_subspace_; 
    /// Number of guess vectors to build
    int nguess_;    

    // => Iteration values <= //

    /// The current subspace size
    int nsubspace_;
    /// The number of converged roots
    int nconverged_;

    // => State values <= //
 
    /// Current eigenvectors (nroots)
    std::vector<boost::shared_ptr<Vector> > c_;
    /// Current eigenvalues (nroots)
    std::vector<std::vector<double> > E_;     
    /// B vectors (nsubspace)
    std::vector<boost::shared_ptr<Vector> > b_;
    /// Sigma vectors (nsubspace)
    std::vector<boost::shared_ptr<Vector> > s_;
    /// Delta Subspace Hamiltonian (preconditioner)
    SharedMatrix A_;
    /// Delta Subspace indices
    std::vector<std::vector<int> > A_inds_;
    /// G_ij Subspace Hamiltonian (nsubspace x nsubspace)
    SharedMatrix G_;
    /// Subspace eigenvectors (nsubspace x nsubspace)
    SharedMatrix a_;
    /// Subspace eigenvalues (nsubspace)   
    boost::shared_ptr<Vector> l_; 
    /// Residual vectors (nroots)
    std::vector<boost::shared_ptr<Vector> > r_;
    /// Residual vector 2-norms (nroots)
    std::vector<double> n_;
    /// Correction vectors (nroots)
    std::vector<boost::shared_ptr<Vector> > d_;
    /// Diagonal of Hamiltonian 
    boost::shared_ptr<Vector> diag_;

    // => Run routines <= // 
 
    // Guess, based on diagonal
    void guess();
    // Compute sigma vectors for the given set of b
    void sigma();
    // Compute subspace Hamiltonian
    void subspaceHamiltonian();
    // Diagonalize subspace Hamiltonian
    void subspaceDiagonalization();
    // Find eigenvectors 
    void eigenvecs();
    // Find eigenvalues
    void eigenvals();   
    // Find residuals, update convergence 
    void residuals();
    // Find correctors 
    virtual void correctors();
    // Orthogonalize/add significant correctors 
    void subspaceExpansion();
    // Collapse subspace if needed
    void subspaceCollapse();

public:

    // => Constructors <= //

    /// Constructor
    DLRSolver(boost::shared_ptr<RHamiltonian> H);  
    /// Destructor
    virtual ~DLRSolver();

    /// Static constructor, uses Options object
    static boost::shared_ptr<DLRSolver> build_solver(Options& options,
        boost::shared_ptr<RHamiltonian> H);

    // => Required Methods <= //
    
    virtual void print_header() const;
    virtual unsigned long int memory_estimate();
    virtual void initialize();
    void solve();
    void finalize();
    
    // => Accessors <= //
    
    /// Eigenvectors, by state
    const std::vector<boost::shared_ptr<Vector> >& eigenvectors() const { return c_; }
    /// Eigenvalues, by state/irrep
    const std::vector<std::vector<double> >& eigenvalues() const { return E_; }

    // => Knobs <= //

    /// Set number of roots (defaults to 1)
    void set_nroot(int nroot) { nroot_ = nroot; }
    /// Set maximum subspace size (defaults to 6)
    void set_max_subspace(double max_subspace) { max_subspace_ = max_subspace; }
    /// Set minimum subspace size, for collapse (defaults to 2)
    void set_min_subspace(double min_subspace) { min_subspace_ = min_subspace; }
    /// Set number of guesses (defaults to 1)
    void set_nguess(int nguess) { nguess_ = nguess; }
    /// Set norm critera for adding vectors to subspace (defaults to 1.0E-6) 
    void set_norm(double norm) { norm_ = norm; }
};

class RayleighRSolver : public DLRSolver {

protected:

    /// Turn an eigenproblem into a linear equations problem 
    boost::shared_ptr<CGRSolver> cg_;

    std::string precondition_steps_;
    int precondition_maxiter_;
    std::string quantity_;
    
    // Find correctors (Like, the most advanced correctors ever)
    void correctors();

public:

    // => Constructors <= //

    /// Constructor
    RayleighRSolver(boost::shared_ptr<RHamiltonian> H);  
    /// Destructor
    virtual ~RayleighRSolver();

    /// Static constructor, uses Options object
    static boost::shared_ptr<RayleighRSolver> build_solver(Options& options,
        boost::shared_ptr<RHamiltonian> H);

    // => Required Methods <= //
    
    void print_header() const;
    void initialize();
    void finalize();

    void set_precondition_maxiter(int maxiter) { precondition_maxiter_ = maxiter; }
    void set_precondition_steps(const std::string& steps) { precondition_steps_ = steps; }
    void set_quantity(const std::string& quantity) { quantity_ = quantity; }
};

class DLRXSolver : public RSolver {

protected:
    
    // => Control parameters <= //

    /// Number of desired roots
    int nroot_;
    /// Required norm for subspace expansion
    double norm_;
    /// Maximum allowed subspace size 
    int max_subspace_; 
    /// Minimum allowed subspace size (after collapse)
    int min_subspace_; 
    /// Number of guess vectors to build
    int nguess_;    

    // => Iteration values <= //

    /// The current subspace size
    int nsubspace_;
    /// The number of converged roots
    int nconverged_;

    // => State values <= //
 
    /// Current eigenvectors (nroots)
    std::vector<boost::shared_ptr<Vector> > c_;
    /// Current eigenvalues (nroots)
    std::vector<std::vector<double> > E_;     
    /// B vectors (nsubspace)
    std::vector<boost::shared_ptr<Vector> > b_;
    /// Sigma vectors (nsubspace)
    std::vector<boost::shared_ptr<Vector> > s_;
    /// G_ij Subspace Hamiltonian (nsubspace x nsubspace)
    SharedMatrix G_;
    /// Subspace eigenvectors (nsubspace x nsubspace)
    SharedMatrix a_;
    /// Subspace eigenvalues (nsubspace)   
    boost::shared_ptr<Vector> l_; 
    /// Residual vectors (nroots)
    std::vector<boost::shared_ptr<Vector> > r_;
    /// Residual vector 2-norms (nroots)
    std::vector<double> n_;
    /// Correction vectors (nroots)
    std::vector<boost::shared_ptr<Vector> > d_;
    /// Diagonal of Hamiltonian 
    boost::shared_ptr<Vector> diag_;

    // => Run routines <= // 
 
    // Guess, based on diagonal
    void guess();
    // Compute sigma vectors
    void sigma();
    // Compute subspace Hamiltonian
    void subspaceHamiltonian();
    // Diagonalize subspace Hamiltonian
    void subspaceDiagonalization();
    // Find eigenvectors 
    void eigenvecs();
    // Find eigenvalues
    void eigenvals();   
    // Find residuals, update convergence 
    void residuals();
    // Find correctors 
    void correctors();
    // Orthogonalize/add significant correctors 
    void subspaceExpansion();
    // Collapse subspace if needed
    void subspaceCollapse();

public:

    // => Constructors <= //

    /// Constructor
    DLRXSolver(boost::shared_ptr<RHamiltonian> H);  
    /// Destructor
    virtual ~DLRXSolver();

    /// Static constructor, uses Options object
    static boost::shared_ptr<DLRXSolver> build_solver(Options& options,
        boost::shared_ptr<RHamiltonian> H);

    // => Required Methods <= //
    
    void print_header() const;
    unsigned long int memory_estimate();
    void initialize();
    void solve();
    void finalize();
    
    // => Accessors <= //
    
    /// Eigenvectors, by state
    const std::vector<boost::shared_ptr<Vector> >& eigenvectors() const { return c_; }
    /// Eigenvalues, by state/irrep
    const std::vector<std::vector<double> >& eigenvalues() const { return E_; }

    // => Knobs <= //

    /// Set number of roots (defaults to 1)
    void set_nroot(int nroot) { nroot_ = nroot; }
    /// Set maximum subspace size (defaults to 6)
    void set_max_subspace(double max_subspace) { max_subspace_ = max_subspace; }
    /// Set minimum subspace size, for collapse (defaults to 2)
    void set_min_subspace(double min_subspace) { min_subspace_ = min_subspace; }
    /// Set number of guesses (defaults to 1)
    void set_nguess(int nguess) { nguess_ = nguess; }
    /// Set norm critera for adding vectors to subspace (defaults to 1.0E-6) 
    void set_norm(double norm) { norm_ = norm; }
};

}
#endif
