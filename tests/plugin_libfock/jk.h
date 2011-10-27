#ifndef JK_H 
#define JK_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class GPUDFJKHelper;
class BasisSet;
class Matrix;
class IntegralFactory;
class ERISieve;
class TwoBodyAOInt;
class Options;
class FittingMetric;
class PSIO;

// => BASE CLASS <= //

/**
 * Class JK 
 *
 * Class to compute Generalized Fock-Matrix Contributions
 * of the form:
 *
 *  J_mn = (mn|ls) C_li^left C_si^right
 *  K_mn = (ml|ns) C_li^left C_si^right
 *  wK_mn = (ml|w|ns) C_li^left C_si^right
 *
 * These matrices appear prominently in Hartree-Fock, CIS,
 * and CPHF theory. This class (and associated subclasses)
 * attempt to provide a uniform framework for quickly 
 * obtaining these matrices using a variety of techniques. 
 *
 * This class is abstract, specific instances must be obtained
 * by constructing an object corresponding to the desired
 * algorithm's subclass of JK, e.g., DFJK or DirectJK.  
 *
 * This class is available for symmetric or nonsymmetric C
 * (this refers to C^left = C^right or not). Symmetric or 
 * nonsymmetric behavior is obtained by using the 
 * JKAlgorithm(C,basis) or JKAlgorithm(C,C,basis) constructors,
 * respectively. Note that nonsymmetric can be up to 2x slower
 * than symmetric for certain algorithms.
 *
 * This class is designed to perform efficiently for any
 * number of C^left/C^right combinations. The constructor
 * takes a std::vector of C^left/C^right, which can be used
 * to provide spin-specialization in SCF, arbitrary numbers
 * of trial vectors in CIS, arbitrary numbers of perturbations
 * in CPHF, etc. The number of C^left/C^right pairs may change from 
 * iteration to iteration, and the second dimension (typically
 * noccpi) may change from pair to pair or iteration to iteration.
 *
 * Tasking for J/K/wK is controlled by JK::set_do_J(), JK::set_do_K(),
 * and JK::set_do_wK() knobs. By default, J and K matrices are built,
 * wK matrices are not built. omega is set using the JK::set_omega()
 * knob. set_omega() may be called between iterations, set_do_X() must
 * be called before init().
 *
 * Spatial symmetry is supported, though most of the algorithms
 * backtransform to C1 for better sieving/scaling properties. The
 * results are available after each call to JK::compute() in the
 * methods JK::J(), JK::K(), JK::wK() and JK::D(), and will have the same
 * size and symmetry structure as the forcing C^left/C^right 
 * vectors. J(), K(), wK(), and D() all return in the USO basis.
 *
 * OpenMP and parallel BLAS/LAPACK threading are targeted
 * where possible.
 *
 * The typical calling convention for a JK instance is:
 * 
 *      // These do not even need to be filled yet, 
 *      // jk just needs the reference
 *      std::vector<boost::shared_ptr<Matrix> > C_left;
 *      std::vector<boost::shared_ptr<Matrix> > C_right;
 *
 *      // Nonymmetric constructor, Algorithm corresponds 
 *      // to Type 
 *      boost::shared_ptr<JKType> jk(new JKType(
 *          C_left, C_right, basis, ...));
 *
 *      // Set any desired knobs
 *      
 *      // 8 GB Memory, 1 G doubles
 *      jk->set_memory(1000000000L);    
 *      // Cutoff of 1.0E-12   
 *      jk->set_cutoff(1.0E-12);
 *      // Do J/K, Not wK (superfluous)
 *      jk->set_do_J(true);
 *      jk->set_do_K(true);
 *      jk->set_do_wK(false);
 *      ...
 *
 *      // Initialize (builds integrals)
 *      jk->initialize();
 *
 *      // Enter iterations or whatever
 *      // In each iteration:
 *      
 *      // Modify your reference to C_left/C_right
 *      // You can do whatever, resize, delete/rebuild,
 *      // or just change the double** of each element.
 *      // All that matters is that you still refer to 
 *      // C_left/C_right
 *      ...
 *
 *      // Let jk compute for the given C_left/C_right
 *      jk->compute();
 *
 *      // In return for the flexibility I give you with 
 *      // C_left/C_right, I ask that you renew your reference
 *      // to the results J/K/D here. I only malloc where needed,
 *      // but if anything at all happens to C_left/C_right, last
 *      // iteration's boost::shared_ptr<Matrix> values might not
 *      // be valid. The std::vector<boost::shard_ptr<Matrix> >
 *      // for J/K/D will still be valid, but the pointers in the
 *      // entries may change.
 *
 *      boost::shared_ptr<Matrix> Jnew = jk->J()[0];
 *      boost::shared_ptr<Matrix> Knew = jk->K()[0];
 *
 *      // After iterations:
 *
 *      // Finalize (frees most memory, including C/D/J/K
 *      // pointers you are no longer holding)
 *      jk->finalize();
 */
class JK {

protected:

    // => Utility Variables <= //

    /// Print flag, defaults to 1
    int print_;     
    /// Debug flag, defaults to 0
    int debug_;     
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    unsigned long int memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, omp_get_max_thread() otherwise)
    int omp_nthread_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;

    // => Tasks <= //
    
    /// Do J matrices? Defaults to true
    bool do_J_;     
    /// Do K matrices? Defaults to true
    bool do_K_;     
    /// Do wK matrices? Defaults to false
    bool do_wK_;    

    /// Omega, defaults to 0.0
    double omega_;

    // => Architecture-Level State Variables (Spatial Symmetry) <= // 

    /// C_left_ = C_right_?
    bool lr_symmetric_;
    /// Pseudo-occupied C matrices, left side
    std::vector<boost::shared_ptr<Matrix> > C_left_;
    /// Pseudo-occupied C matrices, right side 
    std::vector<boost::shared_ptr<Matrix> > C_right_;
    /// Pseudo-density matrices D_ls =  C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > D_;
    /// J matrices: J_mn = (mn|ls) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > J_;
    /// K matrices: K_mn = (ml|ns) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > K_;
    /// wK matrices: wK_mn = (ml|w|ns) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > wK_;

    // => Microarchitecture-Level State Variables (No Spatial Symmetry) <= // 

    /// Primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// AO2USO transformation matrix
    boost::shared_ptr<Matrix> AO2USO_;
    /// Pseudo-occupied C matrices, left side
    std::vector<boost::shared_ptr<Matrix> > C_left_ao_;
    /// Pseudo-occupied C matrices, right side 
    std::vector<boost::shared_ptr<Matrix> > C_right_ao_;
    /// Pseudo-density matrices
    std::vector<boost::shared_ptr<Matrix> > D_ao_;
    /// J matrices: J_mn = (mn|ls) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > J_ao_;
    /// K matrices: K_mn = (ml|ns) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > K_ao_;
    /// wK matrices: wK_mn = (ml|w|ns) C_li^left C_si^right
    std::vector<boost::shared_ptr<Matrix> > wK_ao_;

    // => Per-Iteration Setup/Finalize Routines <= //

    /// Build the pseudo-density D_, before compute_JK()
    void compute_D();
    /// Transform current C_left_/C_right_/D_ to C_left_ao_/C_right_ao_/D_ao_, before compute_JK()
    void USO2AO();
    /// Transform finished J_ao_/K_ao_ to J_/K_, after compute_JK()
    void AO2USO();
    /// Allocate J_/K_ should we be using SOs
    void allocate_JK(); 
    /// Common initialization
    void common_init();

    // => Required Algorithm-Specific Methods <= //

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const = 0;
    /// Setup integrals, files, etc
    virtual void preiterations() = 0; 
    /// Compute J/K for current C/D 
    virtual void compute_JK() = 0;
    /// Delete integrals, files, etc
    virtual void postiterations() = 0; 

    // => Helper Routines <= //

    /// Memory (doubles) used to hold J/K/wK/C/D and ao versions, at current moment
    unsigned long int memory_overhead();

public:
    // => Constructors <= //

    /**
     * Non-Symmetric Constructor
     * @param C_left reference to std::vector that will hold 
     *        left-side C matrices
     * @param C_right reference to std::vector that will hold 
     *        right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    JK(std::vector<boost::shared_ptr<Matrix> >& C_left, 
       std::vector<boost::shared_ptr<Matrix> >& C_right,
       boost::shared_ptr<BasisSet> primary);

    /**
     * Symmetric Constructor
     * @param C_symm reference to std::vector that will hold 
     *        both left- and right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    JK(std::vector<boost::shared_ptr<Matrix> >& C_symm,
       boost::shared_ptr<BasisSet> primary);
    /// Destructor
    virtual ~JK();


    /** 
    * Static instance constructor, used to get prebuilt DFJK/DirectJK objects
    * using knobs in options. 
    * @param options Options reference, with preset parameters
    * @param lr_symmetric left/right symmetric or not? defaults to false
    * @return abstract JK object, tuned in with preset options
    */
    static boost::shared_ptr<JK> build_JK(Options& options, bool lr_symmetric = false);

    // => Knobs <= // 

    /**
     * Cutoff for individual contributions to the J/K matrices
     * Eventually we hope to use Schwarz/MBIE/Density cutoffs,
     * for now just Schwarz
     * @param cutoff ceiling of magnitude of elements to be
     *        ignored if possible
     */
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    /**
     * Maximum memory to use, in doubles (for tensor-based methods,
     * integral generation objects typically ignore this)
     * @param memory maximum number of doubles to allocate
     */
    void set_memory(unsigned long int memory) { memory_ = memory; }
    /**
     * Maximum number of OpenMP threads to use. It may be necessary
     * to clamp this to some value smaller than the total number of
     * cores for machines with a high core-to-memory ratio to avoid
     * running out of memory due to integral generation objects
     * @param omp_nthread Maximum number of threads to use in 
     *        integral generation objects (BLAS/LAPACK can still
     *        run with their original maximum number) 
     */ 
    void set_omp_nthread(int omp_nthread) { omp_nthread_ = omp_nthread; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 1)
    void set_debug(int debug) { debug_ = debug; }
    /**
    * Set to do J tasks
    * @param do_J do J matrices or not,
    *        defaults to true 
    */
    void set_do_J(bool do_J) { do_J_ = do_J; }
    /**
    * Set to do K tasks
    * @param do_K do K matrices or not,
    *        defaults to true 
    */
    void set_do_K(bool do_K) { do_K_ = do_K; }
    /**
    * Set to do wK tasks
    * @param do_wK do wK matrices or not,
    *        defaults to false
    */
    void set_do_wK(bool do_wK) { do_wK_ = do_wK; }
    /**
    * Set the omega value for wK
    * @param omega range-separation parameter
    */
    void set_omega(double omega) { omega_ = omega; }

    // => Computers <= //

    /** 
     * Initialize the integral technology.
     * MUST be called AFTER setting knobs
     * but BEFORE first call of compute()
     */
    void initialize();
    /**
     * Compute D/J/K for the current C
     * Update values in your reference to
     * C_left/C_right BEFORE calling this,
     * renew your references to the matrices
     * in D/J/K AFTER calling this.
     */
    void compute();
    /**
     * Method to clear off memory without
     * totally destroying the object. The 
     * object can be rebuilt later by calling
     * initialize()
     */
    void finalize();

    // => Accessors <= //

    /**
     * Reference to C_left queue. It is YOUR job to
     * allocate and fill this object out 
     */
    std::vector<boost::shared_ptr<Matrix> >& C_left() { return C_left_; } 
    /**
     * Reference to C_right queue. It is YOUR job to
     * allocate and fill this object out. Only fill
     * C_left if symmetric. 
     */
    std::vector<boost::shared_ptr<Matrix> >& C_right() { return C_right_; } 

    /**
     * Reference to J results. The reference to the 
     * std::vector<boost::shared_ptr<Matrix> > is valid
     * throughout the life of the object. However, the 
     * entries (actual boost::shared_ptr<Matrix> pointers)
     * may be changed in each call of compute();
     * @return J vector of J matrices
     */
    const std::vector<boost::shared_ptr<Matrix> >& J() const { return J_; }
    /**
     * Reference to K results. The reference to the 
     * std::vector<boost::shared_ptr<Matrix> > is valid
     * throughout the life of the object. However, the 
     * entries (actual boost::shared_ptr<Matrix> pointers)
     * may be changed in each call of compute();
     * @return K vector of K matrices
     */
    const std::vector<boost::shared_ptr<Matrix> >& K() const { return K_; }
    /**
     * Reference to wK results. The reference to the 
     * std::vector<boost::shared_ptr<Matrix> > is valid
     * throughout the life of the object. However, the 
     * entries (actual boost::shared_ptr<Matrix> pointers)
     * may be changed in each call of compute();
     * @return wK vector of wK matrices
     */
    const std::vector<boost::shared_ptr<Matrix> >& wK() const { return wK_; }
    /**
     * Reference to D results. The reference to the 
     * std::vector<boost::shared_ptr<Matrix> > is valid
     * throughout the life of the object. However, the 
     * entries (actual boost::shared_ptr<Matrix> pointers)
     * may be changed in each call of compute();
     * @return D vector of D matrices
     */
    const std::vector<boost::shared_ptr<Matrix> >& D() const { return D_; }

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const = 0;
};

// => APPLIED CLASSES <= //

/**
 * Class DirectJK
 *
 * JK implementation using sieved, threaded 
 * integral-direct technology
 *
 * Note: This class builds a TwoBodyAOInt for each OpenMP
 * thread, for thread safety. This might be a bad idea if
 * you have a high core-to-memory ratio. Clamp the
 * openmp_nthread value if this fate befalls you.
 */
class DirectJK : public JK {

protected:

    /// Integral objects 
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri_;
    /// Integral factory (must be retained for Spherical Transforms)
    boost::shared_ptr<IntegralFactory> factory_;
    /// ERI Sieve
    boost::shared_ptr<ERISieve> sieve_;

    // => Required Algorithm-Specific Methods <= //

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return true; }
    /// Setup integrals, files, etc
    virtual void preiterations(); 
    /// Compute J/K for current C/D 
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations(); 

    /// Common initialization
    void common_init();

public:
    // => Constructors < = //
        
    /**
     * Non-Symmetric Constructor
     * @param C_left reference to std::vector that will hold 
     *        left-side C matrices
     * @param C_right reference to std::vector that will hold 
     *        right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    DirectJK(std::vector<boost::shared_ptr<Matrix> >& C_left, 
       std::vector<boost::shared_ptr<Matrix> >& C_right,
       boost::shared_ptr<BasisSet> primary);
    /**
     * Symmetric Constructor
     * @param C_symm reference to std::vector that will hold 
     *        left- and right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    DirectJK(std::vector<boost::shared_ptr<Matrix> >& C_symm,
       boost::shared_ptr<BasisSet> primary);
    /// Destructor
    virtual ~DirectJK();

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};

/**
 * Class DFJK
 *
 * JK implementation using sieved, threaded 
 * density-fitted technology
 */
class DFJK : public JK {

protected:

    // => DF-Specific stuff <= //

    /// Auxiliary basis set
    boost::shared_ptr<BasisSet> auxiliary_;
    /// PSIO object
    boost::shared_ptr<PSIO> psio_; 
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;
    /// File number for (Q|mn) tensor
    unsigned int unit_; 
    /// Core or disk?
    bool is_core_;
    /// Maximum number of rows to handle at a time
    int max_rows_;
    /// Maximum number of nocc in C vectors 
    int max_nocc_;
    /// Sieve, must be static throughout the life of the object
    boost::shared_ptr<ERISieve> sieve_;
    
    // => Temps (built/destroyed in compute_JK) <= //

    boost::shared_ptr<Matrix> Qmn_;
    boost::shared_ptr<Vector> J_temp_; 
    boost::shared_ptr<Vector> D_temp_; 
    boost::shared_ptr<Vector> d_temp_; 
    boost::shared_ptr<Matrix> E_left_; 
    boost::shared_ptr<Matrix> E_right_; 
    std::vector<boost::shared_ptr<Matrix> > C_temp_; 
    std::vector<boost::shared_ptr<Matrix> > Q_temp_; 

    // => Required Algorithm-Specific Methods <= //

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return true; }
    /// Setup integrals, files, etc
    virtual void preiterations(); 
    /// Compute J/K for current C/D 
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations(); 

    /// Common initialization
    void common_init();

    bool is_core();
    unsigned long int memory_temp();
    int max_rows();
    int max_nocc();
    void initialize_temps();
    void free_temps();

    virtual void initialize_JK_core();
    virtual void initialize_JK_disk();
    virtual void manage_JK_core();
    virtual void manage_JK_disk();
    virtual void block_J(double** Qmnp, int naux);
    virtual void block_K(double** Qmnp, int naux);

public:
    // => Constructors < = //
        
    /**
     * Non-Symmetric Constructor
     * @param C_left reference to std::vector that will hold 
     *        left-side C matrices
     * @param C_right reference to std::vector that will hold 
     *        right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    DFJK(std::vector<boost::shared_ptr<Matrix> >& C_left, 
       std::vector<boost::shared_ptr<Matrix> >& C_right,
       boost::shared_ptr<BasisSet> primary,
       boost::shared_ptr<BasisSet> auxiliary);

    /**
     * Symmetric Constructor
     * @param C_symm reference to std::vector that will hold 
     *        left- and right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    DFJK(std::vector<boost::shared_ptr<Matrix> >& C_symm,
       boost::shared_ptr<BasisSet> primary,
       boost::shared_ptr<BasisSet> auxiliary);
    /// Destructor
    virtual ~DFJK();

    // => Knobs <= //

    /**
     * Minimum relative eigenvalue to retain in fitting inverse
     * All eigenvectors with \epsilon_i < condition * \epsilon_max
     * will be discarded
     * @param condition, minimum relative eigenvalue allowed,
     *        defaults to 1.0E-12
     */
    void set_condition(double condition) { condition_ = condition; }
    /**
     * Which file number should the (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit(unsigned int unit) { unit_ = unit; }
    
    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};

/**
 * Class GPUDFJK
 *
 * JK implementation using sieved, threaded 
 * density-fitted CUDA technology
 */
class GPUDFJK : public DFJK {

protected:
    /// Eugene overloads these guys. gogo CUDA!
    virtual void block_J(double** Qmnp, int naux);
    virtual void block_K(double** Qmnp, int naux);
    void common_init();
    /// Setup integrals, files, etc
    virtual void preiterations(); 
    /// Delete integrals, files, etc
    virtual void postiterations(); 

public:

    /**
     * Non-Symmetric Constructor
     * @param C_left reference to std::vector that will hold 
     *        left-side C matrices
     * @param C_right reference to std::vector that will hold 
     *        right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    GPUDFJK(std::vector<boost::shared_ptr<Matrix> >& C_left, 
       std::vector<boost::shared_ptr<Matrix> >& C_right,
       boost::shared_ptr<BasisSet> primary,
       boost::shared_ptr<BasisSet> auxiliary);

    /**
     * Symmetric Constructor
     * @param C_symm reference to std::vector that will hold 
     *        left- and right-side C matrices
     * @param primary primary basis set for this system. 
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    GPUDFJK(std::vector<boost::shared_ptr<Matrix> >& C_symm,
       boost::shared_ptr<BasisSet> primary,
       boost::shared_ptr<BasisSet> auxiliary);
    /// Destructor
    virtual ~GPUDFJK();

    /**
     * helper class definied in gpudfjkhelper.h
     */
    boost::shared_ptr<GPUDFJKHelper> helper_;

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};


}
#endif
