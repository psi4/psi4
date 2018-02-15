/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef JK_H
#define JK_H

#include <vector>
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/dimension.h"

namespace psi {
class MinimalInterface;
class BasisSet;
class Matrix;
class ERISieve;
class TwoBodyAOInt;
class Options;
class PSIO;

namespace pk {
class PKManager;
}

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
 * \code
 *      // Constructor, Algorithm corresponds
 *      // to Type
 *      std::shared_ptr<JKType> jk(new JKType(
 *          basis, ...));
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
 *      // Initialize calls your derived class's preiterations member
 *      jk->initialize();
 *
 *      // Enter iterations or whatever
 *      // In each iteration:
 *
 *      // Clear and pack the C_left/C_right arrays
 *      // If C_left == C_right, only pack C_left
 *      // Make sure to declare your handles as references (&)
 *      std::vector<SharedMatrix>& C_left = jk->C_left();
 *      C_left.clear()
 *      C_left.push_back(Cocc)
 *      ...
 *
 *      // Let jk compute for the given C_left/C_right
 *      jk->compute();
 *
 *      // In return for the flexibility I give you with
 *      // C_left/C_right, I ask that you renew your reference
 *      // to the results J/K/D here. I only malloc where needed,
 *      // but if anything at all happens to C_left/C_right, last
 *      // iteration's SharedMatrix values might not
 *      // be valid. The std::vector<boost::shard_ptr<Matrix> >
 *      // for J/K/D will still be valid, but the pointers in the
 *      // entries may change.
 *
 *      SharedMatrix Jnew = jk->J()[0];
 *      SharedMatrix Knew = jk->K()[0];
 *
 *      // After iterations:
 *
 *      // Finalize (frees most memory, including C/D/J/K
 *      // pointers you are no longer holding)
 *      jk->finalize();
 *      \endcode
 *
 *  If you're like me and are trying to implement a new flavor of JK
 *  object, you probably want to know what sorts of things it needs to
 *  do and what hooks are
 *  At least within the confines of the HF code, the JK object
 *  works off of a modified visitor pattern.
 *  The derived wavefunction class, HF, contains a pointer to an instance
 *  of a JK object.  At certain hooks (vide infra) throughout the SCF
 *  the boilerplate HF code calls the overloaded functions of the JK
 *  object.  These are your chances to interact with the HF code.
 *
 *  1) During HF::integrals the JK instance is built by calling
 *     JK::build_JK, which is a static factory function.  This is your
 *     chance to set up your derived JK as you want it, via its
 *     constructor.
 *
 *  2) Also within HF::integrals, JK::initialize is called, which is
 *     a wrapper to an abstract function preiterations().  Best I can
 *     tell this a second chance to setup your derived JK object.
 *
 *  3) The last thing HF::integrals does is call JK::print_header, which
 *     is an abstract function that your derived class is supposed to
 *     implement so that it, you guessed it, prints a header.
 *
 *  4) Now the SCF proceeds until HF::FormG() (an abstract function) is
 *     called.  Because this function is implemented differently for
 *     each flavor of SCF, there is the possibility that different calls
 *     are made to your instance; however, all the calls seem to follow
 *     the same ordering:
 *
 *     4a) A call to JK::C_left, which allows the SCF code to set the
 *         MO coefficients
 *
 *     4b) A call to JK::compute() which:
 *         4b1) Ensures you have both a C_left and C_right
 *         4b2) Ensures you have a density
 *         4b3) Allocates memory in J_ and K_
 *         4b4) Calls compute_JK(), which is the next hook and is where
 *              you build J and K.  The resulting J and K are to be
 *              placed in the J_ and K_ members.
 *     4c) The SCF routine then grabs the J and K via getters to the
 *         base JK
 *  5) At this point the SCF decides if it's looping, which if it does
 *     return to point 4, otherwise...
 *
 *  6) You're off the hook! Do your cleanup in your destructor
 *
 *
 *  General notes:
 *
 *   - Best I can tell, abstract member JK::postiterations is never
 *     called (at least in RHF).  This means if you have actually used
 *     it to tear down your object, you'll need to call it in your
 *     destructor
 *       - Despite not being called, I suspect most JK flavors are fine
 *         and don't cause memory leaks owing to the use of shared
 *         pointers everywhere
 *       - Perhaps JK::~JK() should call it?
 *    - Your derived class is not responsible for memory associated with
 *      the density, MO coefficients, J, or K.  These things are
 *      allocated by the base class into shared pointers and will be
 *      autodestroyed when the pointers go out of scope
 *
 *
 */
class PSI_API JK {
   protected:
    // => Utility Variables <= //

    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    size_t memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, Process::environment.get_n_threads() otherwise)
    int omp_nthread_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;
    /// Whether to all desymmetrization, for cases when it's already been performed elsewhere
    std::vector<bool> input_symmetry_cast_map_;

    // => Tasks <= //

    /// Do J matrices? Defaults to true
    bool do_J_;
    /// Do K matrices? Defaults to true
    bool do_K_;
    /// Do wK matrices? Defaults to false
    bool do_wK_;

    /// Omega, defaults to 0.0
    double omega_;

    /// Left-right symmetric? Determined in each call of compute()
    bool lr_symmetric_;

    // => Architecture-Level State Variables (Spatial Symmetry) <= //

    /// Pseudo-occupied C matrices, left side
    std::vector<SharedMatrix> C_left_;
    /// Pseudo-occupied C matrices, right side
    std::vector<SharedMatrix> C_right_;
    /// Pseudo-density matrices \f$D_{ls}=C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedMatrix> D_;
    /// J matrices: \f$J_{mn}=(mn|ls)C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedMatrix> J_;
    /// K matrices: \f$K_{mn}=(ml|ns)C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedMatrix> K_;
    /// wK matrices: \f$K_{mn}(\omega)=(ml|\omega|ns)C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedMatrix> wK_;

    // => Microarchitecture-Level State Variables (No Spatial Symmetry) <= //

    /// Primary basis set
    std::shared_ptr<BasisSet> primary_;
    /// AO2USO transformation matrix
    SharedMatrix AO2USO_;
    /// Pseudo-occupied C matrices, left side
    std::vector<SharedMatrix> C_left_ao_;
    /// Pseudo-occupied C matrices, right side
    std::vector<SharedMatrix> C_right_ao_;
    /// Pseudo-density matrices
    std::vector<SharedMatrix> D_ao_;
    /// J matrices: J_mn = (mn|ls) C_li^left C_si^right
    std::vector<SharedMatrix> J_ao_;
    /// K matrices: K_mn = (ml|ns) C_li^left C_si^right
    std::vector<SharedMatrix> K_ao_;
    /// wK matrices: wK_mn = (ml|w|ns) C_li^left C_si^right
    std::vector<SharedMatrix> wK_ao_;

    // => Per-Iteration Setup/Finalize Routines <= //

    /// Build the pseudo-density D_, before compute_JK()
    void compute_D();
    /// Transform current C_left_/C_right_/D_ to C_left_ao_/C_right_ao_/D_ao_, before compute_JK()
    void USO2AO();
    /// Transform finished J_ao_/K_ao_ to J_/K_, after compute_JK()
    void AO2USO();
    /// Allocate J_/K_ should we be using SOs
    void allocate_JK();
    /**
     *  Function that sets a number of flags and allocates memory
     *  and sets up AO2USO.
     *
     *  Important note: no other memory is allocated here. That
     *  needs to be done in your derived class' constructor!!!!
     *
     *  Warning: this function is currently shadowed in at least
     *  one derived class, use with care!!!!!
     *
     */
    void common_init();

    // => Required Algorithm-Specific Methods <= //

    /// Setup integrals, files, etc
    virtual void preiterations() = 0;
    /// Compute J/K for current C/D
    virtual void compute_JK() = 0;
    /// Delete integrals, files, etc
    virtual void postiterations() = 0;

    // => Helper Routines <= //

    /// Memory (doubles) used to hold J/K/wK/C/D and ao versions, at current moment
    size_t memory_overhead() const;

   public:
    // => Constructors <= //

    /**
     *   Wrapper to common_init(), see that reference
     *
     *
     *
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     *
     *
     */
    JK(std::shared_ptr<BasisSet> primary);

    /// Destructor
    virtual ~JK();

    /**
    * Static instance constructor, used to get prebuilt DFJK/DirectJK objects
    * using knobs in options.
    * Nmat and sym are options for GTFock
    * sym means that all density matrices will be symmetric
    * @return abstract JK object, tuned in with preset options
    */
    static std::shared_ptr<JK> build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                        Options& options);
    static std::shared_ptr<JK> build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                        Options& options, std::string jk_type);

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const = 0;

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
    void set_memory(size_t memory) { memory_ = memory; }
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
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
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

    /**
     * Virtual method to provide (ia|ia) integrals for
     * SO-basis C_mi and C_na matrices in O(N^4) or less
     * Only available in DF-type JK integrals
     * Throws by default
     */
    virtual SharedVector iaia(SharedMatrix Ci, SharedMatrix Ca);

    // => Accessors <= //

    /**
     * Reference to C_left queue. It is YOUR job to
     * allocate and fill this object out
     */
    std::vector<SharedMatrix>& C_left() { return C_left_; }
    /**
     * Reference to C_right queue. It is YOUR job to
     * allocate and fill this object out. Only fill
     * C_left if symmetric.
     */
    std::vector<SharedMatrix>& C_right() { return C_right_; }

    /**
     * Reference to J results. The reference to the
     * std::vector<SharedMatrix > is valid
     * throughout the life of the object. However, the
     * entries (actual SharedMatrix pointers)
     * may be changed in each call of compute();
     * @return J vector of J matrices
     */
    const std::vector<SharedMatrix>& J() const { return J_; }
    /**
     * Reference to K results. The reference to the
     * std::vector<SharedMatrix > is valid
     * throughout the life of the object. However, the
     * entries (actual SharedMatrix pointers)
     * may be changed in each call of compute();
     * @return K vector of K matrices
     */
    const std::vector<SharedMatrix>& K() const { return K_; }
    /**
     * Reference to wK results. The reference to the
     * std::vector<SharedMatrix > is valid
     * throughout the life of the object. However, the
     * entries (actual SharedMatrix pointers)
     * may be changed in each call of compute();
     * @return wK vector of wK matrices
     */
    const std::vector<SharedMatrix>& wK() const { return wK_; }
    /**
     * Reference to D results. The reference to the
     * std::vector<SharedMatrix > is valid
     * throughout the life of the object. However, the
     * entries (actual SharedMatrix pointers)
     * may be changed in each call of compute();
     * @return D vector of D matrices
     */
    const std::vector<SharedMatrix>& D() const { return D_; }

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const = 0;
};

// => APPLIED CLASSES <= //

/**
 * Class DiskJK
 *
 * JK implementation using disk-based PK
 * integral technology
 */
class DiskJK : public JK {
    /// Absolute AO index to relative SO index
    int* so2index_;
    /// Absolute AO index to irrep
    int* so2symblk_;

    /// Options object
    Options& options_;

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return false; }
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
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    DiskJK(std::shared_ptr<BasisSet> primary, Options& options);
    /// Destructor
    virtual ~DiskJK();

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};

/**
 * Class PKJK
 *
 * JK implementation using disk-based
 * integral technology
 */
class PKJK : public JK {
    /// The PSIO instance to use for I/O
    std::shared_ptr<PSIO> psio_;

    /// Options object
    Options& options_;

    /// The pk file to use for storing the pk batches
    int pk_file_;

    /// The number of threads to be used for integral computation
    int nthreads_;

    /// Class handling the PK integrals
    std::shared_ptr<pk::PKManager> PKmanager_;

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const;
    /// Setup integrals, files, etc
    virtual void preiterations();
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations();

    /// Common initialization
    void common_init();

    /// Total number of SOs
    int nso_;
    /// Number of irreps
    int nirrep_;
    /// Number of so per irrep
    Dimension nsopi_;

   public:
    // => Constructors < = //

    /**
     * Symmetric Constructor
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    PKJK(std::shared_ptr<BasisSet> primary, Options& options);
    /// Destructor
    virtual ~PKJK();

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};

/**
 * Class DirectJK
 *
 * JK implementation using sieved, threaded
 * integral-direct technology
 *
 * Note: This class builds a TwoBodyAOInt for each OpenMP
 * thread, for thread safety. This might be a bad idea if
 * you have a high core-to-memory ratio. Clamp the
 * DF_INTS_NUM_THREADS value if this fate befalls you.
 */
class DirectJK : public JK {
   protected:
    /// Number of threads for DF integrals TODO: DF_INTS_NUM_THREADS
    int df_ints_num_threads_;
    /// ERI Sieve
    std::shared_ptr<ERISieve> sieve_;

    // => Required Algorithm-Specific Methods <= //

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return true; }
    /// Setup integrals, files, etc
    virtual void preiterations();
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations();

    /// Build the J and K matrices for this integral class
    void build_JK(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints, std::vector<std::shared_ptr<Matrix> >& D,
                  std::vector<std::shared_ptr<Matrix> >& J, std::vector<std::shared_ptr<Matrix> >& K);

    /// Common initialization
    void common_init();

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    DirectJK(std::shared_ptr<BasisSet> primary);
    /// Destructor
    virtual ~DirectJK();

    // => Knobs <= //

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};

/** \brief Derived class extending the JK object to GTFock
 *
 *   Unfortunately GTFock needs to know the number of density
 *   matrices and whether they are symmetric at construction.
 *   These two points are not within the design considerations of
 *   the base JK object and so this adds a slight complication if you
 *   want to use GTFock under those circumstances.  To get around
 *   this, you'll need to manually build a GTFockJK
 *   object and pass it into the constructor.  Don't worry
 *   building a GTFockJK object is easy, take a look at
 *   the Hartree-Fock code in HF.cc
 *
 */
class GTFockJK : public JK {
   private:
    /// The actual instance that does the implementing
    std::shared_ptr<MinimalInterface> Impl_;
    int NMats_ = 0;

   protected:
    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return true; }
    /// Setup integrals, files, etc
    virtual void preiterations() {}
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations() {}
    /// I don't fell the need to further clutter the output...
    virtual void print_header() const {}

   public:
    /** \brief Your public interface to GTFock
     *
     *  \param[in] Primary used by the base JK object, but not
     *         by GTFock.  Long term, this should be changed,
     *         but the reality is GTFock under the hood gets
     *         its basis in the same way as JK::build_JK gets
     *         Primary, so this shouldn't be an issue
     *  \param[in] NMats The number of density matrices you are
     *         passing in and consequently the number of Js and Ks
     *         you'll be getting back
     *  \param[in] AreSymm A flag specifying whether the density
     *         matrices you'll be passing in are symmetric.
     */
    GTFockJK(std::shared_ptr<psi::BasisSet> Primary, size_t NMats, bool AreSymm);
    /** \brief Your interface to GTFock that works well with libfock
    *   GTFock needs number of densities and symmetric at initialization
    *   This code calls GTFock once the number of densities was read from jk object
    */
    GTFockJK(std::shared_ptr<psi::BasisSet> Primary);
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
    std::shared_ptr<BasisSet> auxiliary_;
    /// PSIO object
    std::shared_ptr<PSIO> psio_;
    /// Cache action for three-index integrals
    std::string df_ints_io_;
    /// Number of threads for DF integrals
    int df_ints_num_threads_;
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;
    /// File number for (Q|mn) tensor
    size_t unit_;
    /// Core or disk?
    bool is_core_;
    /// Maximum number of rows to handle at a time
    int max_rows_;
    /// Maximum number of nocc in C vectors
    int max_nocc_;
    /// Sieve, must be static throughout the life of the object
    std::shared_ptr<ERISieve> sieve_;

    /// Main (Q|mn) Tensor (or chunk for disk-based)
    SharedMatrix Qmn_;
    /// (Q|P)^-1 (P|mn) for wK (or chunk for disk-based)
    SharedMatrix Qlmn_;
    /// (Q|w|mn) for wK (or chunk for disk-based)
    SharedMatrix Qrmn_;

    // => Temps (built/destroyed in compute_JK) <= //
    std::shared_ptr<Vector> J_temp_;
    std::shared_ptr<Vector> D_temp_;
    std::shared_ptr<Vector> d_temp_;

    SharedMatrix E_left_;
    SharedMatrix E_right_;
    std::vector<SharedMatrix> C_temp_;
    std::vector<SharedMatrix> Q_temp_;

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

    bool is_core() const;
    size_t memory_temp() const;
    int max_rows() const;
    int max_nocc() const;
    void initialize_temps();
    void free_temps();
    void initialize_w_temps();
    void free_w_temps();

    // => J <= //
    virtual void initialize_JK_core();
    virtual void initialize_JK_disk();
    virtual void manage_JK_core();
    virtual void manage_JK_disk();
    virtual void block_J(double** Qmnp, int naux);
    virtual void block_K(double** Qmnp, int naux);

    // => wK <= //
    virtual void initialize_wK_core();
    virtual void initialize_wK_disk();
    virtual void manage_wK_core();
    virtual void manage_wK_disk();
    virtual void block_wK(double** Qlmnp, double** Qrmnp, int naux);
    virtual void rebuild_wK_disk();

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param auxiliary auxiliary basis set for this system.
     */
    DFJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary);

    /// Destructor
    virtual ~DFJK();

    /**
     * Method to provide (ia|ia) integrals for
     * SO-basis C_mi and C_na matrices in O(N^4) or less
     * Only available in DF-type JK integrals
     * Throws by default
     */
    virtual SharedVector iaia(SharedMatrix Ci, SharedMatrix Ca);

    // => Knobs <= //

    /**
     * Minimum relative eigenvalue to retain in fitting inverse
     * All eigenvectors with \epsilon_i < condition * \epsilon_max
     * will be discarded
     * @param condition minimum relative eigenvalue allowed,
     *        defaults to 1.0E-12
     */
    void set_condition(double condition) { condition_ = condition; }
    /**
     * Which file number should the (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit(size_t unit) { unit_ = unit; }
    /**
     * What action to take for caching three-index integrals
     * @param val One of NONE, LOAD, or SAVE
     */
    void set_df_ints_io(const std::string& val) { df_ints_io_ = val; }
    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;
};
/**
 * Class CDJK
 *
 * JK implementation using
 * cholesky decomposition technology
 */
class CDJK : public DFJK {
   protected:
    // the number of cholesky vectors
    long int ncholesky_;

    // => Required Algorithm-Specific Methods <= //

    virtual bool is_core() { return true; }

    // => J <= //
    virtual void initialize_JK_core();
    virtual void initialize_JK_disk();
    virtual void manage_JK_core();

    double cholesky_tolerance_;

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const;

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     * @param cholesky_tolerance tolerance for cholesky decomposition.
     */
    CDJK(std::shared_ptr<BasisSet> primary, double cholesky_tolerance);

    /// Destructor
    virtual ~CDJK();
};
}
#endif
