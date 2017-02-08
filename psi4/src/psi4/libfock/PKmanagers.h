/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef PKMGR_H
#define PKMGR_H

//TODO Const correctness of everything

namespace psi {

// Forward declarations for Psi4
class Options;
class ERISieve;
class AIOHandler;
class BasisSet;

namespace pk {

class PKWorker;

typedef std::shared_ptr<PKWorker> SharedPKWrkr;

/*-
  ijklBasisIterator: Iterator that goes through all two-electron integral
  indices, in increasing order of canonical index ijkl.
  Used to determine what goes in which bucket for PK.
  Will be extended to include sieving.
-*/

class ijklBasisIterator {
private:
    int nbas_;
    int i_, j_, k_, l_;
    bool done_;
    std::shared_ptr<ERISieve> sieve_;
public:
    // Constructor
    ijklBasisIterator(int nbas, std::shared_ptr<ERISieve> sieve) : nbas_(nbas),
        done_(false), sieve_(sieve) {}

    // Iterator functions
    void first();
    void next();
    bool is_done() { return done_; }

    // Accessor functions
    int i() { return i_; }
    int j() { return j_; }
    int k() { return k_; }
    int l() { return l_; }
};

/**
 * Class PKManager
 *
 * Manages the different algorithms available for PK
 * and stores data shared by all algorithms.
 *
 * This class is abstract and specific instances must be
 * obtained by constructing an object corresponding to the
 * desired algorithm, or automatically selected through the
 * static build_PKManager
 *
 * For more details about each algorithm, see PK_workers.h
 */

class PKManager {
private:
    Options& options_;
    /// Integral cutoff to apply
    double cutoff_;
    std::shared_ptr<BasisSet> primary_;
    int nthreads_;
    /// Number of basis functions
    int nbf_;
    /// Number of tasks for integral computation
    size_t ntasks_;

    /// Do wK in addition to regular J and K?
    bool do_wK_;
    /// Value for omeaga
    double omega_;


    /// Sieving object for the integrals
    std::shared_ptr<ERISieve> sieve_;
    /// Size of triangular list of PK pairs
    size_t pk_pairs_;
    /// Total size of the four-index triangular PK
    size_t pk_size_;
    /// Total memory available in doubles
    size_t memory_;

    /// Array of IOBuffer_PK for task handling
    std::vector<SharedPKWrkr> iobuffers_;

    /// Vector of triangular D matrices
    std::vector<double*> D_vec_;
    /// Vector of original D matrices. We keep it around for building exchange
    /// in non-symmetric cases.
    std::vector<SharedMatrix> D_;
    /// Vector of booleans, true if corresponding density matrix is symmetric
    std::vector< bool > symmetric_;
    /// Are all density matrices symmetric?
    bool all_sym_;
    /// Vector of triangular result J/K matrices
    std::vector<double*> JK_vec_;

    /// Setter functions for internal wK options
    void set_wK(bool dowK) { do_wK_ = dowK; }
    void set_omega(double omega_in) { omega_ = omega_in; }

protected:
    /// Setter objects for internal data
    void fill_buffer(SharedPKWrkr tmp) { iobuffers_.push_back(tmp); }

public:

    /// Base constructor
    PKManager(std::shared_ptr<BasisSet> primary, size_t memory,
              Options& options);
    /// Base destructor, does nothing
    virtual ~PKManager() {}

    /// Accessor functions for simple data
    bool do_wk()                            const { return do_wK_; }
    double omega()                          const { return omega_; }
    double cutoff()                         const { return cutoff_; }
    int nthreads()                          const { return nthreads_; }
    int nbf()                               const { return nbf_; }
    std::shared_ptr< ERISieve > sieve()     const { return sieve_; }
    size_t pk_pairs()                       const { return pk_pairs_; }
    size_t pk_size()                        const { return pk_size_; }
    size_t ntasks()                         const { return ntasks_; }
    size_t memory()                         const { return memory_; }
    SharedPKWrkr buffer(int i)              const { return iobuffers_[i]; }
    double* D_glob_vecs(int i)              const { return D_vec_[i]; }
    double* JK_glob_vecs(int i)             const { return JK_vec_[i]; }
    std::shared_ptr< BasisSet > primary() const { return primary_; }
    bool is_sym(int i)                      const { return symmetric_[i]; }
    bool all_sym()                          const { return all_sym_; }
    SharedMatrix original_D(int N)          const { return D_[N]; }

    /// Accessor that returns buffer corresponding to current thread
    SharedPKWrkr get_buffer();
    void set_ntasks(size_t tmp) { ntasks_ = tmp; }

    /**
     * @brief build_PKManager
     * Static instance constructor, used to get a proper
     * instance of PKManager through automatic selection and
     * options provided
     * @return abstract PKmanager object tuned with relevant options
     */
    static std::shared_ptr<PKManager> build_PKManager(std::shared_ptr<PSIO> psio,
                std::shared_ptr<BasisSet> primary, size_t memory, Options &options,
                bool dowK, double omega_in = 0);

    // Base functions needed for the class to work
    /// Pure virtual initialize function: contains batch sizing and file
    /// opening if needed, allocation of thread buffers
    virtual void initialize()=0;
    /// Initialize for wK integrals
    virtual void initialize_wK()=0;
    /// Forming the PK supermatrices
    virtual void form_PK() = 0;
    /// Forming PK supermatrices for wK
    virtual void form_PK_wK() = 0;
    /// Preparing JK computation
    virtual void prepare_JK(std::vector<SharedMatrix> D,std::vector<SharedMatrix> Cl,
                            std::vector<SharedMatrix> Cr)=0;
    /// Cleaning up after JK computation
    virtual void finalize_JK() = 0;

    /// Function to print batch sizing information
    virtual void print_batches();
    /// Some printing for wK integrals;
    virtual void print_batches_wK() {}
    /// Allocating the thread buffers
    virtual void allocate_buffers()=0;
    /// Actually computing the integrals
    virtual void compute_integrals(bool wK = false);
    /// Computing integrals for wK
    virtual void compute_integrals_wK();
    /// Deallocating the thread buffers
    virtual void finalize_PK()=0;


    /// Store the computed integrals in the appropriate buffers
    void integrals_buffering(const double *buffer, unsigned int P, unsigned int Q,
                             unsigned int R, unsigned int S);
    /// Store the computed wK integrals in the appropriate buffers
    void integrals_buffering_wK(const double *buffer, unsigned int P, unsigned int Q,
                             unsigned int R, unsigned int S);

    /// Write the buffers of integrals to PK storage
    virtual void write() = 0;
    /// Write the buffers of wK integrals to PK storage
    virtual void write_wK() = 0;

    /// Actual computation of J and K
    /// Prepare the density matrix
    void form_D_vec(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cl, std::vector<SharedMatrix> Cr);
    /// Forming J, shared_ptr() initializes to null
    virtual void form_J(std::vector<SharedMatrix> J, std::string exch = "",
                        std::vector<SharedMatrix> K = std::vector<SharedMatrix>())=0;
    /// Preparing triangular vector for J/K
    void make_J_vec(std::vector<SharedMatrix> J);
    /// Extracting results from vectors to matrix
    void get_results(std::vector<SharedMatrix> J,std::string exch);
    /// Forming K
    void form_K(std::vector<SharedMatrix> K);
    /// Forming wK
    virtual void form_wK(std::vector<SharedMatrix> wK);
    /// Finalize and delete the density matrix vectors
    void finalize_D();
};

/* PKMgrDisk: Abstract base class to manage PK algorithms using disk I/O
 */

class PKMgrDisk : public PKManager {
private:
    /// The index of the first pair in each batch
    std::vector<size_t> batch_pq_min_;
    /// The index of the last pair in each batch
    std::vector<size_t> batch_pq_max_;
    /// The index of the first integral in each batch
    std::vector<size_t> batch_index_min_;
    /// The index of the last integral in each batch
    std::vector<size_t> batch_index_max_;
    /// Mapping pq indices to the correct batch
    std::vector<int> batch_for_pq_;

    /// Mapping the min and max pq's to p,q pairs
    std::map< size_t, std::pair< int, int> > ind_for_pq_;

    /// Maximum number of batches
    int max_batches_;

    /// PSIO handler (boost pointer for compatibility)
    std::shared_ptr<PSIO> psio_;
    /// AIO Handler
    std::shared_ptr<AIOHandler> AIO_;
    /// PK file number
    int pk_file_;
    /// Is there any pending AIO writing ?
    bool writing_;

public:
    /// Constructor for PKMgrDisk
    PKMgrDisk(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
              size_t memory, Options &options);
    /// Destructor for PKMgrDisk, does nothing
    virtual ~PKMgrDisk() {}

    /// Setter/Getter functions
    void set_writing(bool tmp) { writing_ = tmp; }
    bool writing()  const { return writing_; }
    int pk_file() const { return pk_file_; }
    std::vector< size_t >& batch_ind_min() { return batch_index_min_;}
    std::vector< size_t >& batch_ind_max() { return batch_index_max_;}
    std::vector< size_t >& batch_pq_min() { return batch_pq_min_;}
    std::vector< size_t >& batch_pq_max() { return batch_pq_max_;}
    std::vector< int >& batch_for_pq() { return batch_for_pq_; }
    std::shared_ptr< AIOHandler > AIO() const { return AIO_; }
    std::shared_ptr< PSIO > psio() const { return psio_; }

    /// Finalize the PK file formation
    virtual void finalize_PK();

    /// Initialize sequence for Disk algorithms
    virtual void initialize();
    /// Initialize wK PK supermatrix, has to be called after
    /// initialize()
    virtual void initialize_wK();
    /// Allocating new buffers for wK integrals
    virtual void allocate_buffers_wK()=0;
    /// Prepare the JK formation for disk algorithms
    virtual void prepare_JK(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cl,
                            std::vector<SharedMatrix> Cr);

    /// Determining the batch sizes
    void batch_sizing();
    /// Printing out the batches
    virtual void print_batches();

    /// Opening files for integral storage
    virtual void prestripe_files()=0;
    /// Opening files for wK integral storage
    virtual void prestripe_files_wK() = 0;

    /// Write integrals on disk
    virtual void write();
    /// Write wK integrals on disk
    virtual void write_wK();

    /// Opening the PK file
    void open_PK_file();
    /// Closing the files
    virtual void close_PK_file(bool keep);

    /// Form J from PK supermatrix, shared_ptr() initialized to null
    virtual void form_J(std::vector<SharedMatrix> J, std::string exch = "",
                        std::vector<SharedMatrix> K = std::vector<SharedMatrix>());

    /// Finalize JK matrix formation
    virtual void finalize_JK();
};

/**
 * Reorder algorithm: we divide the integrals in batches, as large as
 * possible. Then each thread computes one batch of ordered integrals,
 * meaning the thread gets assigned a list of shell quartets to go through
 * to obtain the relevant integrals.
 * Because of this division, some shell quartets will be recomputed several times.
 * For low memory, Yoshimine is recommended and selected automatically
 * by the build_PKManager constructor.
 *
 * This routine uses OMP multithreading. A dedicated thread generated by
 * AIOHandler handles all disk I/O
 */

class PKMgrReorder : public PKMgrDisk {
private:
    std::vector<char*> label_J_;
    std::vector<char*> label_K_;
    std::vector<char*> label_wK_;

    size_t max_mem_buf_;

public:
    /// Constructor
    PKMgrReorder(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
                 size_t memory, Options &options);
    /// Destructor
    virtual ~PKMgrReorder() {}

    /// Sequence of operations to form PK for reorder algo
    virtual void form_PK();
    /// Forming PK for wK integrals
    virtual void form_PK_wK();

    /// Pre-striping the PK file
    virtual void prestripe_files();
    /// Pre-striping PK file for wK integrals
    virtual void prestripe_files_wK();

    /// Allocating the buffers for each thread
    virtual void allocate_buffers();
    /// De-allocating buffer space for J/K supermatrices and
    /// allocating space for wK. Allows us to use buffers twice as big
    virtual void allocate_buffers_wK();
    /// Finalize PK: synchronize AIO writing, delete thread
    /// buffers, keep PK file open since we are going to
    /// use it immediately
    virtual void finalize_PK();

};

/** Yoshimine sorting algoritm: the integrals are divided in
 * N batches that can individually fit in memory. Then, all
 * integrals are computed only once, and sorted in N little buffers,
 * one for each batch. As soon as a buffer is full, it gets dumped
 * onto one of N temporary files on disk in IWL format.
 *
 * Once all integrals are computed and pre-sorted in the N
 * IWL files, each file is read and its integrals sorted,
 * then written onto the PK files.
 *
 * This routine takes advantage of OMP parallelization, then
 * each thread has N little buffers. All disk writing is handled
 * by a dedicated thread generated by AIOHandler
 */

class PKMgrYoshimine : public PKMgrDisk {
private:
    /// Files of pre-sorted IWL integral buckets
    int iwl_file_J_;
    int iwl_file_K_;
    int iwl_file_wK_;
    /// Number of integrals per buckets
    //TODO: change the name of all variables for buckets
    size_t ints_per_buf_;

    /// Total size of one IWL buffer on disk in bytes
    size_t iwl_int_size_;

public:
    /// Constructor
    PKMgrYoshimine(std::shared_ptr<PSIO> psio, std::shared_ptr<BasisSet> primary,
                   size_t memory, Options &options);
    /// Destructor
    virtual ~PKMgrYoshimine() {}

    /// Pre-striping the IWL pre-sorted bucket files
    //TODO: Optimize the vastly excessive amount of pre-striping for the K file
    virtual void prestripe_files();
    /// Pre-striping for IWL wK pre-sorted files
    virtual void prestripe_files_wK();

    /// Gather all steps to form PK
    virtual void form_PK();
    /// Steps to form the supermatrix for wK
    virtual void form_PK_wK();

    /// Allocating the buffers for each thread
    virtual void allocate_buffers();
    /// Allocating buffers for wK integrals
    virtual void allocate_buffers_wK();

    /// Computing the integrals
    virtual void compute_integrals(bool wK = false);
    /// computing wK integrals
    virtual void compute_integrals_wK();

    /// Writing of the last partially filled buffers
    virtual void write();
    /// Writing of the last partially filled buffers for wK
    virtual void write_wK();

    /// Reading and sorting integrals to generate PK file
    void sort_ints(bool wK = false);
    /// Reading and sorting wK integrals for PK file
    void sort_ints_wK();

    /// Close the IWL bucket files
    void close_iwl_buckets();
    /// Close the IWL bucket file for wK
    void close_iwl_buckets_wK();

    /// Generate the J PK supermatrix from IWL integrals
    void generate_J_PK(double* twoel_ints, size_t max_size);
    /// Generate the K PK supermatrix from IWL integrals
    void generate_K_PK(double* twoel_ints, size_t max_size);
    /// Generate the wK PK supermatrix from IWL integrals
    void generate_wK_PK(double* twoel_ints, size_t max_size);
};

/* PKMgrInCore: Class to manage in-core PK algorithm */

/** The simplest algorithm: a large buffer is allocated in core
 * and all integrals are written to it, directly sorted.
 *
 * This routine uses OMP multithreading and, for efficiency reasons,
 * each thread has a dedicated buffer and batch of integrals to compute for it.
 * There is thus some recomputation of the shell quartets, but we completely
 * avoid thread communication.
 */

class PKMgrInCore : public PKManager {
private:
    /// Large in core arrays for integral storage
    std::unique_ptr<double []> J_ints_;
    std::unique_ptr<double []> K_ints_;
    std::unique_ptr<double []> wK_ints_;

public:
    /// Constructor for in-core class
    PKMgrInCore(std::shared_ptr<BasisSet> primary, size_t memory,
                Options &options) : wK_ints_(nullptr), PKManager(primary,memory,options) {}
    /// Destructor for in-core class
    virtual ~PKMgrInCore() {}

    /// Initialize sequence for in-core algorithm
    virtual void initialize();
    /// Initialize the wK integrals
    virtual void initialize_wK();
    /// Sequence of steps to form PK matrix
    virtual void form_PK();
    /// Sequence of steps to form wK PK matrix
    virtual void form_PK_wK();
    /// Steps to prepare JK formation
    virtual void prepare_JK(std::vector<SharedMatrix> D,std::vector<SharedMatrix> Cl,
                            std::vector<SharedMatrix> Cr);

    /// Form J matrix, shared_ptr() initializes to null
    virtual void form_J(std::vector<SharedMatrix> J, std::string exch = "",
                        std::vector<SharedMatrix> K = std::vector<SharedMatrix>());
    /// Finalize JK formation
    virtual void finalize_JK();

    /// No disk write, only finalizes integral arrays
    virtual void write();
    /// Finalize integral arrays for wK
    virtual void write_wK();

    /// Printing the algorithm header
    virtual void print_batches();

    /// Allocate the buffer threads
    virtual void allocate_buffers();
    /// Finalize PK, i.e. deallocate buffers
    virtual void finalize_PK();

};

}
}

#endif
