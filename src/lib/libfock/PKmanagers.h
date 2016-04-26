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

#include <libmints/typedefs.h>


namespace psi {

// Forward declarations for Psi4
class Options;
class ERISieve;

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
    boost::shared_ptr<ERISieve> sieve_;
public:
    // Constructor
    ijklBasisIterator(int nbas, boost::shared_ptr<ERISieve> sieve) : nbas_(nbas),
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
 */

class PKManager {
private:
    Options& options_;
    /// Integral cutoff to apply
    double cutoff_;
    boost::shared_ptr<BasisSet> primary_;
    int nthreads_;
    /// Number of basis functions
    int nbf_;
    /// Number of tasks for integral computation
    size_t ntasks_;


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
    /// Vector of triangular result J/K matrices
    std::vector<double*> JK_vec_;

public:

    /// Base constructor
    PKManager(boost::shared_ptr<BasisSet> primary, size_t memory,
              Options& options);
    /// Base destructor
    virtual ~PKManager();

    /// Accessor functions for simple data
    double cutoff() const { return cutoff_; }
    int nthreads()  const { return nthreads_; }
    int nbf()       const { return nbf_; }
    int pk_file()   const { return pk_file_; }
    size_t pk_pairs() const { return pk_pairs_; }
    size_t pk_size()  const { return pk_size_; }
    size_t ntasks()   const { return ntasks_; }
    size_t memory()   const { return memory_; }

    /// Accessor that returns buffer corresponding to current thread
    SharedPKWrkr buffer();


    /**
     * @brief build_PKManager
     * Static instance constructor, used to get a proper
     * instance of PKManager through automatic selection and
     * options provided
     * @return abstract PKmanager object tuned with relevant options
     */
    static std::shared_ptr<PKManager> build_PKManager(boost::shared_ptr<PSIO> psio,
                boost::shared_ptr<BasisSet> primary, size_t memory, Options &options);

    // Base functions needed for the class to work
    /// Pure virtual initialize function: contains batch sizing and file
    /// opening if needed, allocation of thread buffers
    virtual void initialize()=0;
    /// Forming the PK supermatrices
    virtual void form_PK() = 0;
    /// Preparing JK computation
    virtual void prepare_JK(std::vector<SharedMatrix> D)=0;
    /// Cleaning up after JK computation
    virtual void finalize_JK();

    /// Function to print batch sizing information
    virtual void print_batches();
    /// Allocating the thread buffers
    virtual void allocate_buffers()=0;
    /// Actually computing the integrals
    virtual void compute_integrals();
    /// Deallocating the thread buffers
    virtual void finalize_PK()=0;
    /// Get the buffer for the current thread
    SharedPKWrkr get_buffer();


    /// Store the computed integrals in the appropriate buffers
    void integrals_buffering(const double *buffer, unsigned int P, unsigned int Q,
                             unsigned int R, unsigned int S);

    /// Write the buffers of integrals to PK storage
    virtual void write();

    /// Get TOC labels for J or K
    static char* get_label_J(const int batch);
    static char* get_label_K(const int batch);

    /// Actual computation of J and K
    /// Prepare the density matrix
    void form_D_vec(std::vector<SharedMatrix> D);
    /// Forming J
    virtual void form_J(std::vector<SharedMatrix> J, bool exch = false);
    /// Preparing triangular vector for J/K
    void make_J_vec(std::vector<SharedMatrix> J);
    /// Extracting results from vectors to matrix
    void get_results(std::vector<SharedMatrix> J);
    /// Forming K
    void form_K(std::vector<SharedMatrix> K);
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

    /// Maximum number of batches
    int max_batches_;

    /// PSIO handler (boost pointer for compatibility)
    boost::shared_ptr<PSIO> psio_;
    /// AIO Handler
    std::shared_ptr<AIOHandler> AIO_;
    /// PK file number
    int pk_file_;
    /// Is there any pending AIO writing ?
    bool writing_;

public:
    /// Constructor for PKMgrDisk
    PKMgrDisk(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
              size_t memory, Options &options);
    /// Destructor for PKMgrDisk
    virtual ~PKMgrDisk();

    /// Setter/Getter functions
    void set_writing(bool tmp) { writing_ = tmp; }
    bool writing()  const { return writing_; }

    /// Determining the batch sizes
    void batch_sizing();
    /// Printing out the batches
    virtual void print_batches();

    /// Opening files for integral storage
    virtual void prestripe_files()=0;
    /// Opening the PK file
    void open_PK_file();
    /// Closing the files
    virtual void close_PK_file(bool keep);

    /// Form J from PK supermatrix
    virtual void form_J(std::vector<SharedMatrix> J, bool exch);
};

class PKMgrReorder : public PKMgrDisk {
private:
    std::vector<char*> label_J_;
    std::vector<char*> label_K_;

    size_t max_mem_buf_;

public:
    /// Constructor
    PKMgrReorder(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
                 size_t memory, Options &options);
    /// Destructor
    virtual ~PKMgrReorder() {}

    /// Pre-striping the PK file
    virtual void prestripe_files();

    /// Allocating the buffers for each thread
    void allocate_buffers();
    /// Finalize PK: synchronize AIO writing, delete thread
    /// buffers, keep PK file open since we are going to
    /// use it immediately
    virtual void finalize_PK();

};

class PKMgrYoshimine : public PKMgrDisk {
private:
    /// Files of pre-sorted IWL integral buckets
    int iwl_file_J_;
    int iwl_file_K_;
    /// Number of integrals per buckets
    //TODO: change the name of all variables for buckets
    size_t ints_per_buf_;

    /// Total size of one IWL buffer on disk in bytes
    size_t iwl_int_size_;

public:
    /// Constructor
    PKMgrYoshimine(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
                   size_t memory, Options &options);
    /// Destructor
    virtual ~PKMgrYoshimine() {}

    /// Pre-striping the IWL pre-sorted bucket files
    //TODO: Optimize the vastly excessive amount of pre-striping for the K file
    virtual void prestripe_files();

    /// Gather all steps to form PK
    virtual void form_PK();

    /// Computing the integrals
    virtual void compute_integrals();

    /// Writing of the last partially filled buffers
    virtual void write();

    /// Reading and sorting integrals to generate PK file
    void sort_ints();

    /// Close the IWL bucket files
    void close_iwl_buckets();

    /// Generate the J PK supermatrix from IWL integrals
    void generate_J_PK(double* twoel_ints, size_t max_size);
    /// Generate the K PK supermatrix from IWL integrals
    void generate_K_PK(double* twoel_ints, size_t max_size);

};

/* PKMgrInCore: Class to manage in-core PK algorithm */

class PKMgrInCore : public PKManager {
private:
    /// Large in core arrays for integral storage
    std::unique_ptr<double> J_ints_;
    std::unique_ptr<double> K_ints_;

public:
    /// Constructor for in-core class
    PKMgrInCore(boost::shared_ptr<BasisSet> primary, size_t memory,
                Options &options) : PKManager(primary,memory,options) {}
    /// Destructor for in-core class
    virtual ~PKMgrInCore() {}

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
