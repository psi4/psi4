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

#ifndef WRITERS_H
#define WRITERS_H

#include <libiwl/iwl.hpp>
#include <libmints/integral.h>
#include <libmints/typedefs.h>
#include <libmints/sieve.h>
#include <libpsio/aiohandler.h>

namespace psi {

/*-
  IWLWriter functor for use with SO TEIs
-*/
class IWLWriter {
    IWL& writeto_;
    size_t count_;
    int& current_buffer_count_;

    Label *plabel_;
    Value *pvalue_;
public:

    IWLWriter(IWL& writeto);

    void operator()(int i, int j, int k, int l, int , int , int , int , int , int , int , int , double value);
    size_t count() const { return count_; }
};

/*-
  New buffer for doing asynchronous I/O on a single file.
  Multiple instances of the buffer can exist if they share the
  same AIO handler, and they can all write to the same file.
  If I understood AIO handler well, all writes are queued in one thread
  so it should all work out. These buffers have twice the size of the
  IWL buffers so that they can keep storing integrals while
  those in the first buffer are being written.
  -*/

class IWLAsync {
private:
    int itap_;                    /* File number */
    psio_address bufpos_;         /* We need to know where we are writing */
    int ints_per_buf_;            /* max integrals per buffer */
    int bufszc_;                  /* Size of the buffer in bytes */
    int lastbuf_[2];                 /* Is this the last buffer ? */
    int inbuf_[2];                   /* Number of integrals in the current buffer */
    int idx_;                     /* Index of integral in the current buffer */
    bool keep_;                   /* Whether or not to keep the file upon closing */
    Label *labels_[2];               /* Pointer to the array of four integral labels */
    Value *values_[2];               /* Pointer to the actual integral value */
    int whichbuf_;                /* Which one of the two buffers is currently written into */
    boost::shared_ptr<AIOHandler> AIO_;  /* AIO handler for all asynchronous operations */
    boost::shared_ptr<PSIO> psio_;       /* PSIO instance for opening/closing files */

    /* Map indicating whether the half-buffer pointed to by Value*
     * is currently being written to disk and thus should not be erased
     * before a synchronization of the writing threads
     */
    static std::map<Value*, bool> busy_;
    /* Map storing how many bytes are currently in the queue for
     * writing for each file using an aasynchronous IWL buffer.
     * This means that if, like I think, AIO Handler writes
     * everything in the order it is given tasks, we can compute the
     * starting place for the next write from the byte value stored here
     */
    static std::map<int, unsigned long int> bytes_written_;

public:
    // Constructor
    IWLAsync(boost::shared_ptr<PSIO> psio, boost::shared_ptr<AIOHandler> aio, int itap);
    // Destructor
    ~IWLAsync();

    // Accessor functions to data
    int& itap()                      {return itap_; }
    const int& ints_per_buffer()     {return ints_per_buf_; }
    const int& buffer_size()         {return bufszc_; }
    int& last_buffer()               {return lastbuf_[whichbuf_]; }
    int& buffer_count()              {return inbuf_[whichbuf_]; }
    int& index()                     {return idx_; }
    void set_keep_flag(bool flag)    {keep_ = flag; }

    Label get_p();
    Label get_q();
    Label get_r();
    Label get_s();

    Value get_val();

    void set_p(Label input);
    void set_q(Label input);
    void set_r(Label input);
    void set_s(Label input);

    void set_val(Value input);

    void fill_values(Label p, Label q, Label r, Label s, Value val);

    /// Asynchronously write the data to disk
    void put();
    /* Open the file: must be called only once if
     * multiple buffers for one file
     */
    void open_file(int oldfile);
    void flush(int lastbuf);
};

/*-
  IWLAIOWriter  functor to use with SO TEIs, that is writing
  asynchronously to disk as the integrals get computed
-*/

class IWLAIOWriter {
    IWLAsync& writeto_;
    size_t count_;

public:
    IWLAIOWriter(IWLAsync& writeto);

    void operator()(int i, int j, int k, int l, int, int, int, int,
                    int, int, int, int, double value);
    size_t count() const {return count_; }
};

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
public:
    // Constructor
    ijklBasisIterator(int nbas) : nbas_(nbas), done_(false) {}

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

/*- IOBuffer_PK: Class used as an I/O buffer for PK integrals. It stores
   integrals in local buffers. Each buffer contain sequential integrals in
   canonical order. Each buffer is associated with a list of shell quartets that
   are computed by a single thread (forming a single task). Once the task is done,
   the buffer is written to a pre-striped PK file using asynchronous I/O. Both J and
   K buffers are stored in here, and several of them can be present so that we don't
   need to wait on asynchronous I/O before getting a free buffer.
   -*/

class IOBuffer_PK {
private:
    // File to which we write
    int pk_file_;

    // TOC entry labels, need storage
    std::vector< std::vector<char*> > label_J_;
    std::vector< std::vector<char*> > label_K_;
    // Job IDs to wait on only the relevant job for AIO
    std::vector< std::vector<size_t> > jobID_J_;
    std::vector< std::vector<size_t> > jobID_K_;

    // Internal buffer index
    unsigned int buf_;

    // Vector of internal buffers
    std::vector<double*> J_buf_;
    std::vector<double*> K_buf_;

    // Size of internal buffers
    size_t buf_size_;
    // Number of internal buffers for each J and K
    unsigned int nbuf_;

    // Global buffer index, indicates where we are in canonical order
    unsigned int bufidx_;
    // Offset for current buffer
    size_t offset_;

    // Allocating function
    void allocate();
    // Deallocating function
    void deallocate();

    boost::shared_ptr<BasisSet> primary_;
    AOShellCombinationsIterator shelliter_;
    // Which shell quartet to compute ?
    unsigned int P_, Q_, R_, S_;
    // Does this shell quartet contain target integrals for current buffer ?
    bool is_shell_relevant();
    // Are there any shells left ?
    bool shells_left_;

    boost::shared_ptr<AIOHandler> AIO_;
    psio_address dummy_;

    // That should never be copied
    IOBuffer_PK(const IOBuffer_PK &other) {};
    IOBuffer_PK & operator = (IOBuffer_PK &other) {};
public:

    // Creator
    IOBuffer_PK(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
                size_t buf_size, size_t nbuf, int pk_file);
    //Destructor
    ~IOBuffer_PK();

    // Accessor functions to shell indices
    unsigned int P() { return P_; }
    unsigned int Q() { return Q_; }
    unsigned int R() { return R_; }
    unsigned int S() { return S_; }

    // Iterator functions through current buffer's quartets
    void first_quartet(size_t i);
    bool more_work() { return shells_left_; }
    void next_quartet();

    // Function to fill values in
    void fill_values(double val, size_t i, size_t j, size_t k, size_t l);

    // Task is done, write the buffer to file
    void write(std::vector<size_t> &batch_min_ind, std::vector<size_t> &batch_max_ind,
               size_t pk_pairs);

};


/*- PK_integrals: Class that handles integrals for PK supermatrices in general.
   It is used by the PKJK object to compute the number and size of buckets,
   to load integrals in these buckets, and then to retrieve them for contraction
   with density matrices.
   This is adding a level of abstraction so that implementing sieving later may be
   easier.
   Should also make it easier to test different parallelization schemes for integral writing.
-*/

class PK_integrals_old {
protected:
    int nbf_;  // Number of basis functions
    int max_batches_;  // Max number of buckets
    size_t memory_;
    boost::shared_ptr<BasisSet> primary_;

    // Buffers: hold integrals temporarily during their computation, before
    // storage on disk

    int nbuffers_;  // Number of buffers
    size_t buffer_size_;  // Size of each buffer
    double* J_buf_[2];       // Storage for J integrals
    double* K_buf_[2];       // Storage for K integrals
    int bufidx_;             // Which buffer are we filling?
    size_t offset_;          // Offset to write integrals to buffer
    unsigned short int buf_;  // Which one of the buffer pair is in use?
    std::vector<short int> buf_P; // Vector of P shell quartet indices for the current task
    std::vector<short int> buf_Q; // Vector of Q shell quartet indices for the current task
    std::vector<short int> buf_R; // Vector of R shell quartet indices for the current task
    std::vector<short int> buf_S; // Vector of S shell quartet indices for the current task
    void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    // Values of the AIOHandler job IDs
    unsigned long int jobid_J_[2];
    unsigned long int jobid_K_[2];


    // Batches: A batch of PK-ordered integrals, normally as large as memory
    // allows, for contraction with density matrices
    /// The index of the first pair in each batch
    std::vector<size_t> batch_pq_min_;
    /// The index of the last pair in each batch
    std::vector<size_t> batch_pq_max_;
    /// The index of the first integral in each batch
    std::vector<size_t> batch_index_min_;
    /// The index of the last integral in each batch
    std::vector<size_t> batch_index_max_;
    // This address stores the return value of write statements that we are
    // never going to use.
    psio_address dummy_;
    // The labels, need to be stored because we pass a pointer to them
    // to the AIOHandler
    std::vector<char*> label_J_[2];
    std::vector<char*> label_K_[2];
    // Sieving object for increased efficiency
    boost::shared_ptr<ERISieve> sieve_;

    int itap_J_;    // File number for J supermatrix
    int itap_K_;    // File number for K supermatrix
    size_t pk_pairs_;     // Dimension of pq pairs
    size_t pk_size_;      // Total dimension of supermatrix PK
    boost::shared_ptr<AIOHandler> AIO_;  /* AIO handler for all asynchronous operations */
    boost::shared_ptr<PSIO> psio_;       /* PSIO instance for opening/closing files */
    bool writing_; // Are we still doing AIO ?

    // The data for actually contracting density matrices with integrals
    std::vector<double*> D_vec_;

public:
    // Constructor
    PK_integrals_old(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<PSIO> psio,
                 int max_batches_, size_t memory, double cutoff);

    // Destructor
    virtual ~PK_integrals_old();

    // Sizing the buckets
    void batch_sizing();
    // Printing the buckets
    void print_batches();
    // Initialize and allocate buffers
    virtual void allocate_buffers();
    // Delete buffers
    virtual void deallocate_buffers();
    // We fill each buffer using parallel threads, then write it when
    // all threads have joined. We then loop over the next buffer and restart parallel
    // threads. How many of such macroloops do we need ?
    virtual size_t ntasks() { return nbuffers_; }
    // Function that returns the number of quartets to compute in a given task, i.e.
    // the number of quartets to get all necessary integrals for a given ordered PK buffer
    // In addition, in the current implementation, the function allocates and stores the P,Q,
    // R, S indices of the shell in vectors.
    size_t task_quartets();
    // Sort computed integrals from IntegralFactory buffer to our PK buffers
    virtual void integrals_buffering(const double* buffer, int P, int Q, int R, int S);
    // Write the buffers of ordered integrals to disk
    virtual void write();
    // We want to open the PK files for writing
    virtual void open_files(bool old);
    // And then we want to close the files
    virtual void close_files();

    // Getting J and K labels
    static char* get_label_J(size_t i);
    static char* get_label_K(size_t i);

    // Accessor functions
    short int P(size_t idx) { return buf_P[idx]; }
    short int Q(size_t idx) { return buf_Q[idx]; }
    short int R(size_t idx) { return buf_R[idx]; }
    short int S(size_t idx) { return buf_S[idx]; }
    bool writing()          { return writing_; }
    void set_writing(bool tmp) { writing_ = tmp; }

    // Functions for contracting the density matrix with the integrals
    void form_D_vec(std::vector<SharedMatrix> D_ao);
    virtual void form_J(std::vector<SharedMatrix> J, bool exch = false);
    virtual void form_K(std::vector<SharedMatrix> K);
    // Simply deallocate the D vector
    void finalize_D();
};

class PK_integrals : public PK_integrals_old {
private:
    int pk_file_;
    std::vector<IOBuffer_PK*> iobuffers_;
    int nthreads_;
    // Number of tasks per thread we would like for load balancing
    size_t tgt_tasks_;
    // Max memory per buffer, for fast allocation/writing
    size_t max_mem_buf_;
    // Total number of tasks we actually have
    size_t ntasks_;
public:
    // Constructor. We are pre-striping so everything can go to a single file
    // like in the good old times
    // Note however that PSIF_SO_PK actually contains AO integrals now...
    PK_integrals(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<PSIO> psio,
                 int max_batches, size_t memory, double cutoff);

    // Batch sizing is the same
    // Opening files for the first time: we pre-stripe
    virtual void open_files(bool old);
    // Closing file
    virtual void close_files();
    // Batch printing is the same
    // We need a new function for allocating batches
    // Size the batches first, then adapt to have enough tasks for good load balancing
    // It is not clear yet how many tasks are enough, will have to experiment
    virtual void allocate_buffers();
    // Finalize writing and deallocate the buffers
    virtual void deallocate_buffers();
    // Returns the total number of tasks we need
    virtual size_t ntasks() { return ntasks_; }
    // Returns a pointer to the current thread's IO buffer
    IOBuffer_PK* buffer();

    // Buffering the integrals before writing
    virtual void integrals_buffering(const double* buffer, unsigned int P,
                                     unsigned int Q, unsigned int R, unsigned int S);
    // Actually writing the integrals
    virtual void write();

};

}

#endif
