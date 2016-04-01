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
#include <libpsio/aiohandler.h>

namespace psi {

/**
* IWLWriter functor for use with SO TEIs
**/
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

/**
  New buffer for doing asynchronous I/O on a single file.
  Multiple instances of the buffer can exist if they share the
  same AIO handler, and they can all write to the same file.
  If I understood AIO handler well, all writes are queued in one thread
  so it should all work out. These buffers have twice the size of the
  IWL buffers so that they can keep storing integrals while
  those in the first buffer are being written.
  **/

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

/**
* IWLAIOWriter  functor to use with SO TEIs, that is writing
* asynchronously to disk as the integrals get computed
**/

class IWLAIOWriter {
    IWLAsync& writeto_;
    size_t count_;

public:
    IWLAIOWriter(IWLAsync& writeto);

    void operator()(int i, int j, int k, int l, int, int, int, int,
                    int, int, int, int, double value);
    size_t count() const {return count_; }
};

/**
* ijklBasisIterator: Iterator that goes through all two-electron integral
* indices, in increasing order of canonical index ijkl.
* Used to determine what goes in which bucket for PK.
* Will be extended to include sieving.
**/

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

/** PK_integrals: Class that handles integrals for PK supermatrices in general.
 * It is used by the PKJK object to compute the number and size of buckets,
 * to load integrals in these buckets, and then to retrieve them for contraction
 * with density matrices.
 * This is adding a level of abstraction so that implementing sieving later may be
 * easier.
 * Should also make it easier to test different parallelization schemes for integral writing.
**/

class PK_integrals {
private:
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
    void fill_values(double val, int i, int j, int k, int l);
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
    char* label_J_[2];
    char* label_K_[2];

    int itap_J_;    // File number for J supermatrix
    int itap_K_;    // File number for K supermatrix
    size_t pk_pairs_;     // Dimension of pq pairs
    size_t pk_size_;      // Total dimension of supermatrix PK
    boost::shared_ptr<AIOHandler> AIO_;  /* AIO handler for all asynchronous operations */
    boost::shared_ptr<PSIO> psio_;       /* PSIO instance for opening/closing files */

    // The data for actually contracting density matrices with integrals
    std::vector<double*> D_vec_;

public:
    // Constructor
    PK_integrals(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<PSIO> psio,
                 int max_batches_, size_t memory);

    // Destructor
    ~PK_integrals();

    // Sizing the buckets
    void batch_sizing();
    // Printing the buckets
    void print_batches();
    // Initialize and allocate buffers
    void allocate_buffers();
    // Delete buffers
    void deallocate_buffers();
    // We fill each buffer using parallel threads, then write it when
    // all threads have joined. We then loop over the next buffer and restart parallel
    // threads. How many of such macroloops do we need ?
    int buf_ntasks() { return nbuffers_; }
    // Function that returns the number of quartets to compute in a given task, i.e.
    // the number of quartets to get all necessary integrals for a given ordered PK buffer
    // In addition, in the current implementation, the function allocates and stores the P,Q,
    // R, S indices of the shell in vectors.
    size_t task_quartets();
    // Sort computed integrals from IntegralFactory buffer to our PK buffers
    void integrals_buffering(double* buffer, int P, int Q, int R, int S);
    // Write the buffers of ordered integrals to disk
    void write();
    // We want to open the PK files for writing
    void open_files();
    // And then we want to close the files
    void close_files();

    // Accessor functions
    short int P(size_t idx) { return buf_P[idx]; }
    short int Q(size_t idx) { return buf_Q[idx]; }
    short int R(size_t idx) { return buf_R[idx]; }
    short int S(size_t idx) { return buf_S[idx]; }

    // Functions for contracting the density matrix with the integrals
    void form_D_vec(std::vector<SharedMatrix> D_ao);
    void form_J(std::vector<SharedMatrix> J, bool exch = false);
    void form_K(std::vector<SharedMatrix> K);
    // Simply deallocate the D vector
    void finalize_D();
};
}

#endif
