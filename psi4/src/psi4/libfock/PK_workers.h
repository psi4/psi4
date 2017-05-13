/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef PKWRKR_H
#define PKWRKR_H

#if !defined( EXPLICIT_IOFF )
#   define EXPLICIT_IOFF(i) ( (i) * ((i) + 1) / 2 )
#endif

#if !defined( INDEX2 )
#   define INDEX2(i, j) ( (i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i) )
#endif

#if !defined( INDEX4 )
#   define INDEX4(i, j, k, l) ( INDEX2( INDEX2((i), (j)), INDEX2((k), (l)) ) )
#endif

#include "psi4/libiwl/config.h"
#include "psi4/libpsio/config.h"

namespace psi {

class AIOHandler;
class ERISieve;
class BasisSet;
class GaussianShell;

namespace pk {

// Forward declarations just before definitions
class AOShellSieveIterator;
class AOFctSieveIterator;

typedef std::unique_ptr<AOShellSieveIterator> UniqueAOShellIt;
typedef std::shared_ptr<ERISieve> SharedSieve;

/** AOShellSieveIterator provides an iterator over significant shell
 * quartets using an ERISieve object.
 */

class AOShellSieveIterator {
   private:
    // Basis set
    std::shared_ptr<BasisSet> bs_;
    // Sieve object
    SharedSieve sieve_;
    // Vector of significant shell pairs
    const std::vector< std::pair<int, int> >& shell_pairs_;
    // Number of shell pairs
    size_t npairs_;
    // Shell triangular indices
    size_t PQ_, RS_;
    // Shell indices
    int P_, Q_, R_, S_;
    // Are we there yet?
    bool done_;

    // Populate indices from shell_pairs
    void populate_indices();
public:
    /// Constructor
    AOShellSieveIterator(std::shared_ptr< BasisSet > prim, SharedSieve sieve_input);

    /// Iterator functions
    void first();
    void next();
    bool is_done() { return done_; }

    /// Accessor functions
    int p() const { return P_; }
    int q() const { return Q_; }
    int r() const { return R_; }
    int s() const { return S_; }

    /// Sieved iterator over the basis functions of a given shell
    AOFctSieveIterator integrals_iterator();
};

/** AOFctSieveIterator: provides an iterator over significant functions for
 * a specific shell quartet, using an ERISieve object.
 */

class AOFctSieveIterator {
private:
    // Sieve
    std::shared_ptr<ERISieve> sieve_;
    // Integral indices
    int i_, j_, k_, l_;
    // Relative integral indices within shells
    int irel_, jrel_, krel_, lrel_;
    // Gaussian shells of interest
    const GaussianShell& usi_;
    const GaussianShell& usj_;
    const GaussianShell& usk_;
    const GaussianShell& usl_;
    // Number of functions for each shell
    int ni_, nj_, nk_, nl_;
    // First index for each shell
    int fi_, fj_, fk_, fl_;
    // Max index for each shell
    int maxi_, maxj_, maxk_, maxl_;
    // Are we there yet ?
    bool done_;

    // Do we have the same shells in the quartets?
    bool sh_aaaa_;
    bool sh_abab_;
    //Populate the integral indices, and reorder them
    void populate_indices();
    // Increment the integral indices to the next ones for ij
    void increment_bra();
    // Increment integral indices to the next ones for l,k,j,i (in this order)
    void increment_ket();
    // Reorder indices in canonical order
    void reorder_inds();
public:
    /// Constructor
    AOFctSieveIterator(const GaussianShell& s1, const GaussianShell& s2,
                       const GaussianShell& s3, const GaussianShell& s4,
                       std::shared_ptr<ERISieve> siev);

    /// Iterator functions
    void first();
    void next();
    bool is_done() { return done_; }

    int i() const { return i_; }
    int j() const { return j_; }
    int k() const { return k_; }
    int l() const { return l_; }
};


/** IWLAsync_PK: buffer class for pre-sorted PK integrals. Each buffer class
 * has one IWL bucket for each PK bucket. We pass integrals with their labels
 * to this class, and it stores it in the appropriate bucket. When a bucket is full,
 * it is dumped to an IWL file using asynchronous I/O for storage.
 * When we eventually switch out of IWL to something else, this class will have to be
 * reworked, and also the integral reading from PK using IWL objects, but that is all.
 */
class IWLAsync_PK {
private:
    /// File number
    int itap_;
    /// Position in bytes for next write
    size_t* address_;
    /// Integral labels
    Label* labels_[2];
    /// Integral values
    Value* values_[2];
    /// Job Ids for AIO
    size_t JobID_[2];
    /// Number of integrals per buffer
    size_t ints_per_buf_;
    /// Current number of integrals in buffer
    size_t nints_;
    /// Is this the last buffer for PK bucket?
    int lastbuf_;
    /// Are we using buffer 1 or 2?
    int idx_;
    /// The AIO Handler
    std::shared_ptr<AIOHandler> AIO_;

public:
    /// Constructor, also allocates the arrays
    IWLAsync_PK(size_t* address, std::shared_ptr<AIOHandler> AIO, int itap);
    /// Destructor, also deallocates the arrays
    ~IWLAsync_PK();

    /// Filling values in the bucket
    void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    /// Popping a value from the current buffer, also decrements integral count
    void pop_value(double &val, size_t &i, size_t &j, size_t &k, size_t &l);
    /// Actually writing integrals from the buffer to the disk.
    void write();
    /// Filling buffer with dummy values and flushing it. Also indicates
    /// that this is the last buffer.
    void flush();

    /// Accessor functions
    size_t nints() { return nints_; }
    size_t maxints() { return ints_per_buf_; }

};

/** PK_worker: Base class for PK workers, objects which take care
 * of storing integrals in the appropriate buffers/files during
 * the various algorithms handling PK supermatrix formation.
 * This class only provides a common implementation for relevant
 * derived workers.
 */

class PKWorker {
private:

    /// Current basis set
    std::shared_ptr<BasisSet> primary_;
    /// Current sieve
    SharedSieve sieve_;
    /// Are we doing wK?
    bool do_wK_;

    /// AIOHandler
    std::shared_ptr<AIOHandler> AIO_;
    /// File to write to
    int target_file_;

    /// Iterator over basis functions within a shell quartet
    UniqueAOShellIt shelliter_;
    /// Current global index of the buffer
    unsigned int bufidx_;
    /// Current offset
    size_t offset_;
    /// Current max ijkl index in the buffer
    size_t max_idx_;
    /// Size of one buffer
    size_t buf_size_;
    /// Number of buffers in the worker
    unsigned int nbuf_;
    /// Are there any shells left ?
    bool shells_left_;

    /// Indices of the current shell quartet
    unsigned int P_, Q_, R_, S_;

    /// Is the current shell relevant to the current worker ?
    bool is_shell_relevant();

    // This class should never be copied
    PKWorker(const PKWorker &other) {}
    PKWorker & operator = (PKWorker &other) {}

protected:
    /// Setter function for nbuf_
    void set_nbuf(unsigned int tmp) {nbuf_ = tmp; }
    /// Setting the buffer size, changes for wK
    void set_bufsize(size_t len) { buf_size_ = len; }
    /// Setter function for max_idx
    void set_max_idx(size_t tmp) { max_idx_ = tmp; }

public:
    /// Constructor for PKWorker
    PKWorker(std::shared_ptr<BasisSet> primary, SharedSieve sieve,
             std::shared_ptr<AIOHandler> AIO, int target_file,
             size_t buf_size);
    /// Destructor for PKWorker, does nothing
    virtual ~PKWorker() {}

    /// Accessor functions
    std::shared_ptr< AIOHandler > AIO() const { return AIO_; }
    size_t nbuf()                       const { return nbuf_; }
    size_t buf_size()                   const { return buf_size_; }
    size_t max_idx()                    const { return max_idx_; }
    size_t offset()                     const { return offset_; }
    int target_file()                   const { return target_file_; }
    unsigned int bufidx()               const { return bufidx_; }
    unsigned int P()                    const { return P_; }
    unsigned int Q()                    const { return Q_; }
    unsigned int R()                    const { return R_; }
    unsigned int S()                    const { return S_; }
    bool do_wK()                        const { return do_wK_; }
    /// Set do_wK
    void set_do_wK(bool tmp) { do_wK_ = tmp; }

    /// Get TOC labels for J or K
    static char* get_label_J(const int batch);
    static char* get_label_K(const int batch);
    static char* get_label_wK(const int batch);

    /// Get max ijkl index included in current buffer
    /// Overloaded by specific derived classes
    virtual void initialize_task()=0;

    /// Set up the first shell quartet to be computed
    void first_quartet(size_t i);
    /// Is there a shell quartet left to compute ?
    bool more_work() { return shells_left_; }
    /// Get the next shell quartet for the current worker
    void next_quartet();

    /// Reallocate the buffer memory for wK
    virtual void allocate_wK(size_t buf_size, unsigned int buf_per_thread) {
        throw PSIEXCEPTION("Function allocate_wK not implemented for this PK algorithm.\n");
    }
    /// For IWL, we need different arguments
    virtual void allocate_wK(std::shared_ptr<std::vector<size_t>> pos, int wKfile) {
        throw PSIEXCEPTION("Invoked function allocate_wK not implemented for this PK algorithm.\n");
    }

    /// Filling integral values in the relevant buffers
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l) = 0;
    /// Filling integral values in the relevant buffers for wK
    virtual void fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l) = 0;

    /// Writing integral arrays for storage
    virtual void write(std::vector< size_t > min_ind, std::vector< size_t > max_ind,
                       size_t pk_pairs) = 0;
    /// Writing wK integral arrays for storage
    virtual void write_wK(std::vector< size_t > min_ind, std::vector< size_t > max_ind,
                       size_t pk_pairs) {
        throw PSIEXCEPTION("Function write_wK not implemented for current PK algorithm\n");
    }

    /// For in-core: only finalize the integral array
    virtual void finalize_ints(size_t pk_pairs) {
        throw PSIEXCEPTION("Function not implemented for this PK algorithm.\n");
    }
    /// For in-core: only finalize the wK integral array
    virtual void finalize_ints_wK(size_t pk_pairs) {
        throw PSIEXCEPTION("Function not implemented for this PK algorithm.\n");
    }

    /// Functions specific to disk pre-sorting of integrals
    virtual bool pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
        throw PSIEXCEPTION("Function pop_value not implemented for this class\n");
    }
    /// Functions specific to disk pre-sorting of integrals
    virtual bool pop_value_wK(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) {
        throw PSIEXCEPTION("Function pop_value_wK not implemented for this class\n");
    }

    virtual void insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l) {
        throw PSIEXCEPTION("Function insert_value not implemented for this class\n");
    }
    virtual void insert_value_wK(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l) {
        throw PSIEXCEPTION("Function insert_value_wK not implemented for this class\n");
    }

    /// Flush a buffer to disk
    virtual void flush() {
        throw PSIEXCEPTION("Function flush not implemented for this class\n");
    }
    virtual void flush_wK() {
        throw PSIEXCEPTION("Function flush not implemented for this class\n");
    }
};

/** class PKWrkReord: This worker class is associated with the
 * Reorder algorithm, which computes integrals in the order needed
 * to write them directly to the PK supermatrix. The worker provides
 * the shells needed for a given interval of canonical indices ijkl,
 * which are stored in a single buffer. Once the buffer is full, it gets
 * written in the appropriate entries of the PK file.
 * All integrals needed for a buffer are computed by a single thread
 * to avoid using atomic operations on a big buffer.
 * We thus need enough buffers for good load balancing. More buffers
 * also means more integrals recomputed, though.
 * We may want to improve on this aspect.
 */

class PKWrkrReord : public PKWorker {
private:
    /// This class can have a variable number of internal buffers,
    /// thus we need to have a vector of all relevant quantities.
    /// In addition, each buffer can write to more than one PK entries

    /// TOC entry labels
    std::vector< std::vector<char*> > labels_J_;
    std::vector< std::vector<char*> > labels_K_;
    std::vector< std::vector<char*> > labels_wK_;
    /// Job IDs
    std::vector< std::vector<size_t> > jobID_J_;
    std::vector< std::vector<size_t> > jobID_K_;
    std::vector< std::vector<size_t> > jobID_wK_;
    /// The internal buffers themselves
    std::vector< double* > J_bufs_;
    std::vector< double* > K_bufs_;
    std::vector< double* > wK_bufs_;

    /// Dummy psio_address for storage of write return value
    psio_address dummy_;

    /// Internal buffer index
    unsigned int buf_;

    virtual void initialize_task();

public:
    /// Constructor
    PKWrkrReord(std::shared_ptr<BasisSet> primary, SharedSieve sieve,
                std::shared_ptr<AIOHandler> AIO, int target_file,
                size_t buffer_size, unsigned int nbuffer);
    /// Destructor
    ~PKWrkrReord();

    /// Reallocating memory for wK
    /// We make sure the deallocated buffers have been written to disk
    virtual void allocate_wK(size_t buf_size, unsigned int buf_per_thread);

    /// Filling integral values in relevant buffer
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    /// Filling wK integrals in relevant buffer
    virtual void fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l);

    /// Writing values in the appropriate PK file entry
    virtual void write(std::vector<size_t> min_ind,
                       std::vector<size_t> max_ind, size_t pk_pairs);
    /// Writing wK values in the appropriate PK file entry
    virtual void write_wK(std::vector<size_t> min_ind,
                       std::vector<size_t> max_ind, size_t pk_pairs);
};

/** class PKWrkInCore: Computes all integrals for PK supermatrix
 * and stores them in core. To avoid atomic access to the giant array storing the matrix,
 * it is subdivided in nthreads_ buffers. Integrals are appropriately reordered such that
 * each thread computes the shell quartets needed for its part of the buffer.
 * There are only nthreads tasks so this algorithm likely suffers from load imbalance,
 * although preliminary tests in Vtune seem fine.
 */

class PKWrkrInCore : public PKWorker {
private:
    // We need to know the number of workers to properly
    // adjust the size of the last buffer
    int nworkers_;
    size_t pk_size_;
    size_t last_buf_;
    // Memory allocated and deleted outside this worker
    // Pointers to absolute memory start
    double* J_buf0_;
    double* K_buf0_;
    double* wK_buf0_;
    // Pointers to local memory start
    double* J_bufp_;
    double* K_bufp_;
    double* wK_bufp_;

    virtual void initialize_task();

public:
    PKWrkrInCore(std::shared_ptr<BasisSet> primary, SharedSieve sieve,
                 size_t buf_size, size_t lastbuf, double* Jbuf,
                 double* Kbuf, double* wKbuf, int nworkers);

    /// Filling values in the relevant part of the buffer
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    /// Filling values in the relevant part of the buffer for wK
    virtual void fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l);

    /// Finalize the buffer: divide by two diagonal elements
    virtual void finalize_ints(size_t pk_pairs);
    /// Finalize the buffer for wK: divide by two diagonal elements
    virtual void finalize_ints_wK(size_t pk_pairs);

    /// Function write is never used
    virtual void write(std::vector<size_t> min_ind, std::vector<size_t> max_ind, size_t pk_pairs) {
        throw PSIEXCEPTION("Function not implemented for in-core");
    }
};

/** Class for Yoshimine pre-sorting to obtain the PK supermatrix.
 * This class uses little buckets to pre-sort the integrals for each
 * thread. No communication between threads and asynchronous writing
 * to disk.
 */

class PKWrkrIWL : public PKWorker {
private:
    /// File for K IWL batches
    int K_file_;
    /// File for wK IWL batches
    int wK_file_;
    /// Vector mapping pq index to a bucket
    std::vector< int > buf_for_pq_;
    /// Vectors of IWL buffers for storage in the proper pre-sorting bucket
    std::vector< IWLAsync_PK* > IWL_J_;
    std::vector< IWLAsync_PK* > IWL_K_;
    std::vector< IWLAsync_PK* > IWL_wK_;
    /// Pointer to array with addresses in bytes where the next write
    /// should go for each bucket on file.
    std::shared_ptr<std::vector<size_t>> addresses_;
    std::shared_ptr<std::vector<size_t>> addresses_wK_;

public:
    /// Constructor
    PKWrkrIWL(std::shared_ptr<BasisSet> primary, SharedSieve sieve,
              std::shared_ptr<AIOHandler> AIOp, int targetfile, int K_file,
              size_t buf_size, std::vector< int > &bufforpq,
              std::shared_ptr<std::vector<size_t>> pos);
    /// Destructor
    ~PKWrkrIWL();

    /// Preparing for wK pre-sorting to file
    virtual void allocate_wK(std::shared_ptr<std::vector<size_t>> pos, int wKfile);

    /// Filling integrals in the appropriate buffers
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    /// Filling wK integrals in the appropriate buffers
    virtual void fill_values_wK(double val, size_t i, size_t j, size_t k, size_t l);
    /// Popping a value from a buffer to finalize writing
    virtual bool pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l);
    /// Inserting a value back into a buffer
    virtual void insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l);
    /// Popping a wK value from a buffer to finalize writing
    virtual bool pop_value_wK(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l);
    /// Inserting a wK value back into a buffer
    virtual void insert_value_wK(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l);
    /// Flushing all buffers for current worker
    virtual void flush();
    /// Flushing all wK buffers for current worker
    virtual void flush_wK();

    /// Functions that are not used here
    virtual void initialize_task() {
        throw PSIEXCEPTION("initialize_task not implemented for this class\n");
    }
    virtual void write(std::vector<size_t> min_ind, std::vector<size_t> max_ind,
                       size_t pk_pairs) {
        throw PSIEXCEPTION("write not implemented for this class\n");
    }
};

} // End namespace pk
} // End namespace psi

#endif
