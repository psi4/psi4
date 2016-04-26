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

#ifndef PKWRKR_H
#define PKWRKR_H

#if !defined( INDEX2 )
#   define INDEX2(i, j) ( (i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i) )
#endif

#if !defined( INDEX4 )
#   define INDEX4(i, j, k, l) ( INDEX2( INDEX2((i), (j)), INDEX2((k), (l)) ) )
#endif

#include<libiwl/config.h>

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class AIOHandler;
class BasisSet;
class AOShellCombinationsIterator;

namespace pk {

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
    boost::shared_ptr<AIOHandler> AIO_;

public:
    /// Constructor, also allocates the arrays
    IWLAsync_PK(size_t* address, boost::shared_ptr<AIOHandler> AIO, int itap);
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
    boost::shared_ptr<BasisSet> primary_;

    /// AIOHandler
    boost::shared_ptr<AIOHandler> AIO_;
    /// File to write to
    int target_file_;

    /// Iterator over basis functions within a shell quartet
    AOShellCombinationsIterator shelliter_;
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

public:
    /// Constructor for PKWorker
    PKWorker(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
             int target_file, size_t buf_size);
    /// Destructor for PKWorker
    virtual ~PKWorker();

    /// Accessor functions
    size_t nbuf() { return nbuf_; }


    /// Get max ijkl index included in current buffer
    /// Overloaded by specific derived classes
    virtual size_t get_max_idx() = 0;

    /// Set up the first shell quartet to be computed
    void first_quartet(size_t i);
    /// Is there a shell quartet left to compute ?
    bool more_work() { return shells_left_; }
    /// Get the next shell quartet for the current worker
    void next_quartet();

    /// Filling integral values in the relevant buffers
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l) = 0;

    /// Writing integral arrays for storage
    virtual void write(std::vector< size_t > min_ind, std::vector< size_t > max_ind,
                       size_t pk_pairs) = 0;

    /// Functions specific to disk pre-sorting of integrals
    virtual void pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l) = 0;
    virtual void insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l)=0;
    /// Flush a buffer to disk
    virtual void flush();
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
    /// Job IDs
    std::vector< std::vector<size_t> > jobID_J_;
    std::vector< std::vector<size_t> > jobID_K_;
    /// The internal buffers themselves
    std::vector< double* > J_bufs_;
    std::vector< double* > K_bufs_;

    /// Internal buffer index
    unsigned int buf_;

    virtual size_t get_max_idx();

public:
    /// Constructor
    PKWrkrReord(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
                int target_file, size_t buf_size, unsigned int nbuf);
    /// Destructor
    ~PKWrkrReord();

    /// Filling integral values in relevant buffer
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);

    /// Writing values in the appropriate PK file entry
    virtual void write(std::vector<size_t> min_ind,
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
    size_t pk_size_;
    size_t last_buf_;
    // Memory allocated and deleted outside this worker
    double* J_bufp_;
    double* K_bufp_;

    virtual size_t get_max_idx();

public:
    PKWrkrInCore(boost::shared_ptr<BasisSet> primary, size_t buf_size, size_t lastbuf,
                 double* Jbuf, double* Kbuf);

    /// Filling values in the relevant part of the buffer
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);

    /// Finalize the buffer: divide by two diagonal elements
    /// We do not use the min_ind and max_ind vectors but that allows us
    /// to use the same call for both functions
    //TODO: Fix the above
    virtual void write(std::vector<size_t> min_ind,
                       std::vector<size_t> max_ind, size_t pk_pairs);
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
    /// Vector mapping pq index to a bucket
    std::vector< int > buf_for_pq_;
    /// Vectors of IWL buffers for storage in the proper pre-sorting bucket
    std::vector< IWLAsync_PK* > IWL_J_;
    std::vector< IWLAsync_PK* > IWL_K_;
    /// Pointer to array with addresses in bytes where the next write
    /// should go for each bucket on file.
    size_t * addresses_;

public:
    /// Constructor
    PKWrkrIWL(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<AIOHandler> AIO,
              int target_file, int K_file, size_t buf_size, std::vector< int > &bufforpq);
    /// Destructor
    ~PKWrkrIWL();

    /// Filling integrals in the appropriate buffers
    virtual void fill_values(double val, size_t i, size_t j, size_t k, size_t l);
    /// Popping a value from a buffer to finalize writing
    virtual bool pop_value(unsigned int bufid, double &val, size_t &i, size_t &j, size_t &k, size_t &l);
    /// Inserting a value back into a buffer
    virtual void insert_value(unsigned int bufid, double val, size_t i, size_t j, size_t k, size_t l);
    /// Flushing all buffers for current worker
    virtual void flush();



};

} // End namespace pk
} // End namespace psi

#endif
