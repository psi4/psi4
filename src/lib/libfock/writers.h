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
    bool set_keep_flag(bool flag)    {keep_ = flag; }

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

}

#endif
