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
 * @END LICENSE
 */

#ifndef AIOHANDLER_H
#define AIOHANDLER_H

#include <thread>
#include <mutex>
#include <condition_variable>

namespace psi {

class AIOHandler {
private:
    /// What is the job type?
    std::queue<unsigned int> job_;
    /// Unique job ID to check for job completion. Should NEVER be 0.
    std::deque<unsigned long int> jobID_;
    /// Unit number argument
    std::queue<unsigned int> unit_;
    /// Entry Key (80-char) argument
    std::queue<const char*> key_;
    /// Memory buffer argument
    std::queue<char*> buffer_;
    /// Size argument
    std::queue<ULI> size_;
    /// Start address argument
    std::queue<psio_address> start_;
    /// End address pointer argument
    std::queue<psio_address*> end_;
    /// Matrix pointer for discontinuous I/O
    std::queue<double**> matrix_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> row_length_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> col_length_;
    /// Size argument for discontinuous I/O
    std::queue<ULI> col_skip_;
    /// For IWL: number of ints in the buffer
    std::queue<int> nints_;
    /// For IWL: is this the last buffer ?
    std::queue<int> lastbuf_;
    /// For IWL: pointer to current position in file
    std::queue<size_t*> address_;
    /// PSIO object this AIO_Handler is built on
    std::shared_ptr<PSIO> psio_;
    /// Thread this AIO_Handler is currently running on
    std::shared_ptr<std::thread> thread_;
    /// Lock variable
    std::mutex *locked_;
    /// Latest unique job ID
    unsigned long int uniqueID_;
    /// condition variable to wait for a specific job to finish
    std::condition_variable condition_;
public:
    /// AIO_Handlers are constructed around a synchronous PSIO object
    AIOHandler(std::shared_ptr<PSIO> psio);
    /// Destructor
    ~AIOHandler();
    /// When called, synchronize will not return until all requested data has been read or written
    void synchronize();
    /// Asynchronous read, same as PSIO::read, but nonblocking
    unsigned long int read(unsigned int unit, const char *key, char *buffer, ULI size,
              psio_address start, psio_address *end);
    /// Asynchronous write, same as PSIO::write, but nonblocking
    unsigned long int write(unsigned int unit, const char *key, char *buffer, ULI size,
               psio_address start, psio_address *end);
    /// Asynchronous read_entry, same as PSIO::read_entry, but nonblocking
    unsigned long int read_entry(unsigned int unit, const char *key, char *buffer, ULI size);
    /// Asynchronous read_entry, same as PSIO::write_entry, but nonblocking
    unsigned long int write_entry(unsigned int unit, const char *key, char *buffer, ULI size);
    /// Asynchronous read for reading discontinuous disk space
    /// into a continuous chunk of memory, i.e.
    ///
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]  <<--  [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    /// [***]        [----***-----]
    ///
    /// The buffer has dimensions row_length by col_length. The disk space
    /// has dimensions row_length by col_length + col_skip.
    ///
    /// These functions are not necessary for psio, but for aio they are.
    ///
    unsigned long read_discont(unsigned int unit, const char *key, double **matrix,
      ULI row_length, ULI col_length, ULI col_skip, psio_address start);
    /// Same as read_discont, but for writing
    unsigned long write_discont(unsigned int unit, const char *key, double **matrix,
      ULI row_length, ULI col_length, ULI col_skip, psio_address start);

    /// Zero disk
    /// Fills a double precision disk entry with zeros
    /// Total fill size is rows*cols*sizeof(double)
    /// Buffer memory of cols*sizeof(double) is used
    unsigned long zero_disk(unsigned int unit, const char* key, ULI rows, ULI cols);

    /// Write IWL
    /// Write an IWL buffer, thus containing
    /// IWL_INTS_PER_BUF integrals, 4 labels per integral, plus one
    /// integer indicating whether it is the last buffer and one integer
    /// counting the number of integrals in the current buffer
    unsigned long write_iwl(unsigned int unit, const char* key, size_t nints,
                            int lastbuf, char* labels, char* values, size_t labsize,
                            size_t valsize, size_t* address);
    /// Generic function bound to thread internally
    void call_aio();

    /// Function that checks if a job has been completed using the JobID.
    /// The function only returns when the job is completed.
    void wait_for_job(unsigned long int jobid);
};

}

#endif // AIOHANDLER_H
