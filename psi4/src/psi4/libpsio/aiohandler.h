/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef AIOHANDLER_H
#define AIOHANDLER_H

#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>

#include "config.h"

namespace psi {

class PSIO;

class AIOHandler {
   private:
    /// What is the job type?
    std::queue<size_t> job_;
    /// Unique job ID to check for job completion. Should NEVER be 0.
    std::deque<size_t> jobID_;
    /// Unit number argument
    std::queue<size_t> unit_;
    /// Entry Key (80-char) argument
    std::queue<const char *> key_;
    /// Memory buffer argument
    std::queue<char *> buffer_;
    /// Size argument
    std::queue<size_t> size_;
    /// Start address argument
    std::queue<psio_address> start_;
    /// End address pointer argument
    std::queue<psio_address *> end_;
    /// Matrix pointer for discontinuous I/O
    std::queue<double **> matrix_;
    /// Size argument for discontinuous I/O
    std::queue<size_t> row_length_;
    /// Size argument for discontinuous I/O
    std::queue<size_t> col_length_;
    /// Size argument for discontinuous I/O
    std::queue<size_t> col_skip_;
    /// For IWL: number of ints in the buffer
    std::queue<int> nints_;
    /// For IWL: is this the last buffer ?
    std::queue<int> lastbuf_;
    /// For IWL: pointer to current position in file
    std::queue<size_t *> address_;
    /// PSIO object this AIO_Handler is built on
    std::shared_ptr<PSIO> psio_;
    /// Thread this AIO_Handler is currently running on
    std::shared_ptr<std::thread> thread_;
    /// Lock variable
    std::mutex *locked_;
    /// Latest unique job ID
    size_t uniqueID_;
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
    size_t read(size_t unit, const char *key, char *buffer, size_t size, psio_address start, psio_address *end);
    /// Asynchronous write, same as PSIO::write, but nonblocking
    size_t write(size_t unit, const char *key, char *buffer, size_t size, psio_address start, psio_address *end);
    /// Asynchronous read_entry, same as PSIO::read_entry, but nonblocking
    size_t read_entry(size_t unit, const char *key, char *buffer, size_t size);
    /// Asynchronous read_entry, same as PSIO::write_entry, but nonblocking
    size_t write_entry(size_t unit, const char *key, char *buffer, size_t size);
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
    size_t read_discont(size_t unit, const char *key, double **matrix, size_t row_length, size_t col_length,
                        size_t col_skip, psio_address start);
    /// Same as read_discont, but for writing
    size_t write_discont(size_t unit, const char *key, double **matrix, size_t row_length, size_t col_length,
                         size_t col_skip, psio_address start);

    /// Zero disk
    /// Fills a double precision disk entry with zeros
    /// Total fill size is rows*cols*sizeof(double)
    /// Buffer memory of cols*sizeof(double) is used
    size_t zero_disk(size_t unit, const char *key, size_t rows, size_t cols);

    /// Write IWL
    /// Write an IWL buffer, thus containing
    /// IWL_INTS_PER_BUF integrals, 4 labels per integral, plus one
    /// integer indicating whether it is the last buffer and one integer
    /// counting the number of integrals in the current buffer
    size_t write_iwl(size_t unit, const char *key, size_t nints, int lastbuf, char *labels, char *values,
                     size_t labsize, size_t valsize, size_t *address);
    /// Generic function bound to thread internally
    void call_aio();

    /// Function that checks if a job has been completed using the JobID.
    /// The function only returns when the job is completed.
    void wait_for_job(size_t jobid);
};

}  // namespace psi

#endif  // AIOHANDLER_H
