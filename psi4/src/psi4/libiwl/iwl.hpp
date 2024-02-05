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

#ifndef _psi_src_lib_libiwl_iwl_hpp_
#define _psi_src_lib_libiwl_iwl_hpp_

#include <cstdio>
#include "psi4/libpsio/psio.hpp"
#include "config.h"

namespace psi {

class PSI_API IWL {
    int itap_;            /* tape number for input file */
    psio_address bufpos_; /* current page/offset */
    int ints_per_buf_;    /* integrals per buffer */
    int bufszc_;          /* buffer size in characters (bytes) */
    double cutoff_;       /* cutoff value for writing */
    int lastbuf_;         /* is this the last IWL buffer? 1=yes,0=no */
    int inbuf_;           /* how many ints in current buffer? */
    int idx_;             /* index of integral in current buffer */
    Label *labels_;       /* pointer to where integral values begin */
    Value *values_;       /* integral values */
    /*! Instance of libpsio to use */
    PSIO *psio_;
    /*! Flag indicating whether to keep the IWL file or not */
    bool keep_;

   public:
    IWL();
    IWL(PSIO *psio, int itap, double cutoff, int oldfile, int readflag);
    ~IWL();

    // Accessor functions to data
    int &itap() { return itap_; }
    psio_address &buffer_position() { return bufpos_; }
    int &ints_per_buffer() { return ints_per_buf_; }
    int &buffer_size() { return bufszc_; }
    double &cutoff() { return cutoff_; }
    int &last_buffer() { return lastbuf_; }
    int &buffer_count() { return inbuf_; }
    int &index() { return idx_; }
    Label *labels() { return labels_; }
    Value *values() { return values_; }
    bool &keep() { return keep_; }

    void init(PSIO *psio, int itap, double cutoff, int oldfile, int readflag);

    void set_keep_flag(bool k) { keep_ = k; }
    void close();

    void fetch();
    void put();

    static void read_one(PSIO *psio, int itap, const char *label, double *ints, int ntri, int erase, int printflg,
                         std::string OutFileRMR);
    static void write_one(PSIO *psio, int itap, const char *label, int ntri, double *onel_ints);

    int read(int target_pq, double *ints, int *ioff_lt, int *ioff_rt, int mp2, int printflg, std::string OutFileRMR);

    void write(int p, int q, int pq, int pqsym, double *arr, int rmax, int *ioff, int *orbsym, int *firsti, int *lasti,
               int printflag, std::string OutFileRMR);
    void write_matrix(int ptr, int qtr, double **mat, int rfirst, int rlast, int sfirst, int slast, int *reorder,
                      int reorder_offset, int printflag, int *ioff, std::string OutFileRMR);
    void write_value(int p, int q, int r, int s, double value, int printflag, std::string OutFileRMR, int dirac);

    void flush(int lastbuf);
    void to_end();
};
}

#endif
