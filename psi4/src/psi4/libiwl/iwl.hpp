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

#ifndef _psi_src_lib_libiwl_iwl_hpp_
#define _psi_src_lib_libiwl_iwl_hpp_

#include <cstdio>
#include "psi4/libpsio/psio.hpp"
#include "config.h"

namespace psi {

    class IWL {
        int itap_;                   /* tape number for input file */
        psio_address bufpos_;        /* current page/offset */
        int ints_per_buf_;           /* integrals per buffer */
        int bufszc_;                 /* buffer size in characters (bytes) */
        double cutoff_;              /* cutoff value for writing */
        int lastbuf_;                /* is this the last IWL buffer? 1=yes,0=no */
        int inbuf_;                  /* how many ints in current buffer? */
        int idx_;                    /* index of integral in current buffer */
        Label *labels_;              /* pointer to where integral values begin */
        Value *values_;              /* integral values */
        /*! Instance of libpsio to use */
        PSIO *psio_;
        /*! Flag indicating whether to keep the IWL file or not */
        bool keep_;

    public:

        IWL();
        IWL(PSIO *psio, int itap, double cutoff, int oldfile, int readflag);
        ~IWL();

        // Accessor functions to data
        int& itap()                         { return itap_; }
        psio_address& buffer_position()     { return bufpos_; }
        int& ints_per_buffer()              { return ints_per_buf_; }
        int& buffer_size()                  { return bufszc_; }
        double& cutoff()                    { return cutoff_; }
        int& last_buffer()                  { return lastbuf_; }
        int& buffer_count()                 { return inbuf_; }
        int& index()                        { return idx_; }
        Label* labels()                     { return labels_; }
        Value* values()                     { return values_; }
        bool& keep()                        { return keep_; }

        void init(PSIO *psio, int itap, double cutoff, int oldfile, int readflag);

        void set_keep_flag(bool k) { keep_ = k; }
        void close();

        void fetch();
        void put();

        static void read_one(PSIO *psio, int itap, const char *label, double *ints,
            int ntri, int erase, int printflg, std::string OutFileRMR);
        static void write_one(PSIO *psio, int itap, const char *label, int ntri,
            double *onel_ints);
        static void read_two(PSIO *psio, int itap, double *ints, int *ioff,
            int norbs, int nfzc, int nfzv, int printflg, std::string OutFileRMR);
        static void write_two(PSIO *psio, int itap, int nbfso, double *ints,
            int *ioff, double toler, int printflg, std::string OutFileRMR);

        int read(int target_pq, double *ints, int *ioff_lt, int *ioff_rt,
            int mp2, int printflg, std::string OutFileRMR);
        int read_all(double *ints, int *ioff_lt, int *ioff_rt, int no_pq_perm,
            int *ioff, int printflg, std::string OutFileRMR);
        int read_all2(double **ints, int *ioff_lt, int *ioff_rt, int no_pq_perm,
            int *ioff, int printflg, std::string OutFileRMR);
        int read_all_active(double *ints,
            int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
            int fstact, int lstact, int printflg,std::string OutFileRMR);
        int read_all_mp2r12a(double *ints, int *ioff_lt, int *ioff_rt,
            int bra_ket_symm, int *ioff, int printflg, std::string OutFileRMR);
        int read_array(int target_pq, double *ints, int *rlist, int *slist,
            int *size, int *ioff, int printflg, std::string OutFileRMR);
        int read_array2(double *ints, int *plist, int *qlist, int *rlist,
            int *slist, int *size, int *ioff, int printflg, std::string OutFileRMR);

        void write(int p, int q, int pq, int pqsym,
            double *arr, int rmax, int *ioff, int *orbsym,
            int *firsti, int *lasti, int printflag,
            std::string OutFileRMR);
        void write_all(int nbfso, double *ints, int *ioff, int printflg,
            std::string OutFileRMR);
        void write_mp2(int p, int q, int pq,
            int pqsym, double **arr, int rsym, int *firstr, int *lastr,
            int *firsts, int *lasts, int *occ, int *vir, int *ioff,
            int printflag, std::string OutFileRMR);
        void write_mp2r12a(int p, int q, int pq, int pqsym, double **arr,
            int rsym, int *firstr, int *lastr, int *firsts, int *lasts,
            int *occ, int bra_ket_symm, int *ioff, int printflag,
            std::string OutFileRMR);
        void write_array(double *arr, int *p, int *q, int *r, int *s, long int size);
        void write_array_SI(double *arr, short int *p, short int *q, short int *r,
            short int *s, int size);
        void write_array_SI_nocut(double *arr, short int *p, short int *q,
            short int *r, short int *s, int size);
        void write_array2(double *arr, int p, int q, int *rlist, int *slist,
            int size, int printflag, std::string OutFileRMR);
        void write_matrix(int ptr, int qtr, double **mat, int rfirst, int rlast,
            int sfirst, int slast, int *reorder, int reorder_offset,
            int printflag, int *ioff, std::string OutFileRMR);
        void write_matrix2(int ptr, int qtr, double **mat, int rfirst, int rlast,
            int sfirst, int slast, int *reorder, int reorder_offset,
            int printflag, int *ioff, std::string OutFileRMR);
        void write_value(int p, int q, int r, int s, double value, int printflag,
            std::string OutFileRMR, int dirac);
        void write_value_SI(short int p, short int q,
            short int r, short int s, double value, int printflag,
            std::string OutFileRMR, int dirac);

        void flush(int lastbuf);
        void to_end();

        static void sort_buffer(IWL *Inbuf, IWL *Outbuf,
            double *ints, int fpq, int lpq, int *ioff, int *ioff2,
            int nbfso, int elbert, int intermediate, int no_pq_perm,
            int qdim, int add, int printflg, std::string OutFileRMR);
        static void sort_buffer_pk(IWL* inbuf, int out_tape, int is_exch,
                                   double *ints, unsigned int fpq, unsigned int lpq,
                                   int *so2ind, int *so2sym, int*pksymoff, int printflg,
                                   std::string out);
    };

}

#endif
