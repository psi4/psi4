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

/*
** IWL.H
** Header file for Integrals With Labels Library
**
** David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
*/

#ifndef _psi_src_lib_libiwl_iwl_h_
#define _psi_src_lib_libiwl_iwl_h_

#include <cstdio>
#include "psi4/libpsio/psio.h"
#include "config.h"
#include "psi4/psi4-dec.h"
namespace psi {

struct iwlbuf {
    int itap;            /* tape number for input file */
    psio_address bufpos; /* current page/offset */
    int ints_per_buf;    /* integrals per buffer */
    int bufszc;          /* buffer size in characters (bytes) */
    double cutoff;       /* cutoff value for writing */
    int lastbuf;         /* is this the last IWL buffer? 1=yes,0=no */
    int inbuf;           /* how many ints in current buffer? */
    int idx;             /* index of integral in current buffer */
    Label *labels;       /* pointer to where integral values begin */
    Value *values;       /* integral values */
};

void PSI_API iwl_buf_fetch(struct iwlbuf *Buf);
void iwl_buf_put(struct iwlbuf *Buf);

int iwl_rdone(int itap, const char *label, double *ints, int ntri, int erase, int printflg, std::string out);

void PSI_API iwl_buf_init(struct iwlbuf *Buf, int intape, double cutoff, int oldfile, int readflag);
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf);
void PSI_API iwl_buf_close(struct iwlbuf *Buf, int keep);
void iwl_buf_wrt_val(struct iwlbuf *Buf, int p, int q, int r, int s, double value, int printflag, std::string out,
                     int dirac);
}

#endif /* end _psi_src_lib_libiwl_iwl_h */
