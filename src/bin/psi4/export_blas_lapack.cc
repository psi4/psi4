/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <boost/python.hpp>
#include <libmints/psimath.h>

using namespace boost::python;

void export_blas_lapack()
{
    // BLAS Static Wrappers
    def("DGBMV", &psi::PSI_DGBMV, "docstring");
    def("DGEMM", &psi::PSI_DGEMM, "docstring");
    def("DGEMV", &psi::PSI_DGEMV, "docstring");
    def("DGER", &psi::PSI_DGER, "docstring");
    def("DSBMV", &psi::PSI_DSBMV, "docstring");
    def("DSYMM", &psi::PSI_DSYMM, "docstring");
    def("DSYMV", &psi::PSI_DSYMV, "docstring");
    def("DSYR", &psi::PSI_DSYR, "docstring");
    def("DSYR2", &psi::PSI_DSYR2, "docstring");
    def("DSYR2K", &psi::PSI_DSYR2K, "docstring");
    def("DSYRK", &psi::PSI_DSYRK, "docstring");
    def("DTBMV", &psi::PSI_DTBMV, "docstring");
    def("DTBSV", &psi::PSI_DTBSV, "docstring");
    def("DTRMM", &psi::PSI_DTRMM, "docstring");
    def("DTRMV", &psi::PSI_DTRMV, "docstring");
    def("DTRSM", &psi::PSI_DTRSM, "docstring");
    def("DTRSV", &psi::PSI_DTRSV, "docstring");
    def("DROT", &psi::PSI_DROT, "docstring");
    def("DSWAP", &psi::PSI_DSWAP, "docstring");
    def("DSCAL", &psi::PSI_DSCAL, "docstring");
    def("DAXPY", &psi::PSI_DAXPY, "docstring");
    def("DCOPY", &psi::PSI_DCOPY, "docstring");
    def("DDOT", &psi::PSI_DDOT, "docstring");
    def("DNRM2", &psi::PSI_DNRM2, "docstring");
    def("DASUM", &psi::PSI_DASUM, "docstring");
    def("IDAMAX", &psi::PSI_IDAMAX, "docstring");

    // LAPACK static wrappers

    def("DGEEV", &psi::PSI_DGEEV, "docstring");
    def("DSYEV", &psi::PSI_DSYEV, "docstring");
    def("DSYSV", &psi::PSI_DSYSV, "docstring");
    def("DGETRF", &psi::PSI_DGETRF, "docstring");
    def("DGETRS", &psi::PSI_DGETRS, "docstring");
    def("DGETRI", &psi::PSI_DGETRI, "docstring");
    def("DPOTRF", &psi::PSI_DPOTRF, "docstring");
    def("DPOTRS", &psi::PSI_DPOTRS, "docstring");
    def("DPOTRI", &psi::PSI_DPOTRI, "docstring");
    /**
    def("DBDSDC", &psi::PSI_DBDSDC);
    def("DBDSQR", &psi::PSI_DBDSQR);
    def("DDISNA", &psi::PSI_DDISNA);
    def("DGBBRD", &psi::PSI_DGBBRD);
    def("DGBCON", &psi::PSI_DGBCON);
    def("DGBEQU", &psi::PSI_DGBEQU);
    def("DGBRFS", &psi::PSI_DGBRFS);
    def("DGBSV", &psi::PSI_DGBSV);
    def("DGBSVX", &psi::PSI_DGBSVX);
    def("DGBTRF", &psi::PSI_DGBTRF);
    def("DGBTRS", &psi::PSI_DGBTRS);
    def("DGEBAK", &psi::PSI_DGEBAK);
    def("DGEBAL", &psi::PSI_DGEBAL);
    def("DGEBRD", &psi::PSI_DGEBRD);
    def("DGECON", &psi::PSI_DGECON);
    def("DGEEQU", &psi::PSI_DGEEQU);
    def("DGEES", &psi::PSI_DGEES);
    def("DGEESX", &psi::PSI_DGEESX);
    def("DGEEVX", &psi::PSI_DGEEVX);
    def("DGEGS", &psi::PSI_DGEGS);
    def("DGEGV", &psi::PSI_DGEGV);
    def("DGEHRD", &psi::PSI_DGEHRD);
    def("DGELQF", &psi::PSI_DGELQF);
    def("DGELS", &psi::PSI_DGELS);
    def("DGELSD", &psi::PSI_DGELSD);
    def("DGELSS", &psi::PSI_DGELSS);
    def("DGELSX", &psi::PSI_DGELSX);
    def("DGELSY", &psi::PSI_DGELSY);
    def("DGEQLF", &psi::PSI_DGEQLF);
    def("DGEQP3", &psi::PSI_DGEQP3);
    def("DGEQPF", &psi::PSI_DGEQPF);
    def("DGERFS", &psi::PSI_DGERFS);
    def("DGERQF", &psi::PSI_DGERQF);
    def("DGESDD", &psi::PSI_DGESDD);
    def("DGESV", &psi::PSI_DGESV);
    def("DGESVX", &psi::PSI_DGESVX);
    def("DGETRF", &psi::PSI_DGETRF);
    def("DGETRI", &psi::PSI_DGETRI);
    def("DGETRS", &psi::PSI_DGETRS);
    def("DGGBAK", &psi::PSI_DGGBAK);
    def("DGGBAL", &psi::PSI_DGGBAL);
    def("DGGES", &psi::PSI_DGGES);
    def("DGGESX", &psi::PSI_DGGESX);
    def("DGGEV", &psi::PSI_DGGEV);
    def("DGGEVX", &psi::PSI_DGGEVX);
    def("DGGGLM", &psi::PSI_DGGGLM);
    def("DGGHRD", &psi::PSI_DGGHRD);
    def("DGGLSE", &psi::PSI_DGGLSE);
    def("DGGQRF", &psi::PSI_DGGQRF);
    def("DGGRQF", &psi::PSI_DGGRQF);
    def("DGGSVD", &psi::PSI_DGGSVD);
    def("DGGSVP", &psi::PSI_DGGSVP);
    def("DGTCON", &psi::PSI_DGTCON);
    def("DGTRFS", &psi::PSI_DGTRFS);
    def("DGTSV", &psi::PSI_DGTSV);
    def("DGTSVX", &psi::PSI_DGTSVX);
    def("DGTTRF", &psi::PSI_DGTTRF);
    def("DGTTRS", &psi::PSI_DGTTRS);
    def("DHGEQZ", &psi::PSI_DHGEQZ);
    def("DHSEIN", &psi::PSI_DHSEIN);
    def("DHSEQR", &psi::PSI_DHSEQR);
    def("DORGBR", &psi::PSI_DORGBR);
    def("DORGHR", &psi::PSI_DORGHR);
    def("DORGLQ", &psi::PSI_DORGLQ);
    def("DORGQL", &psi::PSI_DORGQL);
    def("DORGQR", &psi::PSI_DORGQR);
    def("DORGRQ", &psi::PSI_DORGRQ);
    def("DORGTR", &psi::PSI_DORGTR);
    def("DORMBR", &psi::PSI_DORMBR);
    def("DORMHR", &psi::PSI_DORMHR);
    def("DORMLQ", &psi::PSI_DORMLQ);
    def("DORMQL", &psi::PSI_DORMQL);
    def("DORMQR", &psi::PSI_DORMQR);
    def("DORMR3", &psi::PSI_DORMR3);
    def("DORMRQ", &psi::PSI_DORMRQ);
    def("DORMRZ", &psi::PSI_DORMRZ);
    def("DORMTR", &psi::PSI_DORMTR);
    def("DPBCON", &psi::PSI_DPBCON);
    def("DPBEQU", &psi::PSI_DPBEQU);
    def("DPBRFS", &psi::PSI_DPBRFS);
    def("DPBSTF", &psi::PSI_DPBSTF);
    def("DPBSV", &psi::PSI_DPBSV);
    def("DPBSVX", &psi::PSI_DPBSVX);
    def("DPBTRF", &psi::PSI_DPBTRF);
    def("DPBTRS", &psi::PSI_DPBTRS);
    def("DPOCON", &psi::PSI_DPOCON);
    def("DPOEQU", &psi::PSI_DPOEQU);
    def("DPORFS", &psi::PSI_DPORFS);
    def("DPOSV", &psi::PSI_DPOSV);
    def("DPOSVX", &psi::PSI_DPOSVX);
    def("DPOTRF", &psi::PSI_DPOTRF);
    def("DPOTRI", &psi::PSI_DPOTRI);
    def("DPOTRS", &psi::PSI_DPOTRS);
    def("DPTCON", &psi::PSI_DPTCON);
    def("DPTEQR", &psi::PSI_DPTEQR);
    def("DPTRFS", &psi::PSI_DPTRFS);
    def("DPTSV", &psi::PSI_DPTSV);
    def("DPTSVX", &psi::PSI_DPTSVX);
    def("DPTTRF", &psi::PSI_DPTTRF);
    def("DPTTRS", &psi::PSI_DPTTRS);
    def("DSBEV", &psi::PSI_DSBEV);
    def("DSBEVD", &psi::PSI_DSBEVD);
    def("DSBEVX", &psi::PSI_DSBEVX);
    def("DSBGST", &psi::PSI_DSBGST);
    def("DSBGV", &psi::PSI_DSBGV);
    def("DSBGVD", &psi::PSI_DSBGVD);
    def("DSBGVX", &psi::PSI_DSBGVX);
    def("DSBTRD", &psi::PSI_DSBTRD);
    def("DSGESV", &psi::PSI_DSGESV);
    def("DSTEBZ", &psi::PSI_DSTEBZ);
    def("DSTEDC", &psi::PSI_DSTEDC);
    def("DSTEGR", &psi::PSI_DSTEGR);
    def("DSTEIN", &psi::PSI_DSTEIN);
    def("DSTEQR", &psi::PSI_DSTEQR);
    def("DSTERF", &psi::PSI_DSTERF);
    def("DSTEV", &psi::PSI_DSTEV);
    def("DSTEVD", &psi::PSI_DSTEVD);
    def("DSTEVR", &psi::PSI_DSTEVR);
    def("DSTEVX", &psi::PSI_DSTEVX);
    def("DSYCON", &psi::PSI_DSYCON);
    def("DSYGST", &psi::PSI_DSYGST);
    def("DSYGV", &psi::PSI_DSYGV);
    def("DSYGVD", &psi::PSI_DSYGVD);
    def("DSYGVX", &psi::PSI_DSYGVX);
    def("DSYRFS", &psi::PSI_DSYRFS);
    def("DSYTRD", &psi::PSI_DSYTRD);
    def("DSYTRF", &psi::PSI_DSYTRF);
    def("DSYTRI", &psi::PSI_DSYTRI);
    def("DSYTRS", &psi::PSI_DSYTRS);
    def("DTBCON", &psi::PSI_DTBCON);
    def("DTBRFS", &psi::PSI_DTBRFS);
    def("DTBTRS", &psi::PSI_DTBTRS);
    def("DTGEVC", &psi::PSI_DTGEVC);
    def("DTGEXC", &psi::PSI_DTGEXC);
    def("DTGSEN", &psi::PSI_DTGSEN);
    def("DTGSJA", &psi::PSI_DTGSJA);
    def("DTGSNA", &psi::PSI_DTGSNA);
    def("DTGSYL", &psi::PSI_DTGSYL);
    def("DTRCON", &psi::PSI_DTRCON);
    def("DTREVC", &psi::PSI_DTREVC);
    def("DTREXC", &psi::PSI_DTREXC);
    def("DTRRFS", &psi::PSI_DTRRFS);
    def("DTRSEN", &psi::PSI_DTRSEN);
    def("DTRSNA", &psi::PSI_DTRSNA);
    def("DTRSYL", &psi::PSI_DTRSYL);
    def("DTRTRI", &psi::PSI_DTRTRI);
    def("DTRTRS", &psi::PSI_DTRTRS);
    def("DTZRQF", &psi::PSI_DTZRQF);
    def("DTZRZF", &psi::PSI_DTZRZF);
    **/
}