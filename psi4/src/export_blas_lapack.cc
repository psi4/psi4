/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include "psi4/pybind11.h"

#include "psi4/libmints/psimath.h"

void export_blas_lapack(py::module& m)
{
    // BLAS Static Wrappers
    m.def("DGBMV", &psi::PSI_DGBMV, "docstring");
    m.def("DGEMM", &psi::PSI_DGEMM, "docstring");
    m.def("DGEMV", &psi::PSI_DGEMV, "docstring");
    m.def("DGER", &psi::PSI_DGER, "docstring");
    m.def("DSBMV", &psi::PSI_DSBMV, "docstring");
    m.def("DSYMM", &psi::PSI_DSYMM, "docstring");
    m.def("DSYMV", &psi::PSI_DSYMV, "docstring");
    m.def("DSYR", &psi::PSI_DSYR, "docstring");
    m.def("DSYR2", &psi::PSI_DSYR2, "docstring");
    m.def("DSYR2K", &psi::PSI_DSYR2K, "docstring");
    m.def("DSYRK", &psi::PSI_DSYRK, "docstring");
    m.def("DTBMV", &psi::PSI_DTBMV, "docstring");
    m.def("DTBSV", &psi::PSI_DTBSV, "docstring");
    m.def("DTRMM", &psi::PSI_DTRMM, "docstring");
    m.def("DTRMV", &psi::PSI_DTRMV, "docstring");
    m.def("DTRSM", &psi::PSI_DTRSM, "docstring");
    m.def("DTRSV", &psi::PSI_DTRSV, "docstring");
    m.def("DROT", &psi::PSI_DROT, "docstring");
    m.def("DSWAP", &psi::PSI_DSWAP, "docstring");
    m.def("DSCAL", &psi::PSI_DSCAL, "docstring");
    m.def("DAXPY", &psi::PSI_DAXPY, "docstring");
    m.def("DCOPY", &psi::PSI_DCOPY, "docstring");
    m.def("DDOT", &psi::PSI_DDOT, "docstring");
    m.def("DNRM2", &psi::PSI_DNRM2, "docstring");
    m.def("DASUM", &psi::PSI_DASUM, "docstring");
    m.def("IDAMAX", &psi::PSI_IDAMAX, "docstring");

    // LAPACK static wrappers

    m.def("DGEEV", &psi::PSI_DGEEV, "docstring");
    m.def("DSYEV", &psi::PSI_DSYEV, "docstring");
    m.def("DSYSV", &psi::PSI_DSYSV, "docstring");
    m.def("DGETRF", &psi::PSI_DGETRF, "docstring");
    m.def("DGETRS", &psi::PSI_DGETRS, "docstring");
    m.def("DGETRI", &psi::PSI_DGETRI, "docstring");
    m.def("DPOTRF", &psi::PSI_DPOTRF, "docstring");
    m.def("DPOTRS", &psi::PSI_DPOTRS, "docstring");
    m.def("DPOTRI", &psi::PSI_DPOTRI, "docstring");
}
