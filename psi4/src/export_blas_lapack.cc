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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

using namespace psi;

void export_blas_lapack(py::module& m)
{
    // BLAS Static Wrappers
    m.def("DGBMV", &PSI_DGBMV, "docstring");
    m.def("DGEMM", &PSI_DGEMM, "docstring");
    m.def("DGEMV", &PSI_DGEMV, "docstring");
    m.def("DGER", &PSI_DGER, "docstring");
    m.def("DSBMV", &PSI_DSBMV, "docstring");
    m.def("DSYMM", &PSI_DSYMM, "docstring");
    m.def("DSYMV", &PSI_DSYMV, "docstring");
    m.def("DSYR", &PSI_DSYR, "docstring");
    m.def("DSYR2", &PSI_DSYR2, "docstring");
    m.def("DSYR2K", &PSI_DSYR2K, "docstring");
    m.def("DSYRK", &PSI_DSYRK, "docstring");
    m.def("DTBMV", &PSI_DTBMV, "docstring");
    m.def("DTBSV", &PSI_DTBSV, "docstring");
    m.def("DTRMM", &PSI_DTRMM, "docstring");
    m.def("DTRMV", &PSI_DTRMV, "docstring");
    m.def("DTRSM", &PSI_DTRSM, "docstring");
    m.def("DTRSV", &PSI_DTRSV, "docstring");
    m.def("DROT", &PSI_DROT, "docstring");
    m.def("DSWAP", &PSI_DSWAP, "docstring");
    m.def("DSCAL", &PSI_DSCAL, "docstring");
    m.def("DAXPY", &PSI_DAXPY, "docstring");
    m.def("DCOPY", &PSI_DCOPY, "docstring");
    m.def("DDOT", &PSI_DDOT, "docstring");
    m.def("DNRM2", &PSI_DNRM2, "docstring");
    m.def("DASUM", &PSI_DASUM, "docstring");
    m.def("IDAMAX", &PSI_IDAMAX, "docstring");

    // LAPACK static wrappers
    m.def("DGEEV", &PSI_DGEEV, "docstring");
    m.def("DSYEV", &PSI_DSYEV, "docstring");
    m.def("DSYSV", &PSI_DSYSV, "docstring");
    m.def("DGETRF", &PSI_DGETRF, "docstring");
    m.def("DGETRS", &PSI_DGETRS, "docstring");
    m.def("DGETRI", &PSI_DGETRI, "docstring");
    m.def("DPOTRF", &PSI_DPOTRF, "docstring");
    m.def("DPOTRS", &PSI_DPOTRS, "docstring");
    m.def("DPOTRI", &PSI_DPOTRI, "docstring");
}
