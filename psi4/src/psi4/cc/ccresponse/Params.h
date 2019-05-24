/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#ifndef CCRESPONSE_PARAMS_H
#define CCRESPONSE_PARAMS_H

/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <string>
namespace psi {
namespace ccresponse {

struct Params {
    int print;          /* Output level control */
    long int memory;    /* Memory available (in bytes) */
    int cachelev;       /* cacheing level for libdpd */
    int ref;            /* reference determinant (0=RHF, 1=ROHF, 2=UHF) */
    double *omega;      /* energy of applied field (a.u) for dynamic polarizabilities */
    int nomega;         /* number of field energies desired */
    int maxiter;        /* maximum number of iterations allowed to converge perturbed amp eqns. */
    double convergence; /* convergence criterion for perturbed wfns */
    int restart;        /* boolean for allowing a restart from on-disk amps */
    int diis;           /* boolean for using DIIS extrapolation */
    std::string prop;   /* user-selected property */
    int local;          /* boolean for simluation of local correlation */
    int analyze;
    int dertype;
    std::string gauge; /* choice of gauge for optical rotation */
    std::string wfn;
    std::string abcd;
    int num_amps;
    int sekino; /* Sekino-Bartlett size-extensive model-III */
    int linear; /* Bartlett size-extensive (?) linear model */
};

}  // namespace ccresponse
}  // namespace psi
#endif

