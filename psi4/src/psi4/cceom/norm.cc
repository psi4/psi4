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

/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm += global_dpd_->buf4_dot_self(CMNEF);
  norm += global_dpd_->buf4_dot_self(Cmnef);
  norm += global_dpd_->buf4_dot_self(CMnEf);
  norm = sqrt(norm);
  return norm;
}

double norm_C_full(double C0, dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += C0 * C0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm += global_dpd_->buf4_dot_self(CMNEF);
  norm += global_dpd_->buf4_dot_self(Cmnef);
  norm += global_dpd_->buf4_dot_self(CMnEf);
  norm = sqrt(norm);
  return norm;
}


double dot_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm += global_dpd_->buf4_dot_self(CMNEF);
  norm += global_dpd_->buf4_dot_self(Cmnef);
  norm += global_dpd_->buf4_dot_self(CMnEf);
  return norm;
}

double dot_C_full(double C0, dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += C0 * C0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm += global_dpd_->buf4_dot_self(CMNEF);
  norm += global_dpd_->buf4_dot_self(Cmnef);
  norm += global_dpd_->buf4_dot_self(CMnEf);
  return norm;
}

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE) {
  double norm = 0.0;
  norm = 2.0 * global_dpd_->file2_dot_self(CME);
  norm += 2.0 * global_dpd_->buf4_dot_self(CMnEf);
  norm -= global_dpd_->buf4_dot(CMnEf, CMnfE);
  norm = sqrt(norm);
  return norm;
}

double norm_C_rhf_full(double C0, dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE) {
  double norm = 0.0;
  norm = C0 * C0;
  norm += 2.0 * global_dpd_->file2_dot_self(CME);
  norm += 2.0 * global_dpd_->buf4_dot_self(CMnEf);
  norm -= global_dpd_->buf4_dot(CMnEf, CMnfE);
  norm = sqrt(norm);
  return norm;
}

double norm_C1(dpdfile2 *CME, dpdfile2 *Cme)
{
  double norm = 0.0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm = sqrt(norm);
  return norm;
}

double norm_C1_full(double C0, dpdfile2 *CME, dpdfile2 *Cme)
{
  double norm = 0.0;
  norm += C0 * C0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm = sqrt(norm);
  return norm;
}

double norm_C1_rhf(dpdfile2 *CME)
{
  double norm = 0.0;
  norm = 2*global_dpd_->file2_dot_self(CME);
  norm = sqrt(norm);
  return norm;
}

double norm_C1_rhf_full(double C0, dpdfile2 *CME)
{
  double norm = 0.0;
  norm += C0 * C0;
  norm += 2*global_dpd_->file2_dot_self(CME);
  norm = sqrt(norm);
  return norm;
}


void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a)
{
  global_dpd_->file2_scm(CME,a);
  global_dpd_->file2_scm(Cme,a);
  global_dpd_->buf4_scm(CMNEF,a);
  global_dpd_->buf4_scm(Cmnef,a);
  global_dpd_->buf4_scm(CMnEf,a);
  return;
}

void scm_C_full(double *C0, dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a)
{
  (*C0) *= a;
  global_dpd_->file2_scm(CME,a);
  global_dpd_->file2_scm(Cme,a);
  global_dpd_->buf4_scm(CMNEF,a);
  global_dpd_->buf4_scm(Cmnef,a);
  global_dpd_->buf4_scm(CMnEf,a);
  return;
}

void scm_C2(dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a)
{
  global_dpd_->buf4_scm(CMNEF,a);
  global_dpd_->buf4_scm(Cmnef,a);
  global_dpd_->buf4_scm(CMnEf,a);
  return;
}

void scm_C1(dpdfile2 *CME, dpdfile2 *Cme, double a)
{
  global_dpd_->file2_scm(CME,a);
  global_dpd_->file2_scm(Cme,a);
  return;
}

void scm_C1_full(double *C0, dpdfile2 *CME, dpdfile2 *Cme, double a)
{
  (*C0) *= a; 
  global_dpd_->file2_scm(CME,a);
  global_dpd_->file2_scm(Cme,a);
  return;
}


}} // namespace psi::cceom
