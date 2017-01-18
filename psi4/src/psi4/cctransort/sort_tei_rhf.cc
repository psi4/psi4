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

#include "psi4/libpsio/psio.hpp"
#include "psi4/libdpd/dpd.h"

namespace psi { namespace cctransort {

void sort_tei_rhf(std::shared_ptr<PSIO> psio, int print)
{
  dpdbuf4 K;

  psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", "i>=j+", "k>=l+", 0, "MO Ints (OO|OO)");
  global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_AINTS, 1);

  psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ab", "cd", "a>=b+", "c>=d+", 0, "MO Ints (VV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "ab", "cd", "B <ab|cd>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_BINTS, 1);

  psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ab", "i>=j+", "a>=b+", 0, "MO Ints (OO|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "ia", "jb", "C <ia|jb>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_CINTS, 1);

  psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "jb", "ia", "jb", 0, "MO Ints (OV|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "ij", "ab", "D <ij|ab>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_DINTS, 1);

  psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ka", "i>=j+", "ka", 0, "MO Ints (OO|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "ai", "jk", "E <ai|jk>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_EINTS, 1);

  psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "bc", "ia", "b>=c+", 0, "MO Ints (OV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "ia", "bc", "F <ia|bc>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "ai", "bc", "F <ai|bc>");
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_CC_FINTS, 1);
}

}} // End namespaces
