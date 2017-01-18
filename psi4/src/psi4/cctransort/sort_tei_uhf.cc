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

namespace psi{ namespace cctransort {

void sort_tei_uhf(std::shared_ptr<PSIO> psio, int print)
{
  dpdbuf4 K;

  psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "KL", "I>=J+", "K>=L+", 0, "MO Ints (OO|OO)");
  global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "IJ", "KL", "A <IJ|KL>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "IJ", "KL", 0, "A <IJ|KL>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", "i>=j+", "k>=l+", 0, "MO Ints (oo|oo)");
  global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "kl", "I>=J+", "k>=l+", 0, "MO Ints (OO|oo)");
  global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "Ij", "Kl", "A <Ij|Kl>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "Ij", "Kl", 0, "A <Ij|Kl>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_AINTS, 1);

  psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "CD", "A>=B+", "C>=D+", 0, "MO Ints (VV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "AB", "CD", "B <AB|CD>");
  global_dpd_->buf4_close(&K);
  if(print > 10) {
    global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "AB", "CD", 0, "B <AB|CD>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ab", "cd", "a>=b+", "c>=d+", 0, "MO Ints (vv|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "ab", "cd", "B <ab|cd>");
  global_dpd_->buf4_close(&K);
  if(print > 10) {
    global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "cd", "A>=B+", "c>=d+", 0, "MO Ints (VV|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "Ab", "Cd", "B <Ab|Cd>");
  global_dpd_->buf4_close(&K);
  if(print > 10) {
    global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "Ab", "Cd", 0, "B <Ab|Cd>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_BINTS, 1);

  psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "AB", "I>=J+", "A>=B+", 0, "MO Ints (OO|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "IA", "JB", "C <IA|JB>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "IA", "JB", 0, "C <IA|JB>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ab", "i>=j+", "a>=b+", 0, "MO Ints (oo|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "ia", "jb", "C <ia|jb>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "ab", "I>=J+", "a>=b+", 0, "MO Ints (OO|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "Ia", "Jb", "C <Ia|Jb>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "Ia", "Jb", 0, "C <Ia|Jb>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "AB", "ij", "A>=B+", "i>=j+", 0, "MO Ints (VV|oo)");
  global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "Ai", "Bj", "C <Ai|Bj>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "Ai", "Bj", 0, "C <Ai|Bj>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_CINTS, 1);

  psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "JB", 0, "MO Ints (OV|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "IJ", "AB", "D <IJ|AB>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "IJ", "AB", 0, "D <IJ|AB>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "jb", "ia", "jb", 0, "MO Ints (ov|ov)");
  global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "ij", "ab", "D <ij|ab>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "jb", 0, "MO Ints (OV|ov)");
  global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "Ij", "Ab", "D <Ij|Ab>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "Ij", "Ab", 0, "D <Ij|Ab>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_DINTS, 1);

  psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "KA", "I>=J+", "KA", 0, "MO Ints (OO|OV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "AI", "JK", "E <AI|JK>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "AI", "JK", 0, "E <AI|JK>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ka", "i>=j+", "ka", 0, "MO Ints (oo|ov)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "ai", "jk", "E <ai|jk>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "KA", "ij", "KA", "i>=j+", 0, "MO Ints (OV|oo)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, qspr, "Ai", "Jk", "E <Ai|Jk>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "Ai", "Jk", 0, "E <Ai|Jk>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IJ", "ka", "I>=J+", "ka", 0, "MO Ints (OO|ov)");
  global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, prqs, "Ij", "Ka", "E <Ij|Ka>");
  global_dpd_->buf4_close(&K);
  if(print > 6) {
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "Ij", "Ka", 0, "E <Ij|Ka>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_EINTS, 1);

  psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "BC", "IA", "B>=C+", 0, "MO Ints (OV|VV)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "IA", "BC", "F <IA|BC>");
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "IA", "BC", 0, "F <IA|BC>");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "AI", "BC", "F <AI|BC>");
  global_dpd_->buf4_close(&K);
  if(print > 8) {
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "IA", "BC", 0, "F <IA|BC>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "bc", "ia", "b>=c+", 0, "MO Ints (ov|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "ia", "bc", "F <ia|bc>");
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "ai", "bc", "F <ai|bc>");
  global_dpd_->buf4_close(&K);
  if(print > 8) {
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "IA", "bc", "IA", "b>=c+", 0, "MO Ints (OV|vv)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "Ia", "Bc", "F <Ia|Bc>");
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ia", "Bc", 0, "F <Ia|Bc>");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "aI", "bC", "F <aI|bC>");
  global_dpd_->buf4_close(&K);
  if(print > 8) {
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ia", "Bc", 0, "F <Ia|Bc>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "aI", "bC", 0, "F <aI|bC>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }

  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "BC", "ia", "B>=C+", "ia", 0, "MO Ints (VV|ov)");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "Ai", "Bc", "F <Ai|Bc>");
  global_dpd_->buf4_close(&K);
  global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ai", "Bc", 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, rspq, "Ab", "Ci", "F <Ab|Ci>");
  global_dpd_->buf4_close(&K);
  if(print > 8) {
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ai", "Bc", 0, "F <Ai|Bc>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "Ab", "Ci", 0, "F <Cb|Ci>");
    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_close(&K);
  }
  psio->close(PSIF_CC_FINTS, 1);
}

}} // End namespaces
