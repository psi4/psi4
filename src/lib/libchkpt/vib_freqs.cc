/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*!
    \file read and write vibrational frequencies
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double *Chkpt::rd_vib_freqs(void)
{
  int natom, nfreqs, rottype;
  double *vib_freqs = NULL;
  char *keyword;
  keyword = build_keyword("Vibrational Frequencies");

  rottype = rd_rottype();
  natom = rd_natom();

  if (rottype == 3) // linear
    nfreqs = 3*natom-5;
  else
    nfreqs = 3*natom-6;

  if (nfreqs > 0) {
    vib_freqs = array<double>(nfreqs);

    psio->read_entry(PSIF_CHKPT, keyword, (char *) vib_freqs, 
      nfreqs*sizeof(double));
  }
  free(keyword);
  return vib_freqs;
}

void Chkpt::wt_vib_freqs(double *vib_freqs)
{
  int natom, nfreqs, rottype;
  char *keyword;
  keyword = build_keyword("Vibrational Frequencies");

  rottype = rd_rottype();
  natom = rd_natom();

  if (rottype == 3) // linear
    nfreqs = 3*natom-5;
  else
    nfreqs = 3*natom-6;

  if (nfreqs > 0) {
    psio->write_entry(PSIF_CHKPT, keyword, (char *) vib_freqs,
      nfreqs*sizeof(double));
  }
  free(keyword);
  return;
}

extern "C" {
/*!
** chkpt_rd_vib_freqs()
** Reads the vibrational frequencies from the checkpoint file.
**
** arguments: none
**
** returns: 
**   double *vib_freqs: An array of the frequencies
*/
  double *chkpt_rd_vib_freqs(void)
  {
    return _default_chkpt_lib_->rd_vib_freqs();
  }

/*!
** chkpt_wt_vib_freqs()
** Writes the vibrational frequencies to the checkpoint file.
**
** \param double *vib_freqs: An array of the frequencies
**
** returns: nothing
*/
  void chkpt_wt_vib_freqs(double *vib_freqs)
  {
    _default_chkpt_lib_->wt_vib_freqs(vib_freqs);
  }
}
