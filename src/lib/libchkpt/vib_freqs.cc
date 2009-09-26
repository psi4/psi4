/*!
    \file read and write vibrational frequencies
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
  #include <libchkpt/chkpt.h>
}
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
