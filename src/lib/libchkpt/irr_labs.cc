/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

char **Chkpt::rd_irr_labs(void)
{
  int i,nirreps;
  char **irr_labs;
  psio_address ptr;
  char *keyword;
  keyword = build_keyword("Irrep labels");

  nirreps = rd_nirreps();

  ptr = PSIO_ZERO;
  irr_labs = (char **)malloc(sizeof(char *)*nirreps);
  for(i=0;i<nirreps;i++) {
    irr_labs[i] = (char *) malloc(4*sizeof(char));
    psio->read(PSIF_CHKPT, keyword, (char *) irr_labs[i],
      4*sizeof(char), ptr, &ptr);
      irr_labs[i][3] = '\0';
    }

  free(keyword);
  return irr_labs;
}

void Chkpt::wt_irr_labs(char **irr_labs)
{
  int i,nirreps;
  psio_address ptr;
  char *keyword;
  keyword = build_keyword("Irrep labels");

  nirreps = rd_nirreps();

  ptr = PSIO_ZERO;
  for(i=0;i<nirreps;i++) {
    psio->write(PSIF_CHKPT, keyword, (char *) irr_labs[i],
      4*sizeof(char), ptr, &ptr);
  }

  free(keyword);
}

char **Chkpt::rd_irr_labs_lowercase(void)
{
  int i,nirreps;
  char **irr_labs;
  psio_address ptr;
  char *keyword;
  keyword = build_keyword("Irrep labels");

  nirreps = rd_nirreps();

  ptr = PSIO_ZERO;
  irr_labs = (char **)malloc(sizeof(char *)*nirreps);
  for(i=0;i<nirreps;i++) {
    irr_labs[i] = (char *) malloc(4*sizeof(char));
    psio->read(PSIF_CHKPT, keyword, (char *) irr_labs[i],
      4*sizeof(char), ptr, &ptr);
      irr_labs[i][0] = std::tolower(irr_labs[i][0]);
      irr_labs[i][1] = std::tolower(irr_labs[i][1]);
      irr_labs[i][2] = std::tolower(irr_labs[i][2]);
      irr_labs[i][3] = '\0';
    }

  free(keyword);
  return irr_labs;
}


extern "C" {
/*!
** chkpt_rd_irr_labs(): Read in the symmetry labels for all irreps in the
**   point group in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: irr_labs =  an array of labels (strings) which denote
**      the irreps for the point group	in which the molecule is considered,
**      _regardless_ of whether there exist any symmetry orbitals which
**      transform as that irrep.
** \ingroup CHKPT
*/
        char **chkpt_rd_irr_labs(void)
        {
                return _default_chkpt_lib_->rd_irr_labs();
        }
        char **chkpt_rd_irr_labs_lowercase(void)
        {
                return _default_chkpt_lib_->rd_irr_labs_lowercase();
        }

/*!
** chkpt_wt_irr_labs(): Write out the symmetry labels for all irreps in the
**   point group in which the molecule is considered.
**
**  arguments:
**  \param irr_labs = an array of labels (strings) which denote
**      the irreps for the point group	in which the molecule is considered,
**      _regardless_ of whether there exist any symmetry orbitals which
**      transform as that irrep.
** \ingroup CHKPT
*/
        void chkpt_wt_irr_labs(char **irr_labs)
        {
                _default_chkpt_lib_->wt_irr_labs(irr_labs);
        }
}

