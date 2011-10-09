/*!
  \file read and write the rotational symmetry number
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_rot_symm_num(void)
{
  int rot_symm_num;
  char *keyword;
  keyword = build_keyword("Rotational Symmetry Number");

  psio->read_entry(PSIF_CHKPT, keyword, (char *) &rot_symm_num, sizeof(int));

  free(keyword);
  return rot_symm_num;
}

void Chkpt::wt_rot_symm_num(int rot_symm_num)
{
  char *keyword;
  keyword = build_keyword("Rotational Symmetry Number");

  psio->write_entry(PSIF_CHKPT, keyword, (char *) &rot_symm_num, sizeof(int));

  free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_rot_symm_num()
** Reads the rotational symmetry number.
**
** returns: rot_symm_num = rotational symmetry number
** \ingroup CHKPT
*/
  int chkpt_rd_rot_symm_num(void)
  {
  return _default_chkpt_lib_->rd_rot_symm_num();
  }

/*!
** void chkpt_wt_rot_symm_num(int)
** Writes the rotational symmetry number.
**
** \param rot_symm_num = rotational symmetry number
**
** \ingroup CHKPT
*/
  void chkpt_wt_rot_symm_num(int rot_symm_num)
  {
  _default_chkpt_lib_->wt_rot_symm_num(rot_symm_num);
  }
}
