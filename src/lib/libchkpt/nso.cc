/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

namespace psi {

int Chkpt::rd_nso(void)
{
	int nso;
	char *keyword;
	keyword = build_keyword("Num. SO");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &nso, sizeof(int));

	free(keyword);
	return nso;
}

void Chkpt::wt_nso(int nso)
{
	char *keyword;
	keyword = build_keyword("Num. SO");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &nso, sizeof(int));

	free(keyword);
}

/*!
** int chkpt_rd_nso()  
** Reads in the total number of SOs.
**
** returns: nso = total number of symmetry-adapted basis functions.
** \ingroup CHKPT
*/
	int chkpt_rd_nso(void)
	{
		return _default_chkpt_lib_->rd_nso();
	}

/*!
** void chkpt_wt_nso(int)  
** Writes out the total number of SOs.
**
** \param nso = total number of symmetry-adapted basis functions.
**
** returns: none
**
** \ingroup CHKPT
*/
	void chkpt_wt_nso(int nso)
	{
		_default_chkpt_lib_->wt_nso(nso);
	}
}
