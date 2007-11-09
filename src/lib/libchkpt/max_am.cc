/*!
  \file max_am.cc
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
#include <libchkpt/chkpt.h>
}
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_max_am(void)
{
	int max_am;
	char *keyword;
	keyword = build_keyword("Max. AM");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

	free(keyword);
	return max_am;
}

void Chkpt::wt_max_am(int max_am)
{
	char *keyword;
	keyword = build_keyword("Max. AM");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

	free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_max_am()
** Reads in the maximum orbital quantum number 
** of AOs in the basis.
**
** returns: max_am = (0 corresponds to s-functions, 
**                    1 - to up to p-functions, etc.)
** \ingroup (CHKPT)
*/
	int chkpt_rd_max_am(void)
	{
		int max_am;
		max_am = _default_chkpt_lib_->rd_max_am();
		return max_am;
	}

/*!
** void chkpt_wt_max_am()
** Writes out the maximum orbital quantum number 
** of AOs in the basis.
**
** arguments: 
**  \param max_am = (0 corresponds to s-functions, 
**                   1 - to up to p-functions, etc.)
**
** returns: none
** \ingroup (CHKPT)
*/
	void chkpt_wt_max_am(int max_am)
	{
		_default_chkpt_lib_->wt_max_am(max_am);
	}
}
