/*!
  \file num_unique_shell.cc
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

int Chkpt::rd_num_unique_shell(void)
{
	int nunique;
	char *keyword;
	keyword = build_keyword("Num. unique shells");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

	free(keyword);
	return nunique;
}

void Chkpt::wt_num_unique_shell(int nunique)
{
	char *keyword;
	keyword = build_keyword("Num. unique shells");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

	free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_num_unique_shell()  
** Reads in the number of symmetry unique shells. 
**
** returns: nunique = number of symmetry unique shells.
** \ingroup (CHKPT)
*/
	int chkpt_rd_num_unique_shell(void)
	{
		return _default_chkpt_lib_->rd_num_unique_shell();
	}

/*!
** void chkpt_wt_num_unique_shell(int)  
** Writes out the number of symmetry unique shells. 
**
** \param nunique = number of symmetry unique shells.
**
** returns: none
** \ingroup (CHKPT)
*/
	void chkpt_wt_num_unique_shell(int nunique)
	{
		_default_chkpt_lib_->wt_num_unique_shell(nunique);
	}
}
