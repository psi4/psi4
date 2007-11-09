/*!
  \file sloc.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
	#include <libchkpt/chkpt.h>
}
#include <libchkpt/chkpt.hpp>

using namespace psi;

int *Chkpt::rd_sloc(void)
{
	int *sloc;
	int nshell;
	char *keyword;
	keyword = build_keyword("First AO per shell");

	nshell = rd_nshell();

	sloc = array<int>(nshell);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

	free(keyword);
	return sloc;
}

void Chkpt::wt_sloc(int *sloc)
{
	int nshell;
	char *keyword;
	keyword = build_keyword("First AO per shell");

	nshell = rd_nshell();

	psio->write_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_sloc():	Read in an array of the numbers of the first AO 
**			from the shells.
**
**  takes no arguments.
**
**  returns: sloc = An array nshell long of the numbers of 
**                  the first AOs from the shells.
**
** \ingroup (CHKPT)
*/
	int *chkpt_rd_sloc(void)
	{
		return _default_chkpt_lib_->rd_sloc();
	}

/*!
** chkpt_wt_sloc():	
**  Writes out an array of the numbers of the first AO from the shells.
**
**  \param sloc = An array nshell long of the numbers of the first AOs 
**                from the shells.
**  returns: none
**
** \ingroup (CHKPT)
*/
	void chkpt_wt_sloc(int *sloc)
	{
		_default_chkpt_lib_->wt_sloc(sloc);
	}
}
