/*!
  \file snuc.c
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

int *Chkpt::rd_snuc(void)
{
	int *snuc;
	int nshell;
	char *keyword;
	keyword = build_keyword("Shell nucleus");

	nshell = rd_nshell();

	snuc = array<int>(nshell);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

	free(keyword);
	return snuc;
}

void Chkpt::wt_snuc(int *snuc)
{
	int nshell;
	char *keyword;
	keyword = build_keyword("Shell nucleus");

	nshell = rd_nshell();

	psio->write_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_snuc(): Reads in array of the nuclei numbers shells belong to.
**
**  takes no arguments.
**
**  returns: snuc = an array of the nuclei numbers to which shells 
**                  belong to.
**
** \ingroup (CHKPT)
*/
	int *chkpt_rd_snuc(void)
	{
		return _default_chkpt_lib_->rd_snuc();
	}

/*!
** chkpt_wt_snuc(): Writes out array of the nuclei numbers shells belong to.
**
**  \param snuc = an array of the nuclei numbers to which shells belong to
**
**  returns: none
**
** \ingroup (CHKPT)
*/
	void chkpt_wt_snuc(int *snuc)
	{
		_default_chkpt_lib_->wt_snuc(snuc);
	}
}
