/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

namespace psi {

int *Chkpt::rd_snumg(void)
{
	int *snumg;
	int nshell;
	char *keyword;
	keyword = build_keyword("Primitives per shell");

	nshell = rd_nshell();
	snumg = array<int>(nshell);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) snumg, nshell*sizeof(int));

	free(keyword);
	return snumg;
}

void Chkpt::wt_snumg(int *snumg)
{
	int nshell;
	char *keyword;
	keyword = build_keyword("Primitives per shell");

	nshell = rd_nshell();

	psio->write_entry(PSIF_CHKPT, keyword, (char *) snumg, nshell*sizeof(int));

	free(keyword);
}

/*!
** chkpt_rd_snumg()
**
** Reads in array of the numbers of the primitive Gaussians in shells.
**
**  takes no arguments.
**
**  returns: 
**    snumg = Reads in array of the numbers of the primitive Gaussians
**            in shells
**
** \ingroup CHKPT
*/
	int *chkpt_rd_snumg(void)
	{
		return _default_chkpt_lib_->rd_snumg();
	}

/*!
** chkpt_wt_snumg()
**
** Writes out array of the numbers of the primitive Gaussians in shells.
**
**  \param snumg = array of the numbers of the primitive Gaussians
**                 in shells
**
** \ingroup CHKPT
*/
	void chkpt_wt_snumg(int *snumg)
	{
		_default_chkpt_lib_->wt_snumg(snumg);
	}
}
