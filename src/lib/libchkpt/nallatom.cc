/*!
  \file nallatom.c
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

int Chkpt::rd_nallatom(void)
{
	int num_allatoms;
	char *keyword;
	keyword = build_keyword("Num. all atoms");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms, sizeof(int));

	free(keyword);
	return num_allatoms;
}

void Chkpt::wt_nallatom(int num_allatoms)
{
	char *keyword;
	keyword = build_keyword("Num. all atoms");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms, sizeof(int));

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_nallatom()  
**
** Reads number of all atoms (including dummy atoms)
**
** takes no arguments.
**
** returns: nallatom = number of all atoms (including dummies).
** \ingroup (CHKPT)
*/
	int chkpt_rd_nallatom(void)
	{
		return _default_chkpt_lib_->rd_nallatom();
	}


/*!
** chkpt_wt_nallatom()
**
** Writes the number of all atoms (including dummy atoms)
**
** \param nallatom  = number of all atoms (including dummies).
**
** returns: none
** \ingroup (CHKPT)
*/
	void chkpt_wt_nallatom(int num_allatoms)
	{
		_default_chkpt_lib_->wt_nallatom(num_allatoms);
	}
}
