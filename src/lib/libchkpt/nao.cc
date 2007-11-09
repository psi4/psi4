/*!
  \file nao.cc
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

int Chkpt::rd_nao(void)
{
	int nao;
	char *keyword;
	keyword = build_keyword("Num. AO");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &nao, sizeof(int));

	free(keyword);
	return nao;
}

void Chkpt::wt_nao(int nao)
{
	char *keyword;
	keyword = build_keyword("Num. AO");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &nao, sizeof(int));

	free(keyword);
}

extern "C" {
/*!
**  int chkpt_rd_nao()
**  Reads in the total number of atomic orbitals.
**
**  returns: nao = total number of atomic orbitals.
**  \ingroup (CHKPT)
*/
	int chkpt_rd_nao(void)
	{
		return _default_chkpt_lib_->rd_nao();
	}

/*!
**  void chkpt_wt_nao(int)
**  Writes out the total number of atomic orbitals.
**
**   \param nao = total number of atomic orbitals.
**
**  returns: none
**  \ingroup (CHKPT)
*/
	void chkpt_wt_nao(int nao)
	{
		_default_chkpt_lib_->wt_nao(nao);
	}
}
