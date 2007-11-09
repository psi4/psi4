/*!
  \file enuc.c
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

double Chkpt::rd_enuc(void)
{
	double enuc;
	char *keyword;
	keyword = build_keyword("Nuclear rep. energy");

	psio->read_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

	free(keyword);
	return enuc;
}

void Chkpt::wt_enuc(double enuc)
{
	char *keyword;
	keyword = build_keyword("Nuclear rep. energy");

	psio->write_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_enuc(): Reads in the nuclear repulsion energy
**
**   takes no arguments.
**
**   returns: double enuc  the nuclear repulsion energy.
**
** \ingroup (CHKPT)
*/
	double chkpt_rd_enuc(void)
	{
  		double enuc;
		enuc = _default_chkpt_lib_->rd_enuc();
  		return enuc;
	}

/*!
** chkpt_wt_enuc(): Writes out the nuclear repulsion energy
**
** \param enuc = the nuclear repulsion energy.
**
** returns: none
** \ingroup (CHKPT)
*/
	void chkpt_wt_enuc(double enuc)
	{
		_default_chkpt_lib_->wt_enuc(enuc);
	}
}
