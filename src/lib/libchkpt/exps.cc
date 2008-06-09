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

double *Chkpt::rd_exps(void)
{
	double *exps;
	int nprim = 0;
	char *keyword;
	keyword = build_keyword("Exponents");

	nprim = rd_nprim();
	exps = array<double>(nprim);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) exps, 
		nprim*sizeof(double));

	free(keyword);
	return exps;
}

void Chkpt::wt_exps(double *exps)
{
	int nprim;
	char *keyword;
	keyword = build_keyword("Exponents");

	nprim = rd_nprim();

	psio->write_entry(PSIF_CHKPT, keyword, (char *) exps, 
		nprim*sizeof(double));

	free(keyword);
}

/*!
** chkpt_rd_exps():	
** Reads in the exponents of the primitive Gaussian functions.
**
** takes no arguments.
**
** returns: double *exps   
** The exponents are returned as an array of doubles.
** \ingroup CHKPT
*/
	double *chkpt_rd_exps(void)
	{
		double *exps;
		exps = _default_chkpt_lib_->rd_exps();
		return exps;
	}

/*!
** chkpt_wt_exps(): 
** Writes out the exponents of the primitive Gaussian functions.
**
** arguments:
**  \param exps = The exponents are returned as an array of doubles.
**
** returns: none
** \ingroup CHKPT
*/
	void chkpt_wt_exps(double *exps)
	{
		_default_chkpt_lib_->wt_exps(exps);
	}
}
