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

double **Chkpt::rd_usotbf(void)
{
	double **usotbf;
	int num_so, i;
	psio_address ptr;
	char *keyword;
	keyword = build_keyword("SO->BF transmat");

	num_so = rd_nso();
	usotbf = matrix<double>(num_so,num_so);

	ptr = PSIO_ZERO;
	for(i=0;i<num_so;i++)
		psio->read(PSIF_CHKPT, keyword, (char *) usotbf[i], (int) num_so*sizeof(double), ptr, &ptr);

	free(keyword);
	return usotbf;
}

void Chkpt::wt_usotbf(double **usotbf)
{
	int num_so, i;
	psio_address ptr;
	char *keyword;
	keyword = build_keyword("SO->BF transmat");

	num_so = rd_nso();

	ptr = PSIO_ZERO;
	for(i=0;i<num_so;i++)
		psio->write(PSIF_CHKPT, keyword, (char *) usotbf[i], (int) num_so*sizeof(double), ptr, &ptr);

	free(keyword);
}

/*!
** chkpt_rd_usotbf(): Reads in the SO to basis functions transformation matrix 
**
** takes no arguments.
**
** returns: usotbf = Read in a num_so by num_so matrix of doubles
**
** \ingroup CHKPT
*/
	double **chkpt_rd_usotbf(void)
	{
		return _default_chkpt_lib_->rd_usotbf();
	}

/*!
** chkpt_wt_usotbf(): Writes out the SO to basis functions transformation 
**                    matrix 
**
** \param usotbf = A num_so by num_so matrix of doubles
**
** returns: none
**
** \ingroup CHKPT
*/
	void chkpt_wt_usotbf(double **usotbf)
	{
		_default_chkpt_lib_->wt_usotbf(usotbf);
	}
}
