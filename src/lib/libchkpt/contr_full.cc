/*!
  \file contr_full.cc
  \ingroup CHKPT
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

double **Chkpt::rd_contr_full(void)
{
	double **contr, *temp_contr;
	int nprim, i, j, ij = 0;
	char *keyword;
	keyword = build_keyword("Contraction coefficients");

	nprim = rd_nprim();

	temp_contr = array<double>(MAXANGMOM*nprim);
	contr = matrix<double>(nprim,MAXANGMOM);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) temp_contr,
		MAXANGMOM*nprim*sizeof(double));

/* Picking non-zero coefficients to the "master" array contr */
	for(i=0,ij=0; i < MAXANGMOM; i++) 
	for(j=0; j < nprim; j++, ij++) {
		contr[j][i] = temp_contr[ij];
	}

	free(temp_contr);
	free(keyword);
	return contr;
}

extern "C" {
/*!
** chkpt_rd_contr_full(): Reads in the normalized contraction coefficients.
**
**  takes no arguments.
**
**  returns: double **contr Normalized contraction coefficients are
**  returned as a matrix of doubles.
** \ingroup (CHKPT)
*/

	double **chkpt_rd_contr_full(void)
	{   
		return _default_chkpt_lib_->rd_contr_full();
	}   
}
