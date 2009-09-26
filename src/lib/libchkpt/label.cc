/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

char *Chkpt::rd_label()
{
	char *label;
	char *keyword;
	keyword = build_keyword("Label");

	label = (char *) malloc(80 * sizeof(char));

	psio->read_entry(PSIF_CHKPT, keyword, (char *) label, 80*sizeof(char));

	free(keyword);
	return label;
}

void Chkpt::wt_label(char *label)
{
	char *keyword;
	keyword = build_keyword("Label");
	
	psio->write_entry(PSIF_CHKPT, keyword, (char*)label, 80*sizeof(char));
	
	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_label():  Reads the main chkpt label.
**
**   takes no arguments.
**
**   returns: pointer to the checkpoint label
** \ingroup CHKPT
*/
	char *chkpt_rd_label(void)
	{
  		char *label;
		label = _default_chkpt_lib_->rd_label();
  		return label;
	}

/*!
** chkpt_wt_label():  Writes the main chkpt label.
**
**  arguments:
**  \param label = The calculation label.
**
**   returns: none
** \ingroup CHKPT
*/

	void chkpt_wt_label(char *label)
	{
		_default_chkpt_lib_->wt_label(label);
	}
}
