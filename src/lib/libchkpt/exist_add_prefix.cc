/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
#include <libchkpt/chkpt.h>
}
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::exist_add_prefix(const char *keyword)
{
	int exists=0;
        char *keyword2;
        keyword2 = build_keyword(keyword);
	if (psio->tocscan(PSIF_CHKPT, keyword2) != NULL)
		exists=1;
        free(keyword2);		
	return exists;
}

extern "C" {
/*!
** chkpt_exist_add_prefix(): Checks to see if entry already exists in chkpt 
** file. This is like chkpt_exist() but it prepends the prefix automatically,
** so it should be ok to call by functions outside the libchkpt library.
**
**   \param keyword = keyword to look for (not including the prefix)
**  
**   returns: 1 if entry exists, 0 otherwise
**        
** \ingroup CHKPT
*/
	int chkpt_exist_add_prefix(const char *keyword)
	{
		return(_default_chkpt_lib_->exist_add_prefix(keyword));
	}
}
