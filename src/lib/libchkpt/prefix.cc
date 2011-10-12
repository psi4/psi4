/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

char *Chkpt::rd_prefix()
{
    char *prefix;

    //prefix = (char *) malloc(CHKPT_PREFIX_LEN*sizeof(char));
    prefix = new char[CHKPT_PREFIX_LEN];

    psio->read_entry(PSIF_CHKPT, "Default prefix", prefix, CHKPT_PREFIX_LEN*sizeof(char));

    return prefix;
}

void Chkpt::wt_prefix(const char *prefix)
{
    psio->write_entry(PSIF_CHKPT, "Default prefix", (char *) prefix, CHKPT_PREFIX_LEN*sizeof(char));
}

void Chkpt::set_prefix(const char *prefix)
{
    ::strncpy(chkpt_prefix, prefix, CHKPT_PREFIX_LEN);
    chkpt_prefix[CHKPT_PREFIX_LEN-1] = '\0';
}

void Chkpt::commit_prefix()
{
    wt_prefix(chkpt_prefix);
}

void Chkpt::reset_prefix()
{
    chkpt_prefix[0] = '\0';
}

char *Chkpt::get_prefix(void)
{
    char *prefix;

    prefix = (char *) malloc(CHKPT_PREFIX_LEN*sizeof(char));

    ::strncpy(prefix,chkpt_prefix,CHKPT_PREFIX_LEN);
    prefix[CHKPT_PREFIX_LEN-1] = '\0';

    return prefix;
}

extern "C" {
/*!
        **  char *chkpt_rd_prefix()
        **  Reads the global default chkpt prefix keyword stored in the CHKPT file.
        **
        **  returns: the prefix string
        ** \ingroup CHKPT
*/
char *chkpt_rd_prefix(void)
{
    char *prefix;
    prefix = _default_chkpt_lib_->rd_prefix();
    return prefix;
}

/*!
        **  void chkpt_wt_prefix()
        **  Writes the global default chkpt prefix keyword.
        **
        **  \param prefix = the prefix string (must be CHKPT_PREFIX_LEN long)
        **
        **  returns: none
        ** \ingroup CHKPT
*/
void chkpt_wt_prefix(const char *prefix)
{
    _default_chkpt_lib_->wt_prefix(prefix);
}


/*!
        **  void chkpt_set_prefix()
        **  Sets the default chkpt prefix in global memory.  After this is set,
        **  it is intended that all chkpt_rd_() and chkpt_wt_() calls will use
        **  this prefix for psio keyword strings.
        **
        **  \param prefix = the prefix string
        **
        **  returns: none
        ** \ingroup CHKPT
*/
void chkpt_set_prefix(const char *prefix)
{
    _default_chkpt_lib_->set_prefix(prefix);
}

/*!
        **  void chkpt_commit_prefix()
        **  Writes the default chkpt prefix from global memory into the chkpt file.
        **
        **  arguments: none
        **
        **  returns: none
        ** \ingroup CHKPT
*/
void chkpt_commit_prefix(void)
{
    _default_chkpt_lib_->commit_prefix();
}

/*!
        **  void chkpt_reset_prefix()
        **  Sets the chkpt prefix in global memory back to its default.  At
        **  present this is a null string.
        **
        **  arguments: none
        **
        **  returns: none
        ** \ingroup CHKPT
*/
void chkpt_reset_prefix(void)
{
    _default_chkpt_lib_->reset_prefix();
}

/*!
        **  char * chkpt_get_prefix()
        **  Returns a copy of the current chkpt prefix default stored
        **  in global memory.
        **
        **  arguments: none
        **
        **  returns: prefix = the current global prefix
        ** \ingroup CHKPT
*/
char *chkpt_get_prefix(void)
{
    char *prefix;
    prefix = _default_chkpt_lib_->get_prefix();
    return prefix;
}
}
