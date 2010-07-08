/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

// void start_io(int argc, char *argv[])
void start_io(Options& options)
{
    keep_chkpt = 0;
    read_chkpt = 0;
    chkpt_mos = 0;
    chkpt_geom = 0;
    dont_project_mos = 0;
    geomdat_geom = 0;
    save_oldcalc = 0;
    overwrite_output = 0;
    no_comshift = 0;
    no_reorient = 0;

    /*--- read MOs from checkpoint file and project onto new basis ---*/
    if (options.get_bool("CHKPT_MOS")) {
        read_chkpt = 1;
        chkpt_mos = 1;
    }
    /*--- read MOs from checkpoint file and save to a separate file ---*/
    if (options.get_bool("SAVE_MOS")) {
        read_chkpt = 1;
        save_oldcalc = 1;
    }
    /*--- don't project MOs but simply keep them ---*/
    if (options.get_bool("NOPROJECT")) {
        dont_project_mos = 1;
    }
    /*--- read geometry from checkpoint file (in findif calculations) ---*/
    if (options.get_bool("CHKPT_GEOM")) {
        read_chkpt = 0;
        chkpt_geom = 1;
        /* preserve the information about the original reference frame
         * so that properties can be rotated back later by other codes */
        keep_ref_frame = 1;
        print_lvl = 0;
        cartOn = 1;
        overwrite_output = 0;
    }
    if (options.get_bool("KEEP_OUTPUT")) {
        /* Don't overwrite the output file. This was added for the new psirb driver module. -Jet 30 Jul 07 */
        overwrite_output = 0;
    }
    if (options.get_bool("NO_COM_SHIFT")) {
        no_comshift = 1;
    }
    if (options.get_bool("NO_REORIENT")) {
        no_reorient = 1;
    }
    if (options.get_bool("KEEP_CHKPT")) {
        keep_chkpt = 1;
    }

    tstart();
}

void stop_io()
{
    tstop();
}

}} // namespace psi::input

