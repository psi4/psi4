/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
#include <libmints/mints.h>
#include "globals.h"

using namespace boost;

namespace psi {

namespace lmp2 {

PsiReturnType lmp2(Options &options) {

    boost::shared_ptr<PSIO> psio_obj(new PSIO);
//    psiopp_ipv1_config(psio_obj);
    boost::shared_ptr<Chkpt> chkpt_obj(new Chkpt(psio_obj, PSIO_OPEN_OLD));

    if (Communicator::world->me() == 0) {
        tstart();
        fprintf(
                outfile,
                "\n**************************** Begin LMP2 *********************************\n\n\n");

        fprintf(outfile, "\t\t\t*************************\n");
        fprintf(outfile, "\t\t\t*                       *\n");
        if (options.get_bool("RI_LMP2")) {
            fprintf(outfile, "\t\t\t*       DF-LMP2         *\n");
        } else {
            fprintf(outfile, "\t\t\t*         LMP2          *\n");
        }
        fprintf(outfile, "\t\t\t*                       *\n");
        fprintf(outfile, "\t\t\t*************************\n");
        fprintf(outfile, "\t\t\tRunning on %d processors\n",
                Communicator::world->nproc());
        fflush(outfile);

    }

    /** Compute the LMP2 energy **/
    {
        LMP2 lmp2_obj(psio_obj, chkpt_obj);

        if (Communicator::world->me() == 0)
            timer_on("GETPARAMS");
        lmp2_obj.get_params(options);
        if (Communicator::world->me() == 0)
            timer_off("GETPARAMS");

        if (Communicator::world->me() == 0)
            timer_on("GETMOINFO");
        lmp2_obj.get_moinfo();
        if (Communicator::world->me() == 0)
            timer_off("GETMOINFO");

        if (Communicator::world->me() == 0)
            timer_on("OPDM");
        lmp2_obj.opdm();
        if (Communicator::world->me() == 0)
            timer_off("OPDM");

        if (Communicator::world->me() == 0)
            timer_on("LOCALIZE");
        lmp2_obj.localize();
        if (Communicator::world->me() == 0)
            timer_off("LOCALIZE");

        if (Communicator::world->me() == 0)
            timer_on("GETFOCK");
        lmp2_obj.get_fock();
        if (Communicator::world->me() == 0)
            timer_off("GETFOCK");

        if (Communicator::world->me() == 0)
            timer_on("DOMAIN");
        lmp2_obj.domains();
        if (Communicator::world->me() == 0)
            timer_off("DOMAIN");

        if (Communicator::world->me() == 0)
            timer_on("PROJECT");
        lmp2_obj.projection();
        if (Communicator::world->me() == 0)
            timer_off("PROJECT");

        if (options.get_bool("RI_LMP2")) {
            lmp2_obj.direct_df_transformation3();
        } else {
            if (Communicator::world->me() == 0)
                timer_on("TRANS");
            lmp2_obj.direct_transformation();
            if (Communicator::world->me() == 0)
                timer_off("TRANS");
        }

        lmp2_obj.allocate_T();

        if (Communicator::world->me() == 0)
            timer_on("ITERATE");
        lmp2_obj.iterate();
        if (Communicator::world->me() == 0)
            timer_off("ITERATE");

    }
    /** LMP2 complete **/

    if (Communicator::world->me() == 0) {
        fprintf(
                outfile,
                "\n**************************** End of LMP2 *********************************\n");

        tstop();
    }

    return Success;
}

}
}

