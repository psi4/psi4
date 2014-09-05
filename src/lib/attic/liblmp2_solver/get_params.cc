/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>


namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {

#ifdef HAVE_MADNESS

void LMP2::params() {

    int errcod, iconv, rconv, fs;
    long int max_bytes;
    char *cachetype = NULL;
    char *junk;


    /* Default glob.reference is RHF */
    reference_ = options_.get_str("REFERENCE");

    ri_lmp2_ = options_.get_bool("DF_LMP2");

    print_ = options_.get_int("PRINT");

    maxiter_ = options_.get_int("MAXITER");

    econv_ = options_.get_double("E_CONVERGENCE");
    std::cout << "e_conv = " << iconv << std::endl;

    escreen_ = options_.get_double("INTS_TOLERANCE");

    screen_int_ = options_.get_bool("SCREEN_INTS");

    rmsconv_ = options_.get_double("R_CONVERGENCE");

    fskip_ = options_.get_double("FOCK_TOLERANCE");

    diis_ = options_.get_bool("DIIS");

    diis_start_ = options_.get_int("DIIS_START_ITER");
    if(diis_start_ < 3) {
    if (me_ == 0) {
      outfile->Printf( "\n\t*** WARNING ***\n");
      outfile->Printf( "\tDIIS_START_ITER can not be less than 3\n");
      outfile->Printf( "\tReseting DIIS_START_ITER to 3\n");
      outfile->Printf( "\t***************\n");
    }
    diis_start_ = 3;
    }
    max_diis_vectors_ = options_.get_int("DIIS_MAX_VECS");

    cutoff_ = options_.get_double("LOCAL_CUTOFF");
    neglect_dp_ = options_.get_bool("NEGLECT_DISTANT_PAIR");
    dp_cutoff_ = options_.get_double("DISTANT_PAIR_CUTOFF");

//    scs_scale_os = options_.get_double("MP2_OS_SCALE");
//    scs_scale_ss = options_.get_double("MP2_SS_SCALE");

    only_print_domains_ = options_.get_bool("DOMAIN_PRINT_EXIT");

    print_params();

}

void LMP2::print_params() const {

    if (me_ == 0) {
        outfile->Printf( "\n  ========> LMP2 Parameters <========\n\n");
        outfile->Printf( "  Ref WFN \t\t= %s\n", reference_.c_str());
        outfile->Printf( "  Processes \t\t= %d\n", nproc_);
        outfile->Printf( "  Max Iter \t\t= %d\n", maxiter_);
        outfile->Printf( "  Energy Conv \t\t= %3.1e\n", econv_);
        outfile->Printf( "  RMS Conv \t\t= %3.1e\n", rmsconv_);
        outfile->Printf( "  Int Screen \t\t= %3.1e\n", escreen_);
        outfile->Printf( "  F-Skip \t\t= %3.1e\n", fskip_);
        outfile->Printf( "  Print \t\t= %d\n", print_);
        outfile->Printf( "  Local Cutoff \t\t= %3.1e\n", cutoff_);
        outfile->Printf( "  RI Approximation \t= %s\n", ri_lmp2_ ? "Yes" : "No");
        outfile->Printf( "  Neglect Distant Pairs = %s\n", neglect_dp_ ? "Yes" : "No");
        if (neglect_dp_)
          outfile->Printf( "  Distace Pair Cutoff \t= %2.1f\n", dp_cutoff_);
        //  outfile->Printf("  Opposite-spin scaled by %10.4lf\n",scs_scale_os);
        //  outfile->Printf("  Same-spin scaled by     %10.4lf\n",scs_scale_ss);
        outfile->Printf( "  DIIS \t\t\t= %s\n", diis_ ? "Yes" : "No");
        if (diis_) {
            outfile->Printf( "  DIIS Start \t\t= %d\n", diis_start_);
            outfile->Printf( "  Max DIIS Matrices \t= %d\n", max_diis_vectors_);
        }
        outfile->Printf( "\n  ===================================\n");
    }

}

#endif // have_madness

}}
