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

    econv_ = options_.get_double("ENERGY_CONV");
    std::cout << "e_conv = " << iconv << std::endl;

    iconv = options_.get_int("SCREENING");
    escreen_ = 1.0*pow(10, (double) -iconv);

    screen_int_ = options_.get_bool("SCREEN_INTS");

    rmsconv_ = options_.get_double("RMS_CONV");

    fs = options_.get_int("FSKIP");
    fskip_ = 1.0*pow(10.0,(double) -fs);

    diis_ = options_.get_bool("DIIS");

    diis_start_ = options_.get_int("DIIS_START_ITER");
    if(diis_start_ < 3) {
    if (me_ == 0) {
      fprintf(outfile, "\n\t*** WARNING ***\n");
      fprintf(outfile, "\tDIIS_START_ITER can not be less than 3\n");
      fprintf(outfile, "\tReseting DIIS_START_ITER to 3\n");
      fprintf(outfile, "\t***************\n");
    }
    diis_start_ = 3;
    }
    max_diis_vectors_ = options_.get_int("DIIS_MAX_VECS");

    cutoff_ = options_.get_double("LOCAL_CUTOFF");
    neglect_dp_ = options_.get_bool("NEGLECT_DP");
    dp_cutoff_ = options_.get_double("DISTANT_PAIR");

//    tol = options_.get_double("INTS_TOLERANCE");

//    scs_scale_os = options_.get_double("SCALE_OS");
//    scs_scale_ss = options_.get_double("SCALE_SS");

    only_print_domains_ = options_.get_bool("DOMAIN_PRINT");

    print_params();

}

void LMP2::print_params() const {

    if (me_ == 0) {
        fprintf(outfile, "\n  ========> LMP2 Parameters <========\n\n");
        fprintf(outfile, "  Ref WFN \t\t= %s\n", reference_.c_str());
        fprintf(outfile, "  Processes \t\t= %d\n", nproc_);
        fprintf(outfile, "  Max Iter \t\t= %d\n", maxiter_);
        fprintf(outfile, "  Energy Conv \t\t= %3.1e\n", econv_);
        fprintf(outfile, "  RMS Conv \t\t= %3.1e\n", rmsconv_);
        fprintf(outfile, "  Int Screen \t\t= %3.1e\n", escreen_);
        fprintf(outfile, "  F-Skip \t\t= %3.1e\n", fskip_);
        fprintf(outfile, "  Print \t\t= %d\n", print_);
        fprintf(outfile, "  Local Cutoff \t\t= %3.1e\n", cutoff_);
        fprintf(outfile, "  RI Approximation \t= %s\n", ri_lmp2_ ? "Yes" : "No");
        fprintf(outfile, "  Neglect Distant Pairs = %s\n", neglect_dp_ ? "Yes" : "No");
        if (neglect_dp_)
          fprintf(outfile, "  Distace Pair Cutoff \t= %2.1f\n", dp_cutoff_);
        //  fprintf(outfile,"  Opposite-spin scaled by %10.4lf\n",scs_scale_os);
        //  fprintf(outfile,"  Same-spin scaled by     %10.4lf\n",scs_scale_ss);
        fprintf(outfile, "  DIIS \t\t\t= %s\n", diis_ ? "Yes" : "No");
        if (diis_) {
            fprintf(outfile, "  DIIS Start \t\t= %d\n", diis_start_);
            fprintf(outfile, "  Max DIIS Matricies \t= %d\n", max_diis_vectors_);
        }
        fprintf(outfile, "\n  ===================================\n");
    }

}

#endif // have_madness

}}
