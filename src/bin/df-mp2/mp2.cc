#include "mp2.h"
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {
namespace dfmp2 {

DFMP2::DFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    Wavefunction(options,psio,chkpt)
{
    common_init();
}
DFMP2::~DFMP2()
{
}
void DFMP2::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    reference_wavefunction_ = Process::environment.reference_wavefunction();
    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("DFMP2: Run SCF first");
    }
    copy(reference_wavefunction_);
    
    energies_["MP2J Energy"] = 0.0;
    energies_["MP2K Energy"] = 0.0;
    energies_["Reference Energy"] = reference_wavefunction_->reference_energy();

    sss_ = options_.get_double("SCALE_SS");
    oss_ = options_.get_double("SCALE_OS");

    // If the user doesn't spec a basis name, pick it yourself
    // TODO: Verify that the basis assign does not messs this up
    if (options_.get_str("RI_BASIS_MP2") == "") {
        molecule_->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "RI_BASIS_MP2");
        fprintf(outfile, "\tNo auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
        fflush(outfile);
    }

    // Form ribasis object and auxiliary basis indexing:
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    ribasis_ = BasisSet::construct(parser, molecule_, "RI_BASIS_MP2");
}
double DFMP2::compute_energy()
{
    print_header();
    form_Aia();
    form_Qia();
    form_energy();
    print_energies();
    
    return energies_["Total Energy"]; 
}
void DFMP2::print_header()
{
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *            %4s Reference, %3d Threads               *\n", options_.get_str("REFERENCE").c_str(), nthread);
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *      Rob Parrish, Justin Turney, Andy Simmonet,      *\n");
    fprintf(outfile, "\t *         Ed Hohenstein, and C. David Sheriill         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");
}
void DFMP2::print_energies()
{
    energies_["Opposite-Spin Energy"] = 0.5*energies_["MP2J Energy"];
    energies_["Same-Spin Energy"] = 0.5*energies_["MP2J Energy"] +  energies_["MP2K Energy"];
    energies_["Correlation Energy"] = energies_["MP2J Energy"] + energies_["MP2K Energy"];
    energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

    energies_["SCS Opposite-Spin Energy"] = 0.5*oss_*energies_["MP2J Energy"];
    energies_["SCS Same-Spin Energy"] = 0.5*sss_*energies_["MP2J Energy"] +  sss_*energies_["MP2K Energy"];
    energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"];
    energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2J Energy",              energies_["MP2J Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2K Energy",              energies_["MP2K Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Same-Spin Energy",         energies_["Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Opposite-Spin Energy",     energies_["Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Correlation Energy",       energies_["Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Total Energy",             energies_["Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ==================> SCS-MP2 Energies <================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Same-Spin Scale",      sss_);
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Opposite-Spin Scale",  oss_);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Same-Spin Energy",     energies_["SCS Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Opposite-Spin Energy", energies_["SCS Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Correlation Energy",   energies_["SCS Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Total Energy",         energies_["SCS Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\n");
    fflush(outfile);
}

RDFMP2::RDFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    DFMP2(options,psio,chkpt)
{
    common_init();
}
RDFMP2::~RDFMP2()
{
}
void RDFMP2::common_init()
{
    Caocc_ = Ca_subset("AO","ACTIVE_OCC");
    Cavir_ = Ca_subset("AO","ACTIVE_VIR");

    eps_aocc_ = epsilon_a_subset("AO","ACTIVE_OCC");
    eps_avir_ = epsilon_a_subset("AO","ACTIVE_VIR");
}
void RDFMP2::print_header()
{
    DFMP2::print_header();

    int focc = frzcpi_.sum();
    int fvir = frzvpi_.sum();
    int aocc = Caocc_->colspi()[0];
    int avir = Cavir_->colspi()[0];
    int occ = focc + aocc;
    int vir = fvir + avir;

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    fprintf(outfile, "\t %7s %7d %7d %7d %7d %7d %7d\n", "PAIRS", focc, occ, aocc, avir, vir, fvir);
    fprintf(outfile, "\t --------------------------------------------------------\n\n");
}
void RDFMP2::form_Aia()
{
}
void RDFMP2::form_Qia()
{
}
void RDFMP2::form_energy()
{
}

UDFMP2::UDFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    DFMP2(options,psio,chkpt)
{
    common_init();
}
UDFMP2::~UDFMP2()
{
}
void UDFMP2::common_init()
{
    Caocc_a_ = Ca_subset("AO","ACTIVE_OCC");
    Cavir_a_ = Ca_subset("AO","ACTIVE_VIR");
    Caocc_b_ = Cb_subset("AO","ACTIVE_OCC");
    Cavir_b_ = Cb_subset("AO","ACTIVE_VIR");

    eps_aocc_a_ = epsilon_a_subset("AO","ACTIVE_OCC");
    eps_avir_a_ = epsilon_a_subset("AO","ACTIVE_VIR");
    eps_aocc_b_ = epsilon_b_subset("AO","ACTIVE_OCC");
    eps_avir_b_ = epsilon_b_subset("AO","ACTIVE_VIR");
}
void UDFMP2::print_header()
{
    DFMP2::print_header();

    int focc_a = frzcpi_.sum();
    int fvir_a = frzvpi_.sum();
    int aocc_a = Caocc_a_->colspi()[0];
    int avir_a = Cavir_a_->colspi()[0];
    int occ_a = focc_a + aocc_a;
    int vir_a = fvir_a + avir_a;

    int focc_b = frzcpi_.sum();
    int fvir_b = frzvpi_.sum();
    int aocc_b = Caocc_b_->colspi()[0];
    int avir_b = Cavir_b_->colspi()[0];
    int occ_b = focc_b + aocc_b;
    int vir_b = fvir_b + avir_b;

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    fprintf(outfile, "\t %7s %7d %7d %7d %7d %7d %7d\n", "ALPHA", focc_a, occ_a, aocc_a, avir_a, vir_a, fvir_a);
    fprintf(outfile, "\t %7s %7d %7d %7d %7d %7d %7d\n", "BETA", focc_b, occ_b, aocc_b, avir_b, vir_b, fvir_b);
    fprintf(outfile, "\t --------------------------------------------------------\n\n");
}
void UDFMP2::form_Aia()
{
}
void UDFMP2::form_Qia()
{
}
void UDFMP2::form_energy()
{
}

}}

