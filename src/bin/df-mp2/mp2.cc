#include "mp2.h"
#include <libmints/mints.h>
#include <libmints/sieve.h>
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

    if (options_.get_str("RI_BASIS_MP2") == "") {
        molecule_->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "RI_BASIS_MP2");
        fprintf(outfile, "\tNo auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
        fflush(outfile);
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    ribasis_ = BasisSet::construct(parser, molecule_, "RI_BASIS_MP2");
}
double DFMP2::compute_energy()
{
    print_header();
    form_singles();
    form_Aia();
    form_Qia();
    form_energy();
    print_energies();
    
    return energies_["Total Energy"]; 
}
void DFMP2::form_singles()
{
    double E_singles_a = 0.0;    
    double E_singles_b = 0.0;    

    SharedMatrix Caocc_a = Ca_subset("SO","ACTIVE_OCC");
    SharedMatrix Caocc_b = Cb_subset("SO","ACTIVE_OCC");
    SharedMatrix Cavir_a = Ca_subset("SO","ACTIVE_VIR");
    SharedMatrix Cavir_b = Cb_subset("SO","ACTIVE_VIR");

    SharedVector eps_aocc_a = epsilon_a_subset("SO","ACTIVE_OCC");
    SharedVector eps_aocc_b = epsilon_b_subset("SO","ACTIVE_OCC");
    SharedVector eps_avir_a = epsilon_a_subset("SO","ACTIVE_VIR");
    SharedVector eps_avir_b = epsilon_b_subset("SO","ACTIVE_VIR");

    SharedMatrix Fia_a(new Matrix("Fia a", Caocc_a->colspi(), Cavir_a->colspi()));  
    SharedMatrix Fia_b(new Matrix("Fia b", Caocc_b->colspi(), Cavir_b->colspi()));  

    double* temp = new double[Fa_->max_nrow() * (ULI) (Cavir_a->max_ncol() > Cavir_b->max_ncol() ?
        Cavir_a->max_ncol() : Cavir_b->max_ncol())];

    // Fia a
    for (int h = 0; h < Caocc_a->nirrep(); h++) {
        
        int nso = Fa_->rowspi()[h];
        int naocc = Caocc_a->colspi()[h];
        int navir = Cavir_a->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        double** Fsop = Fa_->pointer(h);
        double** Fmop = Fia_a->pointer(h);
        double** Cip = Caocc_a->pointer(h); 
        double** Cap = Cavir_a->pointer(h); 

        C_DGEMM('N','N',nso,navir,nso,1.0,Fsop[0],nso,Cap[0],navir,0.0,temp,navir);
        C_DGEMM('T','N',naocc,navir,nso,1.0,Cip[0],naocc,temp,navir,0.0,Fmop[0],navir); 

        double* eps_i = eps_aocc_a->pointer(h);
        double* eps_a = eps_avir_a->pointer(h);

        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_a += Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }  
 
    // Fia b
    for (int h = 0; h < Caocc_b->nirrep(); h++) {
        
        int nso = Fb_->rowspi()[h];
        int naocc = Caocc_b->colspi()[h];
        int navir = Cavir_b->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        double** Fsop = Fb_->pointer(h);
        double** Fmop = Fia_b->pointer(h);
        double** Cip = Caocc_b->pointer(h); 
        double** Cap = Cavir_b->pointer(h); 

        double* eps_i = eps_aocc_b->pointer(h);
        double* eps_a = eps_avir_b->pointer(h);

        C_DGEMM('N','N',nso,navir,nso,1.0,Fsop[0],nso,Cap[0],navir,0.0,temp,navir);
        C_DGEMM('T','N',naocc,navir,nso,1.0,Cip[0],naocc,temp,navir,0.0,Fmop[0],navir); 
        
        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_b += Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }  
 
    delete[] temp;

    energies_["Singles Energy"] = E_singles_a + E_singles_b;

    if (debug_) {
        Fia_a->print(); 
        Fia_b->print(); 
        fprintf(outfile, "  Alpha singles energy = %24.16E\n", E_singles_a);
        fprintf(outfile, "  Beta  singles energy = %24.16E\n\n", E_singles_b);
    }
}
void DFMP2::print_energies()
{
    energies_["Opposite-Spin Energy"] = 0.5*energies_["MP2J Energy"];
    energies_["Same-Spin Energy"] = 0.5*energies_["MP2J Energy"] +  energies_["MP2K Energy"];
    energies_["Correlation Energy"] = energies_["MP2J Energy"] + energies_["MP2K Energy"] + energies_["Singles energy"];
    energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

    energies_["SCS Opposite-Spin Energy"] = 0.5*oss_*energies_["MP2J Energy"];
    energies_["SCS Same-Spin Energy"] = 0.5*sss_*energies_["MP2J Energy"] +  sss_*energies_["MP2K Energy"];
    energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"] + energies_["Singles energy"];
    energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Singles Energy",           energies_["Singles Energy"]);
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
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *            RMP2 Wavefunction, %3d Threads            *\n", nthread);
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *      Rob Parrish, Justin Turney, Andy Simmonet,      *\n");
    fprintf(outfile, "\t *         Ed Hohenstein, and C. David Sheriill         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

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
    /**
    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_,options_.get_double("SCHWARZ_CUTOFF")))

    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    int maxQ = ribasis_->max_function_per_shell();

    ULI Amn_memory = maxQ * (ULI) nso * nso;
    ULI Ami_memory = maxQ * (ULI) nso * naocc;
    
    ULI Aia_memory = (ULI) (0.8 * (memory_ / 8L)) - Amn_memory - Ami_memory;
    ULI max_temp = Aia_memory / (naocc * (ULI) navir);

    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    SharedMatrix Amn(new Matrix("(A|mn) Block", maxQ, nso * (ULI) nso));
    SharedMatrix Ami(new Matrix("(A|mi) Block", maxQ, nso * (ULI) naocc));
    SharedMatrix Aia(new Matrix("(A|ia) Block", max_naux, naocc * (ULI) navir));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;

    int current_naux = 0;
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
    
        int Qstart = ribasis_->shell(Q)->function_index();
        int nQ = ribasis_->shell(Q)->nfunction();


        // Flush the buffer
        if (current_naux + nQ > max_naux) {
            psio_->write(PSIF_DFMP2_AIA,"Aia",(char*) Aiap[0],sizeof(double)*current_naux*naocc*navir,next_AIA,&next_AIA);
            current_naux = 0;
        }

        // Put integrals in the buffer
        current_naux += nQ;
    
    }

    psio_->close(PSIF_DFMP2_AIA,1);
    **/
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
    // Bit of a hack: if in ROHF-MBPT(2), semicanonicalize before the UDFMP2 bit picks up C or epsilon
    if (options_.get_str("REFERENCE") == "ROHF") 
        semicanonicalize();    

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
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *            UMP2 Wavefunction, %3d Threads            *\n", nthread);
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *      Rob Parrish, Justin Turney, Andy Simmonet,      *\n");
    fprintf(outfile, "\t *         Ed Hohenstein, and C. David Sheriill         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

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

RODFMP2::RODFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    UDFMP2(options,psio,chkpt)
{
    common_init();
}
RODFMP2::~RODFMP2()
{
}
void RODFMP2::common_init()
{
}
void RODFMP2::print_header()
{
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *        ROHF-MBPT(2) Wavefunction, %3d Threads        *\n", nthread);
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *      Rob Parrish, Justin Turney, Andy Simmonet,      *\n");
    fprintf(outfile, "\t *         Ed Hohenstein, and C. David Sheriill         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

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

}}

