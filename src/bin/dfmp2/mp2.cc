#include "mp2.h"
#include "corr_grad.h"
#include <lib3index/3index.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libfock/jk.h>
#include <libfock/apps.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <physconst.h>
#include <psifiles.h>

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

    if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
        reference_wavefunction_->semicanonicalize();

    copy(reference_wavefunction_);

    energies_["Singles Energy"] = 0.0;
    energies_["Opposite-Spin Energy"] = 0.0;
    energies_["Same-Spin Energy"] = 0.0;
    energies_["Reference Energy"] = reference_wavefunction_->reference_energy();

    sss_ = options_.get_double("MP2_SS_SCALE");
    oss_ = options_.get_double("MP2_OS_SCALE");

    if (options_.get_str("DF_BASIS_MP2") == "") {
        molecule_->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "DF_BASIS_MP2");
        fprintf(outfile, "\tNo auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
        fflush(outfile);
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    ribasis_ = BasisSet::construct(parser, molecule_, "DF_BASIS_MP2");
}
double DFMP2::compute_energy()
{
    print_header();
    timer_on("DFMP2 Singles");
    form_singles();
    timer_off("DFMP2 Singles");
    timer_on("DFMP2 Aia");
    form_Aia();
    timer_off("DFMP2 Aia");
    timer_on("DFMP2 Qia");
    form_Qia();
    timer_off("DFMP2 Qia");
    timer_on("DFMP2 Energy");
    form_energy();
    timer_off("DFMP2 Energy");
    print_energies();

    return energies_["Total Energy"];
}
SharedMatrix DFMP2::compute_gradient()
{
    print_header();

    timer_on("DFMP2 Singles");
    form_singles();
    timer_off("DFMP2 Singles");

    timer_on("DFMP2 Aia");
    form_Aia();
    timer_off("DFMP2 Aia");

    timer_on("DFMP2 iaQ");
    form_Qia_gradient();
    timer_off("DFMP2 iaQ");

    timer_on("DFMP2 aiQ");
    form_Qia_transpose();
    timer_off("DFMP2 aiQ");

    timer_on("DFMP2 Tij");
    form_Pab();
    timer_off("DFMP2 Tij");

    timer_on("DFMP2 Tab");
    form_Pij();
    timer_off("DFMP2 Tab");

    timer_on("DFMP2 gamma");
    form_gamma();
    timer_off("DFMP2 gamma");

    timer_on("DFMP2 Gia");
    form_G_transpose();
    timer_off("DFMP2 Gia");

    timer_on("DFMP2 AB^x");
    form_AB_x_terms();
    timer_off("DFMP2 AB^x");

    timer_on("DFMP2 Amn^x");
    form_Amn_x_terms();
    timer_off("DFMP2 Amn^x");

    timer_on("DFMP2 L");
    form_L();
    timer_off("DFMP2 L");

    timer_on("DFMP2 P");
    form_P();
    timer_off("DFMP2 P");

    timer_on("DFMP2 W");
    form_W();
    timer_off("DFMP2 W");

    timer_on("DFMP2 Z");
    form_Z();
    timer_off("DFMP2 Z");

    timer_on("DFMP2 grad");
    form_gradient();
    timer_off("DFMP2 grad");

    print_energies();

    print_gradients();

    return gradients_["Total"]; 
}
void DFMP2::form_singles()
{
    double E_singles_a = 0.0;
    double E_singles_b = 0.0;

    SharedMatrix Caocc_a = Ca_subset("SO","ACTIVE_OCC");
    SharedMatrix Cavir_a = Ca_subset("SO","ACTIVE_VIR");
    SharedMatrix Caocc_b = Cb_subset("SO","ACTIVE_OCC");
    SharedMatrix Cavir_b = Cb_subset("SO","ACTIVE_VIR");

    SharedVector eps_aocc_a = epsilon_a_subset("SO","ACTIVE_OCC");
    SharedVector eps_avir_a = epsilon_a_subset("SO","ACTIVE_VIR");
    SharedVector eps_aocc_b = epsilon_b_subset("SO","ACTIVE_OCC");
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
                E_singles_a -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
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
                E_singles_b -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }

    delete[] temp;

    energies_["Singles Energy"] = E_singles_a + E_singles_b;

    if (debug_) {
        Caocc_a->print();
        Cavir_a->print();
        eps_aocc_a->print();
        eps_avir_a->print();
        Caocc_b->print();
        Cavir_b->print();
        eps_aocc_b->print();
        eps_avir_b->print();

        Fia_a->print();
        Fia_b->print();
        fprintf(outfile, "  Alpha singles energy = %24.16E\n", E_singles_a);
        fprintf(outfile, "  Beta  singles energy = %24.16E\n\n", E_singles_b);
    }
}
SharedMatrix DFMP2::form_inverse_metric()
{
    timer_on("DFMP2 Metric");

    int naux = ribasis_->nbf();

    // Load inverse metric from the SCF three-index integral file if it exists
    if (options_.get_str("DF_INTS_IO") == "LOAD") {

        SharedMatrix Jm12(new Matrix("SO Basis Fitting Inverse (Eig)", naux, naux));
        fprintf(outfile,"\t Will attempt to load fitting metric from file %d.\n\n", PSIF_DFSCF_BJ); fflush(outfile);
        psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*) Jm12->pointer()[0], sizeof(double) * naux * naux);
        psio_->close(PSIF_DFSCF_BJ, 1);

        timer_off("DFMP2 Metric");

        return Jm12;

    } else {

        // Form the inverse metric manually
        boost::shared_ptr<FittingMetric> metric(new FittingMetric(ribasis_, true));
        metric->form_eig_inverse(1.0E-10);
        SharedMatrix Jm12 = metric->get_metric();

        // Save inverse metric to the SCF three-index integral file if it exists
        if (options_.get_str("DF_INTS_IO") == "SAVE") {
            fprintf(outfile,"\t Will save fitting metric to file %d.\n\n", PSIF_DFSCF_BJ); fflush(outfile);
            psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
            psio_->write_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*) Jm12->pointer()[0], sizeof(double) * naux * naux);
            psio_->close(PSIF_DFSCF_BJ, 1);
        }

        timer_off("DFMP2 Metric");

        return Jm12;
    }
}
void DFMP2::apply_fitting(SharedMatrix Jm12, unsigned int file, ULI naux, ULI nia)
{
    // Memory constraints
    ULI Jmem = naux * naux;
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
    }
    ULI rem = (doubles - Jmem) / 2L;
    ULI max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<ULI> ia_starts;
    ia_starts.push_back(0);
    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }

    // Tensor blocks
    SharedMatrix Aia(new Matrix("Aia", naux, max_nia));
    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
    double** Aiap = Aia->pointer();
    double** Qiap = Qia->pointer();
    double** Jp   = Jm12->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {

        // Sizing
        ULI ia_start = ia_starts[block];
        ULI ia_stop  = ia_starts[block+1];
        ULI ncols = ia_stop - ia_start;

        // Read Aia
        timer_on("DFMP2 Aia Read");
        for (ULI Q = 0; Q < naux; Q++) {
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(Q*nia+ia_start));
            psio_->read(file,"(A|ia)",(char*)Aiap[Q],sizeof(double)*ncols,next_AIA,&next_AIA);
        }
        timer_off("DFMP2 Aia Read");

        // Apply Fitting
        timer_on("DFMP2 (Q|A)(A|ia)");
        C_DGEMM('T','N',ncols,naux,naux,1.0,Aiap[0],max_nia,Jp[0],naux,0.0,Qiap[0],naux);
        timer_off("DFMP2 (Q|A)(A|ia)");

        // Write Qia
        timer_on("DFMP2 Qia Write");
        psio_->write(file,"(Q|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
        timer_off("DFMP2 Qia Write");

    }
    psio_->close(file, 1);
}
void DFMP2::apply_fitting_grad(SharedMatrix Jm12, unsigned int file, ULI naux, ULI nia)
{
    // Memory constraints
    ULI Jmem = naux * naux;
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
    }
    ULI rem = (doubles - Jmem) / 2L;
    ULI max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<ULI> ia_starts;
    ia_starts.push_back(0);
    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }

    // Tensor blocks
    SharedMatrix Aia(new Matrix("Aia", max_nia, naux));
    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
    double** Aiap = Aia->pointer();
    double** Qiap = Qia->pointer();
    double** Jp   = Jm12->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {

        // Sizing
        ULI ia_start = ia_starts[block];
        ULI ia_stop  = ia_starts[block+1];
        ULI ncols = ia_stop - ia_start;

        // Read Qia
        timer_on("DFMP2 Qia Read");
        psio_->read(file,"(Q|ia)",(char*)Aiap[0],sizeof(double)*ncols*naux,next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        // Apply Fitting
        timer_on("DFMP2 (Q|A)(A|ia)");
        C_DGEMM('N','N',ncols,naux,naux,1.0,Aiap[0],naux,Jp[0],naux,0.0,Qiap[0],naux);
        timer_off("DFMP2 (Q|A)(A|ia)");

        // Write Bia
        timer_on("DFMP2 Bia Write");
        psio_->write(file,"(B|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
        timer_off("DFMP2 Bia Write");

    }
    psio_->close(file, 1);

}
void DFMP2::apply_gamma(unsigned int file, ULI naux, ULI nia)
{
    ULI Jmem = naux * naux;
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 1L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for gamma");
    }
    ULI rem = (doubles - Jmem) / 2L;
    ULI max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<ULI> ia_starts;
    ia_starts.push_back(0);
    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }

    // Tensor blocks
    SharedMatrix Aia(new Matrix("Aia", max_nia, naux));
    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
    SharedMatrix G(new Matrix("g", naux, naux));
    double** Aiap = Aia->pointer();
    double** Qiap = Qia->pointer();
    double** Gp   = G->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {

        // Sizing
        ULI ia_start = ia_starts[block];
        ULI ia_stop  = ia_starts[block+1];
        ULI ncols = ia_stop - ia_start;

        // Read Gia
        timer_on("DFMP2 Gia Read");
        psio_->read(file,"(G|ia)",(char*)Aiap[0],sizeof(double)*ncols*naux,next_AIA,&next_AIA);
        timer_off("DFMP2 Gia Read");

        // Read Cia
        timer_on("DFMP2 Cia Read");
        psio_->read(file,"(B|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
        timer_off("DFMP2 Cia Read");

        // g_PQ = G_ia^P C_ia^Q
        timer_on("DFMP2 g");
        C_DGEMM('T','N',naux,naux,ncols,1.0,Aiap[0],naux,Qiap[0],naux,1.0,Gp[0],naux); 
        timer_off("DFMP2 g");
    }

    psio_->write_entry(file, "G_PQ", (char*) Gp[0], sizeof(double) * naux * naux);

    psio_->close(file,1);
}
void DFMP2::apply_G_transpose(unsigned int file, ULI naux, ULI nia)
{
    // Memory constraints
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    ULI max_nia = (doubles / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<ULI> ia_starts;
    ia_starts.push_back(0);
    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }

    // Tensor blocks
    SharedMatrix Aia(new Matrix("Aia", naux, max_nia));
    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
    double** Aiap = Aia->pointer();
    double** Qiap = Qia->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {

        // Sizing
        ULI ia_start = ia_starts[block];
        ULI ia_stop  = ia_starts[block+1];
        ULI ncols = ia_stop - ia_start;

        // Read Gia
        timer_on("DFMP2 Gia Write");
        psio_->read(file,"(G|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
        timer_off("DFMP2 Gia Write");

        // Transpose
        for (int Q = 0; Q < naux; Q++) {
            C_DCOPY(ncols, &Qiap[0][Q], naux, Aiap[Q], 1);
        } 

        // Write Gia^\dagger
        timer_on("DFMP2 Aia Read");
        for (ULI Q = 0; Q < naux; Q++) {
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(Q*nia+ia_start));
            psio_->write(file,"(G|ia) T",(char*)Aiap[Q],sizeof(double)*ncols,next_AIA,&next_AIA);
        }
        timer_off("DFMP2 Aia Read");
    }
    psio_->close(file, 1);
}
void DFMP2::apply_B_transpose(unsigned int file, ULI naux, ULI naocc, ULI navir)
{
    // Memory constraints
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    ULI rows = doubles / (1L * naocc * naux);
    int max_A = (rows <= 0L ? 1L : rows);
    max_A = (rows > navir ? navir : rows);

    // Block sizing
    std::vector<int> a_starts;
    a_starts.push_back(0);
    for (int a = 0; a < navir; a+= max_A) { 
        if (a + max_A >= navir) {
            a_starts.push_back(navir);
        } else {
            a_starts.push_back(a + max_A);
        }
    }

    // Buffers
    SharedMatrix iaQ(new Matrix("iaQ", max_A * naocc, naux)); 
    double** iaQp = iaQ->pointer();
    
    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < a_starts.size() - 1; block++) {
    
        // Sizing
        int a_start = a_starts[block];
        int a_stop  = a_starts[block+1];
        int na = a_stop - a_start;
            
        // Read iaQ
        for (int a = 0; a < na; a++) {
            for (int i = 0; i < naocc; i++) {
                next_AIA = psio_get_address(PSIO_ZERO, sizeof(double) * (i * navir * naux + (a + a_start) * naux));
                psio_->read(file, "(Q|ia)", (char*) iaQp[a * naocc + i], sizeof(double) * naux, next_AIA, &next_AIA);
            }
        }

        // Write aiQ
        psio_->write(file, "(Q|ai)", (char*) iaQp[0], sizeof(double) * na * naocc * naux, next_QIA, &next_QIA);
    }
    psio_->close(file, 1);
}
void DFMP2::print_energies()
{
    energies_["Correlation Energy"] = energies_["Opposite-Spin Energy"] + energies_["Same-Spin Energy"] + energies_["Singles Energy"];
    energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

    energies_["SCS Opposite-Spin Energy"] = oss_*energies_["Opposite-Spin Energy"];
    energies_["SCS Same-Spin Energy"] = sss_*energies_["Same-Spin Energy"];
    energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"] + energies_["Singles Energy"];
    energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Singles Energy",           energies_["Singles Energy"]);
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

    // LAB TODO: drop DF- in labels to match DF-SCF behavior
    Process::environment.globals["CURRENT ENERGY"] = energies_["Total Energy"];
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = energies_["Correlation Energy"];
    Process::environment.globals["DF-MP2 TOTAL ENERGY"] = energies_["Total Energy"];
    Process::environment.globals["DF-MP2 SINGLES ENERGY"] = energies_["Singles Energy"];
    Process::environment.globals["DF-MP2 SAME-SPIN ENERGY"] = energies_["Same-Spin Energy"];
    Process::environment.globals["DF-MP2 OPPOSITE-SPIN ENERGY"] = energies_["Opposite-Spin Energy"];
    Process::environment.globals["DF-MP2 CORRELATION ENERGY"] = energies_["Correlation Energy"];
    Process::environment.globals["SCS-DF-MP2 TOTAL ENERGY"] = energies_["SCS Total Energy"];
    Process::environment.globals["SCS-DF-MP2 CORRELATION ENERGY"] = energies_["SCS Correlation Energy"];

}
void DFMP2::print_gradients()
{
    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Correlation");
    gradient_terms.push_back("Total");

    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients_.count(gradient_terms[i])) {
                gradients_[gradient_terms[i]]->print_atom_vector(); 
            }
        }
    } else {
        gradients_["Total"]->print();
    }

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
    Cfocc_ = Ca_subset("AO","FROZEN_OCC");
    Caocc_ = Ca_subset("AO","ACTIVE_OCC");
    Cavir_ = Ca_subset("AO","ACTIVE_VIR");
    Cfvir_ = Ca_subset("AO","FROZEN_VIR");

    eps_focc_ = epsilon_a_subset("AO","FROZEN_OCC");
    eps_aocc_ = epsilon_a_subset("AO","ACTIVE_OCC");
    eps_avir_ = epsilon_a_subset("AO","ACTIVE_VIR");
    eps_fvir_ = epsilon_a_subset("AO","FROZEN_VIR");
}
void RDFMP2::print_header()
{
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                          DF-MP2                         \n");
    fprintf(outfile, "\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    fprintf(outfile, "\t              RMP2 Wavefunction, %3d Threads             \n", nthread);
    fprintf(outfile, "\t                                                         \n");
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    fprintf(outfile, "\t           Ed Hohenstein, and C. David Sheriill          \n");
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
    // Schwarz Sieve
    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_,options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // ERI objects
    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("DF_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("DF_INTS_NUM_THREADS");
        }
    #endif

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(),
        basisset_,basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    ULI Amn_cost_per_row = nso * (ULI) nso;
    ULI Ami_cost_per_row = nso * (ULI) naocc;
    ULI Aia_cost_per_row = naocc * (ULI) navir;
    ULI total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    ULI max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(ribasis_->nshell());

    // Tensor blocks
    SharedMatrix Amn(new Matrix("(A|mn) Block", max_naux, nso * (ULI) nso));
    SharedMatrix Ami(new Matrix("(A|mi) Block", max_naux, nso * (ULI) naocc));
    SharedMatrix Aia(new Matrix("(A|ia) Block", max_naux, naocc * (ULI) navir));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    // C Matrices
    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();

    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = ribasis_->shell(Qstart).function_index();
        int nrows  = (Qstop == ribasis_->nshell() ?
                     ribasis_->nbf() -
                     ribasis_->shell(Qstart).function_index() :
                     ribasis_->shell(Qstop).function_index() -
                     ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*) Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        timer_on("DFMP2 (A|mn)");
        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = ribasis_->shell(Q).nfunction();
            int nm = basisset_->shell(M).nfunction();
            int nn = basisset_->shell(N).nfunction();

            int sq =  ribasis_->shell(Q).function_index();
            int sm =  basisset_->shell(M).function_index();
            int sn =  basisset_->shell(N).function_index();

            eri[thread]->compute_shell(Q,0,M,N);

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                        Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                        buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        timer_off("DFMP2 (A|mn)");

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc,nso,1.0,Amnp[0],nso,Caoccp[0],naocc,0.0,Amip[0],naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
        #pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T','N',naocc,navir,nso,1.0,Amip[row],naocc,Cavirp[0],navir,0.0,Aiap[row],navir);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_AIA,"(A|ia)",(char*)Aiap[0],sizeof(double)*nrows*naocc*navir,next_AIA,&next_AIA);
        timer_off("DFMP2 Aia Write");
    }

    psio_->close(PSIF_DFMP2_AIA,1);
}
void RDFMP2::form_Qia()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (ULI) Cavir_->colspi()[0]);
}
void RDFMP2::form_Qia_gradient()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (ULI) Cavir_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (ULI) Cavir_->colspi()[0]);
}
void RDFMP2::form_Qia_transpose()
{
    apply_B_transpose(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0], (ULI) Cavir_->colspi()[0]);
}
void RDFMP2::form_energy()
{
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI Iab_memory = navir * (ULI) navir;
    ULI Qa_memory  = naux  * (ULI) navir;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir,navir)));
    }

    double* eps_aoccp = eps_aocc_->pointer();
    double* eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Qia Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qia Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            timer_off("DFMP2 Qia Read");

            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_ss, e_os)
            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif
                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir,navir,naux,1.0,Qiap[(i-istart)*navir],naux,Qjbp[(j-jstart)*navir],naux,0.0,Iabp[0],navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                        e_ss += (iajb*iajb - iajb*ibja) * denom;
                        e_os += (iajb*iajb) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_AIA,0);

    energies_["Same-Spin Energy"] = e_ss;
    energies_["Opposite-Spin Energy"] = e_os;
}
void RDFMP2::form_Pab()
{
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    doubles -= navir * navir;
    double c = - (double) doubles;
    double b = 4.0 * navir * naux;
    double a = 2.0 * navir * (double) navir;

    int max_i = (int)((-b + sqrt(b*b - 4.0 * a * c)) / (2.0 * a));
    if (max_i <= 0) {
        throw PSIEXCEPTION("Not enough memory in DFMP2");
    }
    max_i = (max_i <= 0 ? 1 : max_i);
    max_i = (max_i > naocc ? naocc : max_i);    

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }

    // 2-Index Tensor blocks
    SharedMatrix Pab(new Matrix("Pab", navir, navir));
    double** Pabp = Pab->pointer();

    // 3-Index Tensor blocks
    SharedMatrix Qia(new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb(new Matrix("Qjb", max_i * (ULI) navir, naux));
    SharedMatrix Gia(new Matrix("Gia", max_i * (ULI) navir, naux));
    SharedMatrix Cjb(new Matrix("Cjb", max_i * (ULI) navir, naux));

    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();
    double** Giap = Gia->pointer();
    double** Cjbp = Cjb->pointer();

    // 4-index Tensor blocks
    SharedMatrix I(new Matrix("I", max_i * (ULI) navir, max_i * (ULI) navir));
    SharedMatrix T(new Matrix("T", max_i * (ULI) navir, max_i * (ULI) navir));
    double** Ip = I->pointer();
    double** Tp = T->pointer();

    ULI nIv = max_i * (ULI) navir; 

    double* eps_aoccp = eps_aocc_->pointer();
    double* eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_address next_AIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Qia Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        // Zero Gamma for current ia
        Gia->zero();

        for (int block_j = 0; block_j < i_starts.size() - 1; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qia Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            timer_off("DFMP2 Qia Read");

            // Read iaC chunk
            timer_on("DFMP2 Cia Read");
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
            psio_->read(PSIF_DFMP2_AIA,"(B|ia)",(char*)Cjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            timer_off("DFMP2 Cia Read");

            // Form the integrals (ia|jb) = B_ia^Q B_jb^Q    
            timer_on("DFMP2 I");
            C_DGEMM('N','T',ni * (ULI) navir, nj * (ULI) navir, naux, 1.0, Qiap[0], naux, Qjbp[0], naux, 0.0, Ip[0], nIv);
            timer_off("DFMP2 I");

            timer_on("DFMP2 T2");
            // Form the T amplitudes t_ia^jb = [2(ia|jb) - (ib|ja)] / (e_a + e_b - e_i - e_j)
            // Form the I amplitudes I_ia^jb = (ia|jb) / (e_a + e_b - e_i - e_j);
            // Form the energy contributions
            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_ss, e_os)
            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i_local = ij / nj;
                ULI j_local = ij % nj;
                ULI i = i_local + istart;
                ULI j = j_local + jstart;

                // Add the MP2 energy contributions and form the T amplitudes in place
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b <= a; b++) {
                        double iajb = Ip[i_local * navir + a][j_local * navir + b];
                        double ibja = Ip[i_local * navir + b][j_local * navir + a];
                        double denom = - 1.0 / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);
                        Tp[i_local * navir + a][j_local * navir + b] = denom * (2.0 * iajb - ibja);
                        Tp[i_local * navir + b][j_local * navir + a] = denom * (2.0 * ibja - iajb);
                        Ip[i_local * navir + a][j_local * navir + b] = denom * (iajb);
                        Ip[i_local * navir + b][j_local * navir + a] = denom * (ibja);

                        e_ss += (iajb*iajb - iajb*ibja) * denom;
                        e_os += (iajb*iajb) * denom;

                        if (a != b) {
                            e_ss += (ibja*ibja - ibja*iajb) * denom;
                            e_os += (ibja*ibja) * denom;
                        }
                    }
                }
            }
            timer_off("DFMP2 T2");

            // Form the Gamma tensor G_ia^P = t_ia^jb C_jb^P
            timer_on("DFMP2 G");
            C_DGEMM('N','N',ni * (ULI) navir, naux, nj * (ULI) navir, 2.0, Tp[0], nIv, Cjbp[0], naux, 1.0, Giap[0], naux); 
            timer_off("DFMP2 G");            

            // Sort the gimp column blocks, if gimp occurred. The idea is to get a contiguous iajb tensor
            if (nj != max_i) {
                unsigned long int counter = 0L;
                for (long int ind = 0L; ind < ni * (ULI) navir; ind++) {
                    ::memmove((void*)&Tp[0][counter],(void*)Tp[ind],sizeof(double) * nj * (ULI) navir);
                    ::memmove((void*)&Ip[0][counter],(void*)Ip[ind],sizeof(double) * nj * (ULI) navir);
                    counter += nj * (ULI) navir;
                }
            }

            // Form the virtual-virtual block of the density matrix
            timer_on("DFMP2 Pab");
            C_DGEMM('T','N',navir,navir,ni*(ULI)nj*navir,2.0,Tp[0],navir,Ip[0],navir,1.0,Pabp[0],navir);
            timer_off("DFMP2 Pab");
        }

        // Write iaG chunk
        timer_on("DFMP2 Gia Write");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->write(PSIF_DFMP2_AIA,"(G|ia)",(char*)Giap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Gia Write");
    }

    // Save the ab block of P
    psio_->write_entry(PSIF_DFMP2_AIA, "Pab", (char*) Pabp[0], sizeof(double) * navir * navir);
    //Pab->print();

    psio_->close(PSIF_DFMP2_AIA,1);

    energies_["Same-Spin Energy"] = e_ss;
    energies_["Opposite-Spin Energy"] = e_os;
}
void RDFMP2::form_Pij()
{
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    doubles -= naocc * naocc;
    double c = - (double) doubles;
    double b = 2.0 * naocc * naux;
    double a = 2.0 * naocc * (double) naocc;

    int max_a = (int)((-b + sqrt(b*b - 4.0 * a * c)) / (2.0 * a));
    if (max_a <= 0) {
        throw PSIEXCEPTION("Not enough memory in DFMP2");
    }
    max_a = (max_a <= 0 ? 1 : max_a);
    max_a = (max_a > navir ? navir : max_a);    

    // Blocks
    std::vector<ULI> a_starts;
    a_starts.push_back(0L);
    for (ULI a = 0; a < navir; a += max_a) {
        if (a + max_a >= navir) {
            a_starts.push_back(navir);
        } else {
            a_starts.push_back(a + max_a);
        }
    }

    // 2-Index Tensor blocks
    SharedMatrix Pij(new Matrix("Pij", naocc, naocc));
    double** Pijp = Pij->pointer();

    // 3-Index Tensor blocks
    SharedMatrix Qia(new Matrix("Qia", max_a * (ULI) naocc, naux));
    SharedMatrix Qjb(new Matrix("Qjb", max_a * (ULI) naocc, naux));

    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    // 4-index Tensor blocks
    SharedMatrix I(new Matrix("I", max_a * (ULI) naocc, max_a * (ULI) naocc));
    SharedMatrix T(new Matrix("T", max_a * (ULI) naocc, max_a * (ULI) naocc));
    double** Ip = I->pointer();
    double** Tp = T->pointer();

    ULI nVi = max_a * (ULI) naocc; 

    double* eps_aoccp = eps_aocc_->pointer();
    double* eps_avirp = eps_avir_->pointer();

    // Loop through pairs of blocks
    psio_address next_AIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    for (int block_a = 0; block_a < a_starts.size() - 1; block_a++) {

        // Sizing
        ULI astart = a_starts[block_a];
        ULI astop  = a_starts[block_a+1];
        ULI na     = astop - astart;

        // Read iaQ chunk
        timer_on("DFMP2 Qai Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(astart * naocc * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ai)",(char*)Qiap[0],sizeof(double)*(na * naocc * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qai Read");

        for (int block_b = 0; block_b < a_starts.size() - 1; block_b++) {

            // Sizing
            ULI bstart = a_starts[block_b];
            ULI bstop  = a_starts[block_b+1];
            ULI nb     = bstop - bstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qai Read");
            if (block_a == block_b) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(na * naocc * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(bstart * naocc * naux));
                psio_->read(PSIF_DFMP2_AIA,"(Q|ai)",(char*)Qjbp[0],sizeof(double)*(nb * naocc * naux),next_AIA,&next_AIA);
            }
            timer_off("DFMP2 Qai Read");

            // Form the integrals (ia|jb) = B_ia^Q B_jb^Q    
            timer_on("DFMP2 I");
            C_DGEMM('N','T',na * (ULI) naocc, nb * (ULI) naocc, naux, 1.0, Qiap[0], naux, Qjbp[0], naux, 0.0, Ip[0], nVi);
            timer_off("DFMP2 I");

            timer_on("DFMP2 T2");
            // Form the T amplitudes t_ia^jb = [2(ia|jb) - (ib|ja)] / (e_a + e_b - e_i - e_j)
            // Form the I amplitudes I_ia^jb = (ia|jb) / (e_a + e_b - e_i - e_j);
            // Form the energy contributions
            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_ss, e_os)
            for (long int ab = 0L; ab < na * nb; ab++) {

                // Sizing
                ULI a_local = ab / na;
                ULI b_local = ab % na;
                ULI a = a_local + astart;
                ULI b = b_local + bstart;

                // Add the MP2 energy contributions and form the T amplitudes in place
                for (int i = 0; i < naocc; i++) {
                    for (int j = 0; j <= i; j++) {
                        double iajb = Ip[a_local * naocc + i][b_local * naocc + j];
                        double ibja = Ip[a_local * naocc + j][b_local * naocc + i];
                        double denom = - 1.0 / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);
                        Tp[a_local * naocc + i][b_local * naocc + j] = denom * (2.0 * iajb - ibja);
                        Tp[a_local * naocc + j][b_local * naocc + i] = denom * (2.0 * ibja - iajb);
                        Ip[a_local * naocc + i][b_local * naocc + j] = denom * (iajb);
                        Ip[a_local * naocc + j][b_local * naocc + i] = denom * (ibja);

                        e_ss += (iajb*iajb - iajb*ibja) * denom;
                        e_os += (iajb*iajb) * denom;

                        if (i != j) {
                            e_ss += (ibja*ibja - ibja*iajb) * denom;
                            e_os += (ibja*ibja) * denom;
                        }
                    }
                }
            }
            timer_off("DFMP2 T2");

            // Sort the gimp column blocks, if gimp occurred. The idea is to get a contiguous iajb tensor
            if (nb != max_a) {
                unsigned long int counter = 0L;
                for (long int ind = 0L; ind < na * (ULI) naocc; ind++) {
                    ::memmove((void*)&Tp[0][counter],(void*)Tp[ind],sizeof(double) * nb * (ULI) naocc);
                    ::memmove((void*)&Ip[0][counter],(void*)Ip[ind],sizeof(double) * nb * (ULI) naocc);
                    counter += nb * (ULI) naocc;
                }
            }

            // Form the virtual-virtual block of the density matrix
            timer_on("DFMP2 Pij");
            C_DGEMM('T','N',naocc,naocc,na*(ULI)nb*naocc,-2.0,Tp[0],naocc,Ip[0],naocc,1.0,Pijp[0],naocc);
            timer_off("DFMP2 Pij");
        }
    }

    // Save the ab block of P
    psio_->write_entry(PSIF_DFMP2_AIA, "Pij", (char*) Pijp[0], sizeof(double) * naocc * naocc);
    //Pij->print();

    psio_->close(PSIF_DFMP2_AIA,1);
}
void RDFMP2::form_gamma()
{
    apply_gamma(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (ULI) Cavir_->colspi()[0]);
}
void RDFMP2::form_G_transpose()
{
    apply_G_transpose(PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_->colspi()[0] * (ULI) Cavir_->colspi()[0]);
}
void RDFMP2::form_AB_x_terms()
{
    // => Sizing <= //

    int natom = basisset_->molecule()->natom();
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();

    // => Gradient Contribution <= //

    gradients_["(A|B)^x"] = SharedMatrix(new Matrix("(A|B)^x Gradient", natom, 3));

    // => Forcing Terms/Gradients <= //

    SharedMatrix V(new Matrix("V", naux, naux));
    double** Vp = V->pointer();
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_->read_entry(PSIF_DFMP2_AIA, "G_PQ", (char*) Vp[0], sizeof(double) * naux * naux);
    psio_->close(PSIF_DFMP2_AIA, 1);

    // => Thread Count <= //    

    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif

    // => Integrals <= //

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(),ribasis_,BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < num_threads; t++) {
        Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));        
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Ktemps;
    for (int t = 0; t < num_threads; t++) {
        Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
    } 

    std::vector<std::pair<int,int> > PQ_pairs;
    for (int P = 0; P < ribasis_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int,int>(P,Q));    
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif 

        Jint[thread]->compute_shell_deriv1(P,0,Q,0);
        const double* buffer = Jint[thread]->buffer();
        
        int nP = ribasis_->shell(P).nfunction();
        int cP = ribasis_->shell(P).ncartesian();
        int aP = ribasis_->shell(P).ncenter();
        int oP = ribasis_->shell(P).function_index();

        int nQ = ribasis_->shell(Q).nfunction();
        int cQ = ribasis_->shell(Q).ncartesian();
        int aQ = ribasis_->shell(Q).ncenter();
        int oQ = ribasis_->shell(Q).function_index();

        int ncart = cP * cQ;
        const double *Px = buffer + 0*ncart;
        const double *Py = buffer + 1*ncart;
        const double *Pz = buffer + 2*ncart;
        const double *Qx = buffer + 3*ncart;
        const double *Qy = buffer + 4*ncart;
        const double *Qz = buffer + 5*ncart;

        double perm = (P == Q ? 1.0 : 2.0);

        double** grad_Kp = Ktemps[thread]->pointer();
    
        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) { 

                double Vval = perm * (0.5 * (Vp[p + oP][q + oQ] + Vp[q + oQ][p + oP]));;
                grad_Kp[aP][0] -= Vval * (*Px);
                grad_Kp[aP][1] -= Vval * (*Py);
                grad_Kp[aP][2] -= Vval * (*Pz);
                grad_Kp[aQ][0] -= Vval * (*Qx);
                grad_Kp[aQ][1] -= Vval * (*Qy);
                grad_Kp[aQ][2] -= Vval * (*Qz);

                Px++;
                Py++;
                Pz++;
                Qx++;
                Qy++;
                Qz++;
            }
        }
    }

    // => Temporary Gradient Reduction <= //

    for (int t = 0; t < num_threads; t++) {
        gradients_["(A|B)^x"]->add(Ktemps[t]);
    }
}
void RDFMP2::form_Amn_x_terms()
{
    // => Sizing <= //

    int natom = basisset_->molecule()->natom();
    int nso = basisset_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int nia = Caocc_->colspi()[0] * Cavir_->colspi()[0];
    int naux = ribasis_->nbf();

    // => ERI Sieve <= //

    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_,options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    int npairs = shell_pairs.size(); 
 
    // => Gradient Contribution <= //

    gradients_["(A|mn)^x"] = SharedMatrix(new Matrix("(A|mn)^x Gradient", natom, 3));

    // => Memory Constraints <= //

    ULI memory = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    int max_rows;
    int maxP = ribasis_->max_function_per_shell();
    ULI row_cost = 0L;
    row_cost += nso * (ULI) nso;
    row_cost += nso * (ULI) naocc;
    row_cost += naocc * (ULI) navir;
    ULI rows = memory / row_cost;
    rows = (rows > naux ? naux : rows);   
    rows = (rows < maxP ? maxP : rows);   
    max_rows = (int) rows;

    // => Block Sizing <= //

    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < ribasis_->nshell(); P++) {
        int nP = ribasis_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(ribasis_->nshell());

    // => Temporary Buffers <= //

    SharedMatrix Gia(new Matrix("Gia", max_rows, naocc * navir));
    SharedMatrix Gmi(new Matrix("Gmi", max_rows, nso * naocc));
    SharedMatrix Gmn(new Matrix("Gmn", max_rows, nso * (ULI) nso));

    double** Giap = Gia->pointer();
    double** Gmip = Gmi->pointer();
    double** Gmnp = Gmn->pointer();

    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();

    // => Thread Count <= //    

    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif

    // => Integrals <= //

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(),basisset_,basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < num_threads; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));        
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Ktemps;
    for (int t = 0; t < num_threads; t++) {
        Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
    } 

    // => PSIO <= //

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;

    // => Master Loop <= //
    
    for (int block = 0; block < Pstarts.size() - 1; block++) {
    
        // > Sizing < //
    
        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = ribasis_->shell(Pstart).function_index();
        int pstop  = (Pstop == ribasis_->nshell() ? naux : ribasis_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > G_ia^P -> G_mn^P < //
        
        psio_->read(PSIF_DFMP2_AIA, "(G|ia) T", (char*) Giap[0], sizeof(double) * np * nia, next_AIA, &next_AIA);
        
        #pragma omp parallel for num_threads(num_threads)
        for (int p = 0; p < np; p++) {
            C_DGEMM('N','T', nso, naocc, navir, 1.0, Cavirp[0], navir, Giap[p], navir, 0.0, Gmip[p], naocc);
        } 

        C_DGEMM('N','T', np * (ULI) nso, nso, naocc, 1.0, Gmip[0], naocc, Caoccp[0], naocc, 0.0, Gmnp[0], nso);

        // > Integrals < //
        #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;
         
            eri[thread]->compute_shell_deriv1(P,0,M,N); 

            const double* buffer = eri[thread]->buffer();

            int nP = ribasis_->shell(P).nfunction();
            int cP = ribasis_->shell(P).ncartesian();
            int aP = ribasis_->shell(P).ncenter();
            int oP = ribasis_->shell(P).function_index() - pstart;

            int nM = basisset_->shell(M).nfunction();
            int cM = basisset_->shell(M).ncartesian();
            int aM = basisset_->shell(M).ncenter();
            int oM = basisset_->shell(M).function_index();

            int nN = basisset_->shell(N).nfunction();
            int cN = basisset_->shell(N).ncartesian();
            int aN = basisset_->shell(N).ncenter();
            int oN = basisset_->shell(N).function_index();

            int ncart = cP * cM * cN;
            const double *Px = buffer + 0*ncart;
            const double *Py = buffer + 1*ncart;
            const double *Pz = buffer + 2*ncart;
            const double *Mx = buffer + 3*ncart;
            const double *My = buffer + 4*ncart;
            const double *Mz = buffer + 5*ncart;
            const double *Nx = buffer + 6*ncart;
            const double *Ny = buffer + 7*ncart;
            const double *Nz = buffer + 8*ncart;

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Kp = Ktemps[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        
                        double Jval = 2.0 * perm * (0.5 * (Gmnp[p + oP][(m + oM) * nso + (n + oN)] + Gmnp[p + oP][(n + oN) * nso + (m + oM)])); 
                        grad_Kp[aP][0] += Jval * (*Px);
                        grad_Kp[aP][1] += Jval * (*Py);
                        grad_Kp[aP][2] += Jval * (*Pz);
                        grad_Kp[aM][0] += Jval * (*Mx);
                        grad_Kp[aM][1] += Jval * (*My);
                        grad_Kp[aM][2] += Jval * (*Mz);
                        grad_Kp[aN][0] += Jval * (*Nx);
                        grad_Kp[aN][1] += Jval * (*Ny);
                        grad_Kp[aN][2] += Jval * (*Nz);

                        Px++; 
                        Py++; 
                        Pz++; 
                        Mx++; 
                        My++; 
                        Mz++; 
                        Nx++; 
                        Ny++; 
                        Nz++; 
                    }                   
                }
            }
        }
    }

    // => Temporary Gradient Reduction <= //

    for (int t = 0; t < num_threads; t++) {
        gradients_["(A|mn)^x"]->add(Ktemps[t]);
    }

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_L()
{
    // => Sizing <= //

    int nso = basisset_->nbf();
    int naocc = Caocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int nia = Caocc_->colspi()[0] * Cavir_->colspi()[0];
    int naux = ribasis_->nbf();

    // => ERI Sieve <= //

    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_,options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    int npairs = shell_pairs.size(); 
 
    // => Memory Constraints <= //

    ULI memory = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    memory -= naocc * nso;
    memory -= navir * nso;
    memory -= naocc * navir;
    int max_rows;
    int maxP = ribasis_->max_function_per_shell();
    ULI row_cost = 0L;
    row_cost += nso * (ULI) nso;
    row_cost += nso * (ULI) naocc;
    row_cost += nso * (ULI) navir;
    row_cost += naocc * (ULI) navir;
    ULI rows = memory / row_cost;
    rows = (rows > naux ? naux : rows);   
    rows = (rows < maxP ? maxP : rows);   
    max_rows = (int) rows;

    // => Block Sizing <= //

    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < ribasis_->nshell(); P++) {
        int nP = ribasis_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(ribasis_->nshell());

    // => Temporary Buffers <= //

    SharedMatrix Gia(new Matrix("Gia", max_rows, naocc * navir));
    SharedMatrix Gim(new Matrix("Pim", max_rows, nso * naocc));
    SharedMatrix Gam(new Matrix("Pam", max_rows, nso * navir));
    SharedMatrix Gmn(new Matrix("Pmn", max_rows, nso * (ULI) nso));

    double** Giap = Gia->pointer();
    double** Gimp = Gim->pointer();
    double** Gamp = Gam->pointer();
    double** Gmnp = Gmn->pointer();

    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();

    double* temp = new double[naocc * navir];

    // => Targets <= //

    SharedMatrix Lmi(new Matrix("Lma", nso, naocc));
    SharedMatrix Lma(new Matrix("Lmi", nso, navir));
    double** Lmip = Lmi->pointer();
    double** Lmap = Lma->pointer();

    // => Thread Count <= //    

    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif

    // => Integrals <= //

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(),basisset_,basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < num_threads; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri()));        
    }

    // => PSIO <= //

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;

    // => Master Loop <= //
    
    for (int block = 0; block < Pstarts.size() - 1; block++) {
    
        // > Sizing < //
    
        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = ribasis_->shell(Pstart).function_index();
        int pstop  = (Pstop == ribasis_->nshell() ? naux : ribasis_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > G_ia^P Read < //
        
        psio_->read(PSIF_DFMP2_AIA, "(G|ia) T", (char*) Giap[0], sizeof(double) * np * nia, next_AIA, &next_AIA);
        
        // > Integrals < //
        Gmn->zero();
        #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;
         
            eri[thread]->compute_shell(P,0,M,N); 

            const double* buffer = eri[thread]->buffer();

            int nP = ribasis_->shell(P).nfunction();
            int oP = ribasis_->shell(P).function_index() - pstart;

            int nM = basisset_->shell(M).nfunction();
            int oM = basisset_->shell(M).function_index();

            int nN = basisset_->shell(N).nfunction();
            int oN = basisset_->shell(N).function_index();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Gmnp[p + oP][(m + oM) * nso + (n + oN)] =     
                        Gmnp[p + oP][(n + oN) * nso + (m + oM)] =     
                        (*buffer++);
                    }                   
                }
            }
        }

        // => L_ma <= //
        
        #pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('T','N',naocc,nso,nso,1.0,Caoccp[0],naocc,Gmnp[p],nso,0.0,Gimp[p],nso);
        }

        C_DGEMM('T','N',nso,navir,naocc * (ULI) np, 1.0, Gimp[0], nso, Giap[0], navir, 1.0, Lmap[0], navir);

        // Sort G_P^ia to G_P^ai
        for (int p = 0; p < np; p++) {
            ::memcpy((void*) temp, (void*) Giap[p], sizeof(double) * naocc * navir);
            for (int i = 0; i < naocc; i++) {
                C_DCOPY(navir, &temp[i * navir], 1, &Giap[p][i], naocc);
            }
        }

        // => L_mi <= //
        
        #pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('T','N',navir,nso,nso,1.0,Cavirp[0],navir,Gmnp[p],nso,0.0,Gamp[p],nso);
        }

        C_DGEMM('T','N',nso,naocc,navir * (ULI) np, 1.0, Gamp[0], nso, Giap[0], naocc, 1.0, Lmip[0], naocc);
    }

    delete[] temp;

    //Lmi->print();
    //Lma->print();

    psio_->write_entry(PSIF_DFMP2_AIA, "Lmi", (char*) Lmip[0], sizeof(double) * nso * naocc);
    psio_->write_entry(PSIF_DFMP2_AIA, "Lma", (char*) Lmap[0], sizeof(double) * nso * navir);

    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_P()
{
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir; 
    
    // => Tensors <= //

    SharedMatrix Pij(new Matrix("Pij", naocc, naocc));
    SharedMatrix Pab(new Matrix("Pab", navir, navir));
    SharedMatrix PIj(new Matrix("PIj", nfocc, naocc));
    SharedMatrix PAb(new Matrix("PAb", nfvir, navir));
    SharedMatrix Ppq(new Matrix("Ppq", nmo, nmo)); 

    double** Pijp = Pij->pointer();
    double** Pabp = Pab->pointer();
    double** PIjp = PIj->pointer();
    double** PAbp = PAb->pointer();
    double** Ppqp = Ppq->pointer();
    
    SharedMatrix Lmi(new Matrix("Lmi", nso, naocc));
    SharedMatrix Lma(new Matrix("Lma", nso, navir));

    double** Lmip = Lmi->pointer();
    double** Lmap = Lma->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "Pij", (char*) Pijp[0], sizeof(double) * naocc * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "Pab", (char*) Pabp[0], sizeof(double) * navir * navir);
    psio_->read_entry(PSIF_DFMP2_AIA, "Lmi", (char*) Lmip[0], sizeof(double) * nso * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "Lma", (char*) Lmap[0], sizeof(double) * nso * navir);


    // => Occ-Occ/Virt-Virt <= //
    
    for (int i = 0; i < naocc; i++) {
        ::memcpy((void*) &Ppqp[nfocc + i][nfocc], (void*) Pijp[i], sizeof(double) * naocc);
    }
    
    for (int a = 0; a < navir; a++) {
        ::memcpy((void*) &Ppqp[nfocc + naocc + a][nfocc + naocc], (void*) Pabp[a], sizeof(double) * navir);
    }

    // => Frozen-Core/Virt <= //

    if (nfocc) {
        double** Cfoccp = Cfocc_->pointer();
        double*  eps_foccp = eps_focc_->pointer();
        double*  eps_aoccp = eps_aocc_->pointer();

        C_DGEMM('T','N',nfocc,naocc,nso,1.0,Cfoccp[0],nfocc,Lmip[0],naocc,0.0,PIjp[0],naocc);
        for (int i = 0; i < naocc; i++) {
            for (int J = 0; J < nfocc; J++) {
                PIjp[J][i] /= (eps_aoccp[i] - eps_foccp[J]); 
            }
        }

        for (int J = 0; J < nfocc; J++) {
            C_DCOPY(naocc, PIjp[J], 1, &Ppqp[J][nfocc], 1);
            C_DCOPY(naocc, PIjp[J], 1, &Ppqp[nfocc][J], nmo);
        }
    }

    if (nfvir) {
        double** Cfvirp = Cfvir_->pointer();
        double*  eps_fvirp = eps_fvir_->pointer();
        double*  eps_avirp = eps_avir_->pointer();

        C_DGEMM('T','N',nfvir,navir,nso,1.0,Cfvirp[0],nfvir,Lmap[0],navir,0.0,PAbp[0],navir);
        for (int b = 0; b < navir; b++) {
            for (int A = 0; A < nfvir; A++) {
                PAbp[A][b] /= -(eps_avirp[b] - eps_fvirp[A]); 
            }
        }

        for (int B = 0; B < nfocc; B++) {
            C_DCOPY(navir, PAbp[B], 1, &Ppqp[nfocc + naocc + navir + B][nfocc + naocc], 1);
            C_DCOPY(navir, PAbp[B], 1, &Ppqp[nfocc + naocc][nfocc + naocc + navir + B], nmo);
        }
    }

    //Pij->print();
    //Pab->print();
    //PIj->print();
    //PAb->print();
    //Ppq->print();

    // => Write out <= //

    psio_->write_entry(PSIF_DFMP2_AIA, "P", (char*) Ppqp[0], sizeof(double) * nmo * nmo);
    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_W()
{
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir; 
    
    // => Tensors <= //

    SharedMatrix Wpq1(new Matrix("Wpq1", nmo, nmo)); 
    double** Wpq1p = Wpq1->pointer();
    
    SharedMatrix Ppq(new Matrix("Ppq", nmo, nmo)); 
    double** Ppqp = Ppq->pointer();

    SharedMatrix Lmi(new Matrix("Lmi", nso, naocc));
    SharedMatrix Lma(new Matrix("Lma", nso, navir));
    SharedMatrix Lia(new Matrix("Lia", naocc + nfocc, navir + nfvir));

    double** Lmip = Lmi->pointer();
    double** Lmap = Lma->pointer();
    double** Liap = Lia->pointer();

    double** Cfoccp = Cfocc_->pointer();
    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();
    double** Cfvirp = Cfvir_->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P", (char*) Ppqp[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "Lmi", (char*) Lmip[0], sizeof(double) * nso * naocc);
    psio_->read_entry(PSIF_DFMP2_AIA, "Lma", (char*) Lmap[0], sizeof(double) * nso * navir);

    // => Term 1 <= //

    // RMP scales?

    // > Occ-Occ < //
    C_DGEMM('T','N',naocc,naocc,nso,-0.5,Caoccp[0],naocc,Lmip[0],naocc,0.0,&Wpq1p[nfocc][nfocc],nmo);
    if (nfocc) {
        C_DGEMM('T','N',nfocc,naocc,nso,-0.5,Cfoccp[0],nfocc,Lmip[0],naocc,0.0,&Wpq1p[0][nfocc],nmo);
    }   
 
    // > Vir-Vir < //
    C_DGEMM('T','N',navir,navir,nso,-0.5,Cavirp[0],navir,Lmap[0],navir,0.0,&Wpq1p[nfocc + naocc][nfocc + naocc],nmo);
    if (nfvir) {
        C_DGEMM('T','N',nfvir,navir,nso,-0.5,Cfvirp[0],nfvir,Lmap[0],navir,0.0,&Wpq1p[nfocc + naocc + navir][nfocc + navir],nmo);
    }   

    // > Occ-Vir < //
    C_DGEMM('T','N',naocc,navir,nso,-0.5,Caoccp[0],naocc,Lmap[0],navir,0.0,&Wpq1p[nfocc][nfocc + naocc], nmo);
    if (nfocc) {
        C_DGEMM('T','N',nfocc,navir,nso,-0.5,Cfoccp[0],nfocc,Lmap[0],navir,0.0,&Wpq1p[0][nfocc + naocc], nmo);
    }   

    // > Vir-Occ < //
    C_DGEMM('T','N',navir,naocc,nso,-0.5,Cavirp[0],navir,Lmip[0],naocc,0.0,&Wpq1p[nfocc + naocc][nfocc], nmo);
    if (nfvir) {
        C_DGEMM('T','N',nfvir,naocc,nso,-0.5,Cfvirp[0],nfvir,Lmip[0],naocc,0.0,&Wpq1p[nfocc + naocc + navir][nfocc], nmo);
    }   

    // => Lia (L contributions) <= //
    
    // Multiply by 2 to remove factor of 0.5 applied above
    for (int i = 0; i < (nfocc + naocc); i++) {
        for (int a = 0; a < (nfvir + navir); a++) {
            Liap[i][a] = 2.0 * (Wpq1p[i][a + naocc + nfocc] - Wpq1p[a + naocc + nfocc][i]);
        }
    }

    // > Symmetrize the result < //
    Wpq1->hermitivitize();
    Wpq1->scale(2.0); // Spin Integration
     
    // => Write-out <= //

    //Lia->print();
    //Wpq1->print();    

    psio_->write_entry(PSIF_DFMP2_AIA, "Lia", (char*) Liap[0], sizeof(double) * (naocc + nfocc) * (navir + nfvir));
    psio_->write_entry(PSIF_DFMP2_AIA, "W", (char*) Wpq1p[0], sizeof(double) * nmo * nmo);
    psio_->close(PSIF_DFMP2_AIA, 1);
}
void RDFMP2::form_Z()
{
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir; 
    int nocc = nfocc + naocc;    
    int nvir = nfvir + navir;    

    // => Tensors <= //

    SharedMatrix Wpq1(new Matrix("Wpq1", nmo, nmo)); 
    double** Wpq1p = Wpq1->pointer();
    SharedMatrix Wpq2(new Matrix("Wpq2", nmo, nmo)); 
    double** Wpq2p = Wpq2->pointer();
    SharedMatrix Wpq3(new Matrix("Wpq3", nmo, nmo)); 
    double** Wpq3p = Wpq3->pointer();
    
    SharedMatrix Ppq(new Matrix("Ppq", nmo, nmo)); 
    double** Ppqp = Ppq->pointer();

    SharedMatrix Lia(new Matrix("Lia", naocc + nfocc, navir + nfvir));
    double** Liap = Lia->pointer();

    SharedMatrix Cocc = Ca_subset("AO", "OCC");
    SharedMatrix Cvir = Ca_subset("AO", "VIR");
    SharedMatrix C    = Ca_subset("AO", "ALL");

    double** Cfoccp = Cfocc_->pointer();
    double** Caoccp = Caocc_->pointer();
    double** Cavirp = Cavir_->pointer();
    double** Cfvirp = Cfvir_->pointer();
    
    double** Coccp = Cocc->pointer();
    double** Cvirp = Cvir->pointer();
    double** Cp    = C->pointer();

    SharedVector eps_occ = epsilon_a_subset("AO", "OCC");
    SharedVector eps_vir = epsilon_a_subset("AO", "VIR");
    SharedVector eps     = epsilon_a_subset("AO", "ALL");
    double* epsp = eps->pointer();

    // => CPHF/JK Object <= //

    boost::shared_ptr<RCPHF> cphf(new RCPHF());
    cphf->set_C(C);
    cphf->set_Caocc(Cocc);
    cphf->set_Cavir(Cvir);
    cphf->set_eps_aocc(eps_occ);
    cphf->set_eps_avir(eps_vir);
    cphf->preiterations();

    boost::shared_ptr<JK> jk = cphf->jk();
    std::vector<SharedMatrix> & Cl = jk->C_left();
    std::vector<SharedMatrix> & Cr = jk->C_right();
    const std::vector<SharedMatrix> &J = jk->J();
    const std::vector<SharedMatrix> &K = jk->K();

    // Make sure the JK object doesn't try to backtransform our matrices
    jk->set_allow_desymmetrization(false);

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P", (char*) Ppqp[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "W", (char*) Wpq1p[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "Lia", (char*) Liap[0], sizeof(double) * (naocc + nfocc) * (navir + nfvir));

    // => Lia += 1/2 A_pqia P_pq (unrelaxed) <= //

    // > Factor the unrelaxed P^(2) (hopefully low rank) < //
    std::pair<SharedMatrix, SharedMatrix> factor1 = Ppq->partial_square_root(options_.get_double("DFMP2_P2_TOLERANCE"));
    SharedMatrix P1 = factor1.first;
    SharedMatrix N1 = factor1.second;
    double** P1p = P1->pointer();
    double** N1p = N1->pointer();
    //Ppq->print();
    //P1->print();
    //N1->print();

    // > Back-transform the transition orbitals < //
    SharedMatrix P1AO(new Matrix("P AO", nso, P1->colspi()[0]));
    SharedMatrix N1AO(new Matrix("N AO", nso, N1->colspi()[0])); 
    double** P1AOp = P1AO->pointer();
    double** N1AOp = N1AO->pointer();

    if (P1->colspi()[0]) {
        C_DGEMM('N','N',nso,P1->colspi()[0],nmo,1.0,Cp[0],nmo,P1p[0],P1->colspi()[0],0.0,P1AOp[0],P1->colspi()[0]);
    }

    if (N1->colspi()[0]) {
        C_DGEMM('N','N',nso,N1->colspi()[0],nmo,1.0,Cp[0],nmo,N1p[0],N1->colspi()[0],0.0,N1AOp[0],N1->colspi()[0]);
    }
    
    // > Form the J/K-like matrices (P,N contributions are separable) < //
    Cl.clear();
    Cr.clear();
    Cl.push_back(P1AO);
    Cl.push_back(N1AO);
    
    jk->compute();

    SharedMatrix J1 = J[0];
    SharedMatrix K1 = K[0]; 
    J1->subtract(J[1]);
    K1->subtract(K[1]);
    double** J1p = J1->pointer();
    double** K1p = K1->pointer();

    //J1->print();
    //K1->print();

    SharedMatrix AP(new Matrix("A_mn^ls P_ls^(2)", nso, nso));
    double** APp = AP->pointer();

    J1->scale(2.0);
    K1->scale(1.0);
    AP->add(J1);
    AP->subtract(K1);

    //AP->print();

    // > Form the contribution to Lia from the J/K-like matrices < //

    SharedMatrix T(new Matrix("T", nocc, nso));
    double** Tp = T->pointer();

    // L_ia += -1.0 (spin) C_mi { [ 4(mn|pq) - (mq|pn) - (mq|nq)] P_pq } C_na
    C_DGEMM('T','N',nocc,nso,nso,1.0,Coccp[0],nocc,APp[0],nso,0.0,Tp[0],nso);
    C_DGEMM('N','N',nocc,nvir,nso,1.0,Tp[0],nso,Cvirp[0],nvir,1.0,Liap[0],nvir);

    //Lia->print();

    // => (\delta_ij \delta_ab (\epsilon_a - \epsilon_i) + A_ia,jb) Z_jb = L_ia <= //

    std::map<std::string, SharedMatrix>&b = cphf->b();
    std::map<std::string, SharedMatrix>&x = cphf->x();

    b["Lia"] = Lia;
    cphf->compute_energy();
    SharedMatrix Zia = x["Lia"];
    Zia->scale(-1.0);

    // > Add Pia and Pai into the OPDM < // 
    SharedMatrix dPpq(new Matrix("dP", nmo, nmo));
    double** dPpqp = dPpq->pointer();
    double** Ziap = Zia->pointer();
    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nvir; a++) {
            dPpqp[i][a + nocc] = dPpqp[a + nocc][i] = Ziap[i][a];
        }
    }

    Ppq->add(dPpq);

    //Zia->print();
 
    // => Wik -= 1/2 A_pqik P_pq (relaxed) <= //

    // > Factor the relaxation contribution dP^(2) (hopefully low rank) < //
    std::pair<SharedMatrix, SharedMatrix> factor2 = dPpq->partial_square_root(options_.get_double("DFMP2_P2_TOLERANCE"));
    SharedMatrix P2 = factor2.first;
    SharedMatrix N2 = factor2.second;
    double** P2p = P2->pointer();
    double** N2p = N2->pointer();
    //dPpq->print();
    //P2->print();
    //N2->print();

    // > Back-transform the transition orbitals < //
    SharedMatrix P2AO(new Matrix("P AO", nso, P2->colspi()[0]));
    SharedMatrix N2AO(new Matrix("N AO", nso, N2->colspi()[0])); 
    double** P2AOp = P2AO->pointer();
    double** N2AOp = N2AO->pointer();

    if (P2->colspi()[0]) {
        C_DGEMM('N','N',nso,P2->colspi()[0],nmo,1.0,Cp[0],nmo,P2p[0],P2->colspi()[0],0.0,P2AOp[0],P2->colspi()[0]);
    }

    if (N2->colspi()[0]) {
        C_DGEMM('N','N',nso,N2->colspi()[0],nmo,1.0,Cp[0],nmo,N2p[0],N2->colspi()[0],0.0,N2AOp[0],N2->colspi()[0]);
    }
    
    // > Form the J/K-like matrices (P,N contributions are separable) < //
    Cl.clear();
    Cr.clear();
    Cl.push_back(P2AO);
    Cl.push_back(N2AO);
    
    jk->compute();
    
    SharedMatrix J2 = J[0]; 
    SharedMatrix K2 = K[0]; 
    J2->subtract(J[1]);
    K2->subtract(K[1]);
    double** J2p = J2->pointer();
    double** K2p = K2->pointer();

    J2->scale(2.0);
    K2->scale(1.0);
    AP->add(J2);
    AP->subtract(K2);

    // > Form the contribution to Lia from the J/K-like matrices < //

    // W_ik += +1.0 (spin) C_mi { [ 4(mn|pq) - (mq|pn) - (mq|nq)] P_pq } C_nk
    C_DGEMM('T','N',nocc,nso,nso,1.0,Coccp[0],nocc,APp[0],nso,0.0,Tp[0],nso);

    // occ-occ term
    C_DGEMM('N','N',nocc,nocc,nso,-1.0,Tp[0],nso,Coccp[0],nocc,0.0,&Wpq3p[0][0],nmo);

    // (RMP) what the frak is this interloper term doing in the plugin? Are the indices correct? why does this get 0.5, and the above get 1.0?
    C_DGEMM('N','N',nocc,nvir,nso,-0.5,Tp[0],nso,Cvirp[0],nvir,0.0,&Wpq3p[0][nocc],nmo);
    C_DGEMM('T','T',nvir,nocc,nso,-0.5,Cvirp[0],nvir,Tp[0],nso,0.0,&Wpq3p[nocc][0],nmo);

    // => W Term 2 <= //

    for (int p = 0; p < nmo; p++) {
        for (int q = 0; q < nmo; q++) {
            Wpq2p[p][q] = -0.5 * (epsp[p] + epsp[q]) * Ppqp[p][q];
        }
    }

    // => Final W <= //

    //Wpq1->print();
    //Wpq2->print();
    //Wpq3->print();

    Wpq1->add(Wpq2);
    Wpq1->add(Wpq3);
    Wpq1->set_name("Wpq");
        
    //Wpq1->print();
        
    psio_->write_entry(PSIF_DFMP2_AIA,"W",(char*) Wpq1p[0], sizeof(double) * nmo * nmo);

    // => Final P <= //
    
    //Ppq->print();

    psio_->write_entry(PSIF_DFMP2_AIA,"P",(char*) Ppqp[0], sizeof(double) * nmo * nmo);

    // => Finalize <= //

    psio_->close(PSIF_DFMP2_AIA,1);
    cphf->postiterations();
}
void RDFMP2::form_gradient()
{
    // => Sizing <= //

    int nso = basisset_->nbf();
    int nfocc = Cfocc_->colspi()[0];
    int navir = Cavir_->colspi()[0];
    int naocc = Caocc_->colspi()[0];
    int nfvir = Cfvir_->colspi()[0];
    int nmo = nfocc + naocc + navir + nfvir; 
    int nocc = nfocc + naocc;    
    int nvir = nfvir + navir;    

    // => Tensors <= //

    SharedMatrix W(new Matrix("W", nmo, nmo)); 
    double** Wp = W->pointer();
    
    SharedMatrix P2(new Matrix("P", nmo, nmo)); 
    double** P2p = P2->pointer();

    SharedMatrix Cocc = Ca_subset("AO", "OCC");
    SharedMatrix C    = Ca_subset("AO", "ALL");

    double** Coccp = Cocc->pointer();
    double** Cp    = C->pointer();

    SharedVector eps = epsilon_a_subset("AO", "ALL");
    double* epsp = eps->pointer();

    // => Read-in <= //

    psio_->open(PSIF_DFMP2_AIA, 1);
    psio_->read_entry(PSIF_DFMP2_AIA, "P", (char*) P2p[0], sizeof(double) * nmo * nmo);
    psio_->read_entry(PSIF_DFMP2_AIA, "W", (char*) Wp[0], sizeof(double) * nmo * nmo);

    // => Dress for SCF <= //

    SharedMatrix P2F(P2->clone());
    double** P2Fp = P2F->pointer();
    P2F->scale(2.0);

    //W->zero();
    //P2->zero();
    W->scale(-1.0);
    for (int i = 0; i < nocc; i++) {
        Wp[i][i] += 2.0 * epsp[i];
        P2p[i][i] += 2.0;
        P2Fp[i][i] += 2.0;
    }

    //P2->print();
    //W->print();

    psio_->write_entry(PSIF_DFMP2_AIA, "P", (char*) P2p[0], sizeof(double) * nmo * nmo);
    psio_->write_entry(PSIF_DFMP2_AIA, "W", (char*) Wp[0], sizeof(double) * nmo * nmo);

    // => Factorize the P matrix <= //

    P2F->scale(0.5);
    std::pair<SharedMatrix, SharedMatrix> factor = P2F->partial_square_root(options_.get_double("DFMP2_P_TOLERANCE"));
    P2F->scale(2.0);

    SharedMatrix P1 = factor.first;
    SharedMatrix N1 = factor.second;
    double** P1p = P1->pointer();
    double** N1p = N1->pointer();

    // => Back-transform <= //
    
    SharedMatrix T1(new Matrix("T", nmo, nso));
    SharedMatrix PAO(new Matrix("P AO", nso, nso));
    SharedMatrix PFAO(new Matrix("PF AO", nso, nso));
    SharedMatrix WAO(new Matrix("W AO", nso, nso));
    SharedMatrix P1AO(new Matrix("P1 AO", nso, P1->colspi()[0]));
    SharedMatrix N1AO(new Matrix("N1 AO", nso, N1->colspi()[0]));

    double** T1p = T1->pointer();
    double** PAOp = PAO->pointer();
    double** PFAOp = PFAO->pointer();
    double** WAOp = WAO->pointer();
    double** P1AOp = P1AO->pointer();
    double** N1AOp = N1AO->pointer();

    C_DGEMM('N','T',nmo,nso,nmo,1.0,P2p[0],nmo,Cp[0],nmo,0.0,T1p[0],nso);
    C_DGEMM('N','N',nso,nso,nmo,1.0,Cp[0],nmo,T1p[0],nmo,0.0,PAOp[0],nso);
    
    C_DGEMM('N','T',nmo,nso,nmo,1.0,P2Fp[0],nmo,Cp[0],nmo,0.0,T1p[0],nso);
    C_DGEMM('N','N',nso,nso,nmo,1.0,Cp[0],nmo,T1p[0],nmo,0.0,PFAOp[0],nso);
    
    C_DGEMM('N','T',nmo,nso,nmo,1.0,Wp[0],nmo,Cp[0],nmo,0.0,T1p[0],nso);
    C_DGEMM('N','N',nso,nso,nmo,1.0,Cp[0],nmo,T1p[0],nmo,0.0,WAOp[0],nso);
    
    //P2->print();
    //PAO->print();

    if (P1->colspi()[0]) {
        C_DGEMM('N','N',nso,P1->colspi()[0],nmo,1.0,Cp[0],nmo,P1p[0],P1->colspi()[0],0.0,P1AOp[0],P1->colspi()[0]);
    }
    
    if (N1->colspi()[0]) {
        C_DGEMM('N','N',nso,N1->colspi()[0],nmo,1.0,Cp[0],nmo,N1p[0],N1->colspi()[0],0.0,N1AOp[0],N1->colspi()[0]);
    }

    // => Random extra drivers for the JK gradients <= //

    SharedMatrix D(PAO->clone()); 
    double** Dp = D->pointer();
    C_DGEMM('N','T',nso,nso,nocc,2.0,Coccp[0],nocc,Coccp[0],nocc,0.0,Dp[0],nso);    

    //D->print();

    SharedMatrix Ds(D->clone());
    Ds->scale(0.5);

    SharedMatrix PFAOs(PFAO->clone());
    PFAOs->scale(0.5);

    // => Gogo Gradients <= //

    std::vector<std::string> gradient_terms;
    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("Coulomb");
    gradient_terms.push_back("Exchange");
    gradient_terms.push_back("Correlation");
    gradient_terms.push_back("Total");

    // => Sizings <= //
    int natom = molecule_->natom();
     
    // => Nuclear Gradient <= //
    gradients_["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
    gradients_["Nuclear"]->set_name("Nuclear Gradient");

    // => Kinetic Gradient <= //
    timer_on("Grad: T");
    {
        gradients_["Kinetic"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Kinetic"]->set_name("Kinetic Gradient");
        gradients_["Kinetic"]->zero();
        double** Tp = gradients_["Kinetic"]->pointer();

        // Kinetic derivatives
        boost::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
        const double* buffer = Tint->buffer();   

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Tint->compute_shell_deriv1(P,Q);
                                
                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();
 
                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);
               
                // Px 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][0] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Py 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][1] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Pz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][2] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][0] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][1] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][2] += perm * PAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        } 
    }
    timer_off("Grad: T");

    // => Potential Gradient <= //
    timer_on("Grad: V");
    {
        gradients_["Potential"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Potential"]->set_name("Potential Gradient");
        gradients_["Potential"]->zero();

        // Thread count
        int threads = 1;
        #ifdef _OPENMP
            threads = omp_get_max_threads();
        #endif

        // Potential derivatives
        std::vector<boost::shared_ptr<OneBodyAOInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        for (int t = 0; t < threads; t++) { 
            Vint.push_back(boost::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
            Vtemps.push_back(SharedMatrix(gradients_["Potential"]->clone()));
        }
       
        // Lower Triangle
        std::vector<std::pair<int,int> > PQ_pairs;
        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {
                PQ_pairs.push_back(std::pair<int,int>(P,Q));
            }
        }

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

            int P = PQ_pairs[PQ].first;
            int Q = PQ_pairs[PQ].second;

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            Vint[thread]->compute_shell_deriv1(P,Q);
            const double* buffer = Vint[thread]->buffer();
                            
            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();
 
            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            double perm = (P == Q ? 1.0 : 2.0);
                
            double** Vp = Vtemps[thread]->pointer();
            
            for (int A = 0; A < natom; A++) {
                const double* ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
                const double* ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
                const double* ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        double Vval = perm * PAOp[p + oP][q + oQ];
                        Vp[A][0] += Vval * (*ref0++);
                        Vp[A][1] += Vval * (*ref1++);
                        Vp[A][2] += Vval * (*ref2++);
                    }
                }
            }
        } 
    
        for (int t = 0; t < threads; t++) { 
            gradients_["Potential"]->add(Vtemps[t]);
        }
    }
    timer_off("Grad: V");

    //gradients_["One-Electron"] = SharedMatrix(gradients_["Nuclear"]->clone());
    //gradients_["One-Electron"]->set_name("One-Electron Gradient");
    //gradients_["One-Electron"]->zero();
    //gradients_["One-Electron"]->add(gradients_["Kinetic"]);
    //gradients_["One-Electron"]->add(gradients_["Potential"]);
    //gradients_["One-Electron"]->print();

    // => Overlap Gradient <= //
    timer_on("Grad: S");
    {
        gradients_["Overlap"] = SharedMatrix(gradients_["Nuclear"]->clone());
        gradients_["Overlap"]->set_name("Overlap Gradient");
        gradients_["Overlap"]->zero();
        double** Sp = gradients_["Overlap"]->pointer();

        // Overlap derivatives
        boost::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
        const double* buffer = Sint->buffer();   

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Sint->compute_shell_deriv1(P,Q);
                                
                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();
 
                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);
               
                // Px 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][0] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Py 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][1] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Pz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][2] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][0] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Qy 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][1] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
               
                // Qz 
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][2] -= perm * WAOp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        } 
    }
    timer_off("Grad: S");

    // => Two-Electron Gradient <= //
    
    timer_on("Grad: JK");
    
    boost::shared_ptr<CorrGrad> jk = CorrGrad::build_CorrGrad();
    jk->set_memory((ULI) (options_.get_double("SCF_MEM_SAFETY_FACTOR") * memory_ / 8L));

    jk->set_Ca(Cocc);
    jk->set_Cb(Cocc);
    jk->set_La(P1AO); 
    jk->set_Lb(P1AO); 
    jk->set_Ra(N1AO);
    jk->set_Rb(N1AO);
    jk->set_Da(Ds);
    jk->set_Db(Ds);
    jk->set_Dt(D);
    jk->set_Pa(PFAOs);
    jk->set_Pb(PFAOs);
    jk->set_Pt(PFAO); 

    jk->print_header();
    jk->compute_gradient();    

    std::map<std::string, SharedMatrix>& jk_gradients = jk->gradients();
    gradients_["Coulomb"] = jk_gradients["Coulomb"];
    gradients_["Exchange"] = jk_gradients["Exchange"];
    gradients_["Exchange"]->scale(-1.0);

    //gradients_["Separable TPDM"] = SharedMatrix(gradients_["Nuclear"]->clone());
    //gradients_["Separable TPDM"]->set_name("Separable TPDM Gradient");
    //gradients_["Separable TPDM"]->zero();
    //gradients_["Separable TPDM"]->add(gradients_["Coulomb"]);
    //gradients_["Separable TPDM"]->add(gradients_["Exchange"]);
    //gradients_["Separable TPDM"]->print();

    timer_off("Grad: JK");

    // => Correlation Gradient (Previously computed) <= //

    SharedMatrix correlation = SharedMatrix(gradients_["Nuclear"]->clone());
    correlation->zero();
    correlation->add(gradients_["(A|mn)^x"]); 
    correlation->add(gradients_["(A|B)^x"]); 
    gradients_["Correlation"] = correlation;
    gradients_["Correlation"]->set_name("Correlation Gradient");

    // => Total Gradient <= //
    SharedMatrix total = SharedMatrix(gradients_["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients_.count(gradient_terms[i])) {
            total->add(gradients_[gradient_terms[i]]); 
        }
    }

    gradients_["Total"] = total; 
    gradients_["Total"]->set_name("Total Gradient");

    // => Finalize <= //

    psio_->close(PSIF_DFMP2_AIA,1);
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
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                          DF-MP2                         \n");
    fprintf(outfile, "\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    fprintf(outfile, "\t              UMP2 Wavefunction, %3d Threads             \n", nthread);
    fprintf(outfile, "\t                                                         \n");
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    fprintf(outfile, "\t           Ed Hohenstein, and C. David Sheriill          \n");
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
    // Schwarz Sieve
    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_,options_.get_double("INTS_TOLERANCE")));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // ERI objects
    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("DF_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("DF_INTS_NUM_THREADS");
        }
    #endif

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(),
        basisset_,basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc_a = Caocc_a_->colspi()[0];
    int navir_a = Cavir_a_->colspi()[0];
    int naocc_b = Caocc_b_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
    int navir = (navir_a > navir_b ? navir_a : navir_b);
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    ULI Amn_cost_per_row = nso * (ULI) nso;
    ULI Ami_cost_per_row = nso * (ULI) naocc;
    ULI Aia_cost_per_row = naocc * (ULI) navir;
    ULI total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    ULI max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(ribasis_->nshell());

    // Tensor blocks
    SharedMatrix Amn(new Matrix("(A|mn) Block", max_naux, nso * (ULI) nso));
    SharedMatrix Ami(new Matrix("(A|mi) Block", max_naux, nso * (ULI) naocc));
    SharedMatrix Aia(new Matrix("(A|ia) Block", max_naux, naocc * (ULI) navir));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    // C Matrices
    double** Caoccap = Caocc_a_->pointer();
    double** Cavirap = Cavir_a_->pointer();
    double** Caoccbp = Caocc_b_->pointer();
    double** Cavirbp = Cavir_b_->pointer();

    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_NEW);
    psio_address next_QIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = ribasis_->shell(Qstart).function_index();
        int nrows  = (Qstop == ribasis_->nshell() ?
                     ribasis_->nbf() -
                     ribasis_->shell(Qstart).function_index() :
                     ribasis_->shell(Qstop).function_index() -
                     ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*) Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        timer_on("DFMP2 (A|mn)");
        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = ribasis_->shell(Q).nfunction();
            int nm = basisset_->shell(M).nfunction();
            int nn = basisset_->shell(N).nfunction();

            int sq =  ribasis_->shell(Q).function_index();
            int sm =  basisset_->shell(M).function_index();
            int sn =  basisset_->shell(N).function_index();

            eri[thread]->compute_shell(Q,0,M,N);

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                        Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                        buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        timer_off("DFMP2 (A|mn)");

        // => Alpha Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc_a,nso,1.0,Amnp[0],nso,Caoccap[0],naocc_a,0.0,Amip[0],naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
        #pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T','N',naocc_a,navir_a,nso,1.0,Amip[row],naocc,Cavirap[0],navir_a,0.0,&Aiap[0][row*(ULI)naocc_a*navir_a],navir_a);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_AIA,"(A|ia)",(char*)Aiap[0],sizeof(double)*nrows*naocc_a*navir_a,next_AIA,&next_AIA);
        timer_off("DFMP2 Aia Write");

        // => Beta Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("DFMP2 (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc_b,nso,1.0,Amnp[0],nso,Caoccbp[0],naocc_b,0.0,Amip[0],naocc);
        timer_off("DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("DFMP2 (A|mi)C_na");
        #pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T','N',naocc_b,navir_b,nso,1.0,Amip[row],naocc,Cavirbp[0],navir_b,0.0,&Aiap[0][row*(ULI)naocc_b*navir_b],navir_b);
        }
        timer_off("DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_QIA,"(A|ia)",(char*)Aiap[0],sizeof(double)*nrows*naocc_b*navir_b,next_QIA,&next_QIA);
        timer_off("DFMP2 Aia Write");
    }

    psio_->close(PSIF_DFMP2_AIA,1);
    psio_->close(PSIF_DFMP2_QIA,1);
}
void UDFMP2::form_Qia()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (ULI) Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (ULI) Cavir_b_->colspi()[0]);
}
void UDFMP2::form_Qia_gradient()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (ULI) Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (ULI) Cavir_b_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_AIA, ribasis_->nbf(), Caocc_a_->colspi()[0] * (ULI) Cavir_a_->colspi()[0]);
    apply_fitting_grad(Jm12, PSIF_DFMP2_QIA, ribasis_->nbf(), Caocc_b_->colspi()[0] * (ULI) Cavir_b_->colspi()[0]);
}
void UDFMP2::form_Qia_transpose()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_energy()
{
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    /* => AA Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_a_->colspi()[0];
    int navir = Cavir_a_->colspi()[0];

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI Iab_memory = navir * (ULI) navir;
    ULI Qa_memory  = naux  * (ULI) navir;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir,navir)));
    }

    double* eps_aoccp = eps_aocc_a_->pointer();
    double* eps_avirp = eps_avir_a_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Qia Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qia Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            timer_off("DFMP2 Qia Read");

            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_ss)
            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif
                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir,navir,naux,1.0,Qiap[(i-istart)*navir],naux,Qjbp[(j-jstart)*navir],naux,0.0,Iabp[0],navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                        e_ss += 0.5*(iajb*iajb - iajb*ibja) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_AIA,1);

    /* End AA Terms */ }

    /* => BB Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_b_->colspi()[0];
    int navir = Cavir_b_->colspi()[0];

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI Iab_memory = navir * (ULI) navir;
    ULI Qa_memory  = naux  * (ULI) navir;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir,navir)));
    }

    double* eps_aoccp = eps_aocc_b_->pointer();
    double* eps_avirp = eps_avir_b_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Qia Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            timer_on("DFMP2 Qia Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            timer_off("DFMP2 Qia Read");

            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_ss)
            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif
                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir,navir,naux,1.0,Qiap[(i-istart)*navir],naux,Qjbp[(j-jstart)*navir],naux,0.0,Iabp[0],navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                        e_ss += 0.5*(iajb*iajb - iajb*ibja) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_QIA,1);

    /* End BB Terms */ }

    /* => AB Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc_a = Caocc_a_->colspi()[0];
    int navir_a = Cavir_a_->colspi()[0];
    int naocc_b = Caocc_b_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
    int navir = (navir_a > navir_b ? navir_a : navir_b);

    // Thread considerations
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Memory
    ULI Iab_memory = navir_a * (ULI) navir_b;
    ULI Qa_memory  = naux  * (ULI) navir_a;
    ULI Qb_memory  = naux  * (ULI) navir_b;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (Qa_memory + Qb_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts_a;
    i_starts_a.push_back(0L);
    for (ULI i = 0; i < naocc_a; i += max_i) {
        if (i + max_i >= naocc_a) {
            i_starts_a.push_back(naocc_a);
        } else {
            i_starts_a.push_back(i + max_i);
        }
    }
    std::vector<ULI> i_starts_b;
    i_starts_b.push_back(0L);
    for (ULI i = 0; i < naocc_b; i += max_i) {
        if (i + max_i >= naocc_b) {
            i_starts_b.push_back(naocc_b);
        } else {
            i_starts_b.push_back(i + max_i);
        }
    }

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir_a, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir_b, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir_a,navir_b)));
    }

    double* eps_aoccap = eps_aocc_a_->pointer();
    double* eps_avirap = eps_avir_a_->pointer();
    double* eps_aoccbp = eps_aocc_b_->pointer();
    double* eps_avirbp = eps_avir_b_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts_a.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts_a[block_i];
        ULI istop  = i_starts_a[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        timer_on("DFMP2 Qia Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir_a * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir_a * naux),next_AIA,&next_AIA);
        timer_off("DFMP2 Qia Read");

        for (int block_j = 0; block_j < i_starts_b.size() - 1; block_j++) {

            // Sizing
            ULI jstart = i_starts_b[block_j];
            ULI jstop  = i_starts_b[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk
            timer_on("DFMP2 Qia Read");
            next_QIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir_b * naux));
            psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir_b * naux),next_QIA,&next_QIA);
            timer_off("DFMP2 Qia Read");

            #pragma omp parallel for schedule(dynamic) num_threads(nthread) reduction(+: e_os)
            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;

                // Which thread is this?
                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif
                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir_a,navir_b,naux,1.0,Qiap[(i-istart)*navir_a],naux,Qjbp[(j-jstart)*navir_b],naux,0.0,Iabp[0],navir_b);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir_a; a++) {
                    for (int b = 0; b < navir_b; b++) {
                        double iajb = Iabp[a][b];
                        double denom = - 1.0 / (eps_avirap[a] + eps_avirbp[b] - eps_aoccap[i] - eps_aoccbp[j]);
                        e_os += (iajb*iajb) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_AIA,0);
    psio_->close(PSIF_DFMP2_QIA,0);

    /* End BB Terms */ }

    energies_["Same-Spin Energy"] = e_ss;
    energies_["Opposite-Spin Energy"] = e_os;
}
void UDFMP2::form_Pab()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_Pij()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_gamma()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_G_transpose()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_AB_x_terms()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_Amn_x_terms()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_L()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_P()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_W()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_Z()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
}
void UDFMP2::form_gradient()
{
    throw PSIEXCEPTION("UDFMP2: Gradients not yet implemented");
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
    fprintf(outfile, "\t                          DF-MP2                         \n");
    fprintf(outfile, "\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    fprintf(outfile, "\t          ROHF-MBPT(2) Wavefunction, %3d Threads         \n", nthread);
    fprintf(outfile, "\t                                                         \n");
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    fprintf(outfile, "\t           Ed Hohenstein, and C. David Sheriill          \n");
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

