#include "mp2.h"
#include <lib3index/3index.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
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
    
    if (options_.get_str("REFERENCE") == "ROHF")
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
    fprintf(outfile, "\t                          DF-MP2                         \n");
    fprintf(outfile, "\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    fprintf(outfile, "\t              RMP2 Wavefunction, %3d Threads             \n", nthread);
    fprintf(outfile, "\t                                                         \n");
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonet,       \n");
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
        int nQ = ribasis_->shell(Q)->nfunction();
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
        int qoff   = ribasis_->shell(Qstart)->function_index();
        int nrows  = (Qstop == ribasis_->nshell() ? 
                     ribasis_->nbf() - 
                     ribasis_->shell(Qstart)->function_index() : 
                     ribasis_->shell(Qstop)->function_index() - 
                     ribasis_->shell(Qstart)->function_index()); 

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
        
            int nq = ribasis_->shell(Q)->nfunction();       
            int nm = basisset_->shell(M)->nfunction();       
            int nn = basisset_->shell(N)->nfunction();       
        
            int sq =  ribasis_->shell(Q)->function_index();       
            int sm =  basisset_->shell(M)->function_index();       
            int sn =  basisset_->shell(N)->function_index();       

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
    ULI Iab_memory = naocc * (ULI) navir;
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
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonet,       \n");
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
        int nQ = ribasis_->shell(Q)->nfunction();
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
        int qoff   = ribasis_->shell(Qstart)->function_index();
        int nrows  = (Qstop == ribasis_->nshell() ? 
                     ribasis_->nbf() - 
                     ribasis_->shell(Qstart)->function_index() : 
                     ribasis_->shell(Qstop)->function_index() - 
                     ribasis_->shell(Qstart)->function_index()); 

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
        
            int nq = ribasis_->shell(Q)->nfunction();       
            int nm = basisset_->shell(M)->nfunction();       
            int nn = basisset_->shell(N)->nfunction();       
        
            int sq =  ribasis_->shell(Q)->function_index();       
            int sm =  basisset_->shell(M)->function_index();       
            int sn =  basisset_->shell(N)->function_index();       

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
    fprintf(outfile, "\t        Rob Parrish, Justin Turney, Andy Simmonet,       \n");
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

