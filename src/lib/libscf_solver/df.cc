#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psiconfig.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"
#include "df.h"
#include <lib3index/3index.h>

//MKL Header
#include <psiconfig.h>
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

#include <libmints/mints.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

// ==> Constructors / Common Init <== //

DFHF::DFHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt) :
    primary_(basis), psio_(psio), options_(opt), unit_(PSIF_DFSCF_BJ), omega_(0.0)
{
    common_init();
}
DFHF::DFHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt, double omega) :
    primary_(basis), psio_(psio), options_(opt), unit_(PSIF_DFSCF_BJ), omega_(omega)
{
    common_init();
}
DFHF::~DFHF()
{
}
void DFHF::common_init()
{
    print_ = 1;
    if (options_["PRINT"].has_changed())
        print_ = options_.get_int("PRINT");

    if (print_) {
        fprintf(outfile, " DFHF: Density-Fitted SCF Algorithms.\n");
        if (omega_ > 0.0) {
            fprintf(outfile, "   Long-Range Omega = %11.8E\n",omega_);
        }
        fprintf(outfile, "   by Rob Parrish\n\n");
    }

    is_initialized_ = false;
    is_disk_initialized_ = false;
    is_jk_ = false;
    restricted_ = false;

    is_disk_ = false;
    memory_ = Process::environment.get_memory() / 8L;
    memory_ = (unsigned long int) (0.7 * memory_);

    // Build auxiliary basis from options
    zero_ = BasisSet::zero_ao_basis_set();

    // If the user doesn't spec a basis name, pick it yourself
    // TODO: Verify that the basis assign does not messs this up
    if (options_.get_str("RI_BASIS_SCF") == "") {
        primary_->molecule()->set_basis_all_atoms(options_.get_str("BASIS") + "-JKFIT", "RI_BASIS_SCF");
        fprintf(outfile, "  No auxiliary basis selected, defaulting to %s-JKFIT\n\n", options_.get_str("BASIS").c_str());
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, primary_->molecule(), "RI_BASIS_SCF");
    parser.reset();

    if (print_) {
        fprintf(outfile, "  ==> Auxiliary Basis: %s <==\n\n", options_.get_str("RI_BASIS_SCF").c_str());
        auxiliary_->print_by_level(outfile, print_);
    }

    timer_on("Schwarz");
    schwarz_ = boost::shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, options_.get_double("SCHWARZ_CUTOFF")));
    timer_off("Schwarz");

    if (print_ > 1) {
        int ntri_shell = schwarz_->get_nshell_pairs();
        int ntri_shell_naive = primary_->nshell() * (primary_->nshell() + 1) / 2;
        int ntri_fun = schwarz_->get_nfun_pairs();
        int ntri_fun_naive = primary_->nbf() * (primary_->nbf() + 1) / 2;

        double ntri_shell_savings = 1.0 - ntri_shell / (double) ntri_shell_naive;
        double ntri_fun_savings   = 1.0 - ntri_fun / (double) ntri_fun_naive;

        fprintf(outfile, "  Schwarz Cutoff is %8.3E\n", options_.get_double("SCHWARZ_CUTOFF"));
        fprintf(outfile, "   -Shell Pairs:    %12d of %12d remain, %3.0f%% savings\n",
            ntri_shell, ntri_shell_naive, 100.0*ntri_shell_savings);
        fprintf(outfile, "   -Function Pairs: %12d of %12d remain, %3.0f%% savings\n\n",
            ntri_fun, ntri_fun_naive, 100.0*ntri_fun_savings);
    }

    if (omega_ > 0.0) {
        Jinv_ = boost::shared_ptr<FittingMetric>(new FittingMetric(auxiliary_, true, omega_));
    } else {
        Jinv_ = boost::shared_ptr<FittingMetric>(new FittingMetric(auxiliary_, true));
    }

    // Make a memory decision here
    is_disk_ = false;
    memory_ = Process::environment.get_memory() / 8L;
    memory_ = (unsigned long int) (0.7 * memory_);

    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    // Two is for buffer space in fitting
    is_disk_ = (three_memory + 2*two_memory < memory_ ? false : true);

    if (print_) {
        fprintf(outfile, "  Using %s Algorithm.\n\n", (is_disk_ ? "Disk" : "Core"));
    }
    if (print_ > 1) {
        fprintf(outfile, "    ----------------------------------------------\n");
        fprintf(outfile, "     Quantity         Doubles              MiB   \n");
        fprintf(outfile, "    ----------------------------------------------\n");
        fprintf(outfile, "     %-6s     %15ld  %15ld\n", "Memory", memory_, (ULI) (memory_*8.0 / 1.0E6));
        fprintf(outfile, "     %-6s     %15ld  %15ld\n", "(A|B)", two_memory, (ULI) (two_memory*8.0 / 1.0E6));
        fprintf(outfile, "     %-6s     %15ld  %15ld\n", "(A|mn)", three_memory, (ULI) (three_memory*8.0 / 1.0E6));
        fprintf(outfile, "    ----------------------------------------------\n\n");
    }

    if (print_) {
        if (options_.get_str("RI_INTS_IO") == "LOAD") {
            fprintf(outfile, "  Will attempt to load (Q|mn) integrals from File %d.\n\n", unit_);
        } else if (options_.get_str("RI_INTS_IO") == "SAVE") {
            fprintf(outfile, "  Will save (Q|mn) integrals to File %d.\n\n", unit_);
        }
    }

    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
    AO2USO_ = SharedMatrix(pet->aotoso());
}
void DFHF::set_unit(unsigned int unit)
{
    unit_ = unit;
    if (disk_iter_.get()) {
        disk_iter_->set_unit(unit);
    }
}
void DFHF::set_memory(unsigned long int memory)
{
    memory_ = memory;

    // Make a memory decision here
    is_disk_ = false;
    memory_ = Process::environment.get_memory() / 8L;
    memory_ = (unsigned long int) (0.7 * memory_);

    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    // Two is for buffer space in fitting
    is_disk_ = (three_memory + 2*two_memory < memory_ ? false : true);

    if (print_) {
        fprintf(outfile, "  Changed my mind. Using %s Algorithm.\n\n", (is_disk_ ? "Disk" : "Core"));
    }
}

// ==> Initialize Integrals (Pre-iterations) <== //

void DFHF::initialize()
{
    if (is_initialized_) return;
    is_initialized_ = true;

    if (is_jk_) {
        if (is_disk_)
            initialize_JK_disk();
        else
            initialize_JK_core();
    } else {
        if (is_disk_)
            initialize_JK_disk();
        else
            initialize_J_core();
    }
}

void DFHF::initialize_JK_disk()
{
    // Try to load
    if (options_.get_str("RI_INTS_IO") == "LOAD") {
        Jinv_.reset();
        return;
    }

    int nso = primary_->nbf();
    int nshell = primary_->nshell();
    int naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    int ntri = schwarz_->get_nfun_pairs();
    int nshellpairs = schwarz_->get_nshell_pairs();
    int* schwarz_shell_pairs = schwarz_->get_schwarz_shells();
    int* schwarz_fun_pairs = schwarz_->get_schwarz_funs();
    long int* schwarz_shell_pairs_r = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs_r = schwarz_->get_schwarz_funs_reverse();

    // ==> Memory Sizing <== //
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI buffer_memory = memory_ - 2*two_memory; // Two is for buffer space in fitting

    //fprintf(outfile, "Buffer memory = %ld words\n", buffer_memory);

    //fprintf(outfile,"Schwarz Shell Pairs:\n");
    //for (int MN = 0; MN < nshellpairs; MN++) {
    //    fprintf(outfile,"  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);
    //}

    //fprintf(outfile,"Schwarz Function Pairs:\n");
    //for (int MN = 0; MN < ntri; MN++) {
    //    fprintf(outfile,"  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);
    //}

    //fprintf(outfile,"Schwarz Reverse Shell Pairs:\n");
    //for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    fprintf(outfile,"  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    //fprintf(outfile,"Schwarz Reverse Function Pairs:\n");
    //for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    fprintf(outfile,"  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell
    boost::shared_ptr<IntVector> MN_mem(new IntVector("Memory per MN pair", nshell * (nshell + 1) / 2));
    int *MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < ntri; mn++) {
        int m = schwarz_fun_pairs[2*mn];
        int n = schwarz_fun_pairs[2*mn + 1];

        int M = primary_->function_to_shell(m);
        int N = primary_->function_to_shell(n);

        MN_memp[M * (M + 1) / 2 + N] += naux;
    }

    //MN_mem->print(outfile);

    // Figure out exactly how much memory per M row
    ULI* M_memp = new ULI[nshell];
    memset(static_cast<void*>(M_memp), '\0', nshell*sizeof(ULI));

    for (int M = 0; M < nshell; M++) {
        for (int N = 0; N <= M; N++) {
            M_memp[M] += MN_memp[M * (M + 1) / 2 + N];
        }
    }

    //fprintf(outfile,"  # Memory per M row #\n\n");
    //for (int M = 0; M < nshell; M++)
    //    fprintf(outfile,"   %3d: %10ld\n", M+1,M_memp[M]);
    //fprintf(outfile,"\n");

    // Find and check the minimum required memory for this problem
    ULI min_mem = naux*(ULI) ntri;
    for (int M = 0; M < nshell; M++) {
        if (min_mem > M_memp[M])
            min_mem = M_memp[M];
    }

    if (min_mem > buffer_memory) {
        std::stringstream message;
        message << "SCF::DF: Disk based algorithm requires 2 (A|B) fitting metrics and an (A|mn) chunk on core." << std::endl;
        message << "         This is 2Q^2 + QNP doubles, where Q is the auxiliary basis size, N is the" << std::endl;
        message << "         primary basis size, and P is the maximum number of functions in a primary shell." << std::endl;
        message << "         For this problem, that is " << ((8L*(min_mem + 2*two_memory))) << " bytes before taxes,";
        message << ((80L*(min_mem + 2*two_memory) / 7L)) << " bytes after taxes. " << std::endl;

        throw PSIEXCEPTION(message.str());
    }

    // ==> Reduced indexing by M <== //

    // Figure out the MN start index per M row
    boost::shared_ptr<IntVector> MN_start(new IntVector("MUNU start per M row", nshell));
    int* MN_startp = MN_start->pointer();

    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (schwarz_shell_pairs[2*MN] == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    boost::shared_ptr<IntVector> mn_start(new IntVector("munu start per M row", nshell));
    int* mn_startp = mn_start->pointer();

    mn_startp[0] = schwarz_fun_pairs[0];
    int m_index = 1;
    for (int mn = 0; mn < ntri; mn++) {
        if (primary_->function_to_shell(schwarz_fun_pairs[2*mn]) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row
    boost::shared_ptr<IntVector> MN_col(new IntVector("MUNU cols per M row", nshell));
    int* MN_colp = MN_col->pointer();

    for (int M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1];
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];

    // Figure out the mn columns per M row
    boost::shared_ptr<IntVector> mn_col(new IntVector("munu cols per M row", nshell));
    int* mn_colp = mn_col->pointer();

    for (int M = 1; M < nshell; M++) {
        mn_colp[M - 1] = mn_startp[M] - mn_startp[M - 1];
    }
    mn_colp[nshell - 1] = ntri - mn_startp[nshell - 1];

    //MN_start->print(outfile);
    //MN_col->print(outfile);
    //mn_start->print(outfile);
    //mn_col->print(outfile);

    // ==> Block indexing <== //
    // Sizing by block
    std::vector<int> MN_start_b;
    std::vector<int> MN_col_b;
    std::vector<int> mn_start_b;
    std::vector<int> mn_col_b;

    // Determine MN and mn block starts
    // also MN and mn block cols
    int nblock = 1;
    ULI current_mem = 0L;
    MN_start_b.push_back(0);
    mn_start_b.push_back(0);
    MN_col_b.push_back(0);
    mn_col_b.push_back(0);
    for (int M = 0; M < nshell; M++) {
        if (current_mem + M_memp[M] > buffer_memory) {
            MN_start_b.push_back(MN_startp[M]);
            mn_start_b.push_back(mn_startp[M]);
            MN_col_b.push_back(0);
            mn_col_b.push_back(0);
            nblock++;
            current_mem = 0L;
        }
        MN_col_b[nblock - 1] += MN_colp[M];
        mn_col_b[nblock - 1] += mn_colp[M];
        current_mem += M_memp[M];
    }

    //fprintf(outfile,"Block, MN start, MN cols, mn start, mn cols\n");
    //for (int block = 0; block < nblock; block++) {
    //    fprintf(outfile,"  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block], mn_col_b[block]);
    //}
    //fflush(outfile);

    // Full sizing not required any longer
    MN_mem.reset();
    MN_start.reset();
    MN_col.reset();
    mn_start.reset();
    mn_col.reset();
    delete[] M_memp;

    // ==> Buffer allocation <== //
    int max_cols = 0;
    for (int block = 0; block < nblock; block++) {
        if (max_cols < mn_col_b[block])
            max_cols = mn_col_b[block];
    }

    // Primary buffer
    Qmn_ = SharedMatrix(new Matrix("(Q|mn) (Disk Chunk)", naux, max_cols));
    // Fitting buffer
    SharedMatrix Amn (new Matrix("(Q|mn) (Buffer)",naux,naux));
    double** Qmnp = Qmn_->pointer();
    double** Amnp = Amn->pointer();

    // ==> Prestripe/Jinv <== //
    psio_->open(unit_,PSIO_OPEN_NEW);
    boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Dispatch the prestripe
    timer_on("(Q|mn) Prestripe");
    aio->zero_disk(unit_,"(Q|mn) Integrals",naux,ntri);

    // Form the J symmetric inverse
    timer_on("(Q|A)^-1/2");
    if (options_.get_str("FITTING_TYPE") == "EIG") {
        Jinv_->form_eig_inverse();
    } else {
        throw PSIEXCEPTION("Fitting Metric type is not implemented.");
    }
    double** Jinvp = Jinv_->get_metric()->pointer();
    timer_off("(Q|A)^-1/2");

    // Synch up
    aio->synchronize();
    timer_off("(Q|mn) Prestripe");

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif

    // ==> ERI initialization <== //
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        if (omega_ > 0.0) {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        } else {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        }
        buffer[Q] = eri[Q]->buffer();
    }

    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block];
        int mn_start_val = mn_start_b[block];
        int MN_col_val = MN_col_b[block];
        int mn_col_val = mn_col_b[block];

        // ==> (A|mn) integrals <== //
        timer_on("(A|mn)");
        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[2*MUNU + 0];
            int NU = schwarz_shell_pairs[2*MUNU + 1];
            int nummu = primary_->shell(MU)->nfunction();
            int numnu = primary_->shell(NU)->nfunction();
            int mu = primary_->shell(MU)->function_index();
            int nu = primary_->shell(NU)->function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P)->nfunction();
                int p = auxiliary_->shell(P)->function_index();
                eri[rank]->compute_shell(P,0,MU,NU);
                for (int dm = 0; dm < nummu; dm++) {
                    int omu = mu + dm;
                    for (int dn = 0; dn < numnu;  dn++) {
                        int onu = nu + dn;
                        if (omu >= onu && schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] >= 0) {
                            int delta = schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] - mn_start_val;
                            for (int dp = 0; dp < nump; dp ++) {
                                int op = p + dp;
                                Qmnp[op][delta] = buffer[rank][dp*nummu*numnu + dm*numnu + dn];
                            }
                        }
                    }
                }
            }
        }
        timer_off("(A|mn)");

        // ==> (Q|mn) fitting <== //
        timer_on("(Q|mn)");
        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val)
                cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }
        timer_off("(Q|mn)");

        // ==> Disk striping <== //
        psio_address addr;
        timer_on("(Q|mn) Write");
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"(Q|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }
        timer_off("(Q|mn) Write");
    }

    // ==> Close out <== //
    psio_->close(unit_,1);
    Qmn_.reset();
    delete[] eri;
}

void DFHF::initialize_JK_core()
{

    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = ((ULI)primary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif
    int rank = 0;

    Qmn_ = SharedMatrix(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Try to load
    if (options_.get_str("RI_INTS_IO") == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_, 1);
        Jinv_.reset();
        return;
    }

    //Get a TEI for each thread
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        if (omega_ > 0.0) {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        } else {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        }
        buffer[Q] = eri[Q]->buffer();
    }

    long int* schwarz_shell_pairs = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs = schwarz_->get_schwarz_funs_reverse();
    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    int index;
    //The integrals (A|mn)
    timer_on("(A|mn)");
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
            //fprintf(outfile,"  Thread %d doing MU = %d",rank,MU); fflush(outfile);
        #endif
        nummu = primary_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell)->nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell)->function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    timer_off("(A|mn)");

//    Qmn_->print();

    delete []buffer;
    delete []eri;

    timer_on("(A|B)^-1/2");
    if (options_.get_str("FITTING_TYPE") == "EIG") {
        Jinv_->form_eig_inverse();
    } else {
        throw PSIEXCEPTION("Fitting Metric type is not implemented.");
    }
    double** Jinvp = Jinv_->get_metric()->pointer();
    timer_off("(A|B)^-1/2");

    ULI max_cols = (memory_-three_memory-two_memory) / auxiliary_->nbf();
    if (max_cols < 1)
        max_cols = 1;
    if (max_cols > ntri)
        max_cols = ntri;
    SharedMatrix temp(new Matrix("Qmn buffer", auxiliary_->nbf(), max_cols));
    double** tempp = temp->pointer();

    int nblocks = ntri / max_cols;
    if ((ULI)nblocks*max_cols != ntri) nblocks++;

    int ncol = 0;
    int col = 0;
    timer_on("(Q|mn)");
    for (int block = 0; block < nblocks; block++) {
        ncol = max_cols;
        if (col + ncol > ntri)
            ncol = ntri - col;

        C_DGEMM('N','N',auxiliary_->nbf(), ncol, auxiliary_->nbf(), 1.0,
            Jinvp[0], auxiliary_->nbf(), &Qmnp[0][col], ntri, 0.0,
            tempp[0], max_cols);

        for (int Q = 0; Q < auxiliary_->nbf(); Q++) {
            C_DCOPY(ncol, tempp[Q], 1, &Qmnp[Q][col], 1);
        }

        col += ncol;
    }
    timer_off("(Q|mn)");

    Jinv_.reset();

    //Qmn_->print();

    // Save if needed
    if (options_.get_str("RI_INTS_IO") == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_, 1);
    }
}
void DFHF::initialize_J_core()
{
    // TODO respect pivoting
    timer_on("chol(A|B)");
    Jinv_->form_cholesky_factor();
    timer_off("chol(A|B)");

    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = (ULI)primary_->nbf()*ntri;
    ULI two_memory = (ULI)auxiliary_->nbf()*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif
    int rank = 0;

    Qmn_ = SharedMatrix(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Save if needed
    if (options_.get_str("RI_INTS_IO") == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(A|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_, 1);
        return;
    }

    //Get a TEI for each thread
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        if (omega_ > 0.0) {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        } else {
            eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        }
        buffer[Q] = eri[Q]->buffer();
    }

    long int* schwarz_shell_pairs = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs = schwarz_->get_schwarz_funs_reverse();
    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    int index;
    //The integrals (A|mn)
    timer_on("(A|mn)");
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
            //fprintf(outfile,"  Thread %d doing MU = %d",rank,MU); fflush(outfile);
        #endif
        nummu = primary_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell)->nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell)->function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    timer_off("(A|mn)");

    delete []buffer;
    delete []eri;

    // No transformation!

    // Save if needed
    if (options_.get_str("RI_INTS_IO") == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "(A|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_, 1);
        return;
    }
}

// ==> USO/AO Transformations (each iteration) <== //

void DFHF::USO2AO()
{
    if (Da_->nirrep() == 1) {
        Da_ao_ = Da_;
        Db_ao_ = Db_;
        Ca_ao_ = Ca_;
        Cb_ao_ = Cb_;
        Ja_ao_ = Ja_;
        Ka_ao_ = Ka_;
        Kb_ao_ = Kb_;
        nalpha_ = nalphapi_[0];
        nbeta_ = 0;
        if (!restricted_)
            nbeta_ = nbetapi_[0];
        return;
    }

    timer_on("USO2AO");
    // Sizing
    nalpha_ = 0;
    nbeta_ = 0;
    for (int h = 0; h < Da_->nirrep(); h++) {
        nalpha_ += nalphapi_[h];
        if (!restricted_)
            nbeta_ += nbetapi_[h];
    }

    int nao = primary_->nbf();
    int nmo = 0;
    if (is_jk_)
        for (int h = 0; h < Ca_->nirrep(); h++) nmo += Ca_->colspi()[h];

    // Allocation
    Da_ao_ = SharedMatrix(new Matrix("D Alpha (AO Basis)", nao, nao));
    Ja_ao_ = SharedMatrix(new Matrix("J (AO Basis)", nao, nao));
    if (!restricted_)
        Db_ao_ = SharedMatrix(new Matrix("D Beta (AO Basis)", nao, nao));

    if (is_jk_) {
        Ca_ao_ = SharedMatrix(new Matrix("C Alpha (AO Basis)", nao, nmo));
        Ka_ao_ = SharedMatrix(new Matrix("K Alpha (AO Basis)", nao, nao));
        if (!restricted_) {
            Cb_ao_ = SharedMatrix(new Matrix("C Beta (AO Basis)", nao, nmo));
            Kb_ao_ = SharedMatrix(new Matrix("K Beta (AO Basis)", nao, nao));
        }
    }

    // Filling
    double** D_ao = Da_ao_->pointer();
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nso = Da_->colspi()[h];
        if (nso == 0) continue;

        double** D_so = Da_->pointer(h);
        double** U = AO2USO_->pointer(h);
        double** Temp = block_matrix(nso,nao);

        C_DGEMM('N','T',nso,nao,nso,1.0,D_so[0],nso,U[0],nso,0.0,Temp[0],nao);
        C_DGEMM('N','N',nao,nao,nso,1.0,U[0],nso,Temp[0],nao,1.0,D_ao[0],nao);

        free_block(Temp);
    }

    if (!restricted_) {
        double** D_ao = Db_ao_->pointer();
        for (int h = 0; h < Da_->nirrep(); h++) {
            int nso = Da_->colspi()[h];
            if (nso == 0) continue;

            double** D_so = Db_->pointer(h);
            double** U = AO2USO_->pointer(h);
            double** Temp = block_matrix(nso,nao);

            C_DGEMM('N','T',nso,nao,nso,1.0,D_so[0],nso,U[0],nso,0.0,Temp[0],nao);
            C_DGEMM('N','N',nao,nao,nso,1.0,U[0],nso,Temp[0],nao,1.0,D_ao[0],nao);

            free_block(Temp);
        }
    }

    if (is_jk_) {
        double** C_ao = Ca_ao_->pointer();
        int counter = 0;
        for (int h = 0; h < Da_->nirrep(); h++) {
            int nso = Da_->colspi()[h];
            int nmopi = Ca_->colspi()[h];
            int nalpha = nalphapi_[h];
            if (nso == 0 || nmopi == 0 || nalpha == 0) continue;

            double** C_so = Ca_->pointer(h);
            double** U = AO2USO_->pointer(h);

            C_DGEMM('N','N',nao,nalpha,nso,1.0,U[0],nso,C_so[0],nmopi,0.0,&C_ao[0][counter],nmo);

            counter += nalpha;
        }
        if (!restricted_) {
            double** C_ao = Cb_ao_->pointer();
            counter = 0;
            for (int h = 0; h < Da_->nirrep(); h++) {
                int nso = Da_->colspi()[h];
                int nmopi = Ca_->colspi()[h];
                int nbeta = nbetapi_[h];
                if (nso == 0 || nmopi == 0 || nbeta == 0) continue;

                double** C_so = Cb_->pointer(h);
                double** U = AO2USO_->pointer(h);

                C_DGEMM('N','N',nao,nbeta,nso,1.0,U[0],nso,C_so[0],nmopi,0.0,&C_ao[0][counter],nmo);

                counter += nbeta;
            }
        }
    }
    timer_off("USO2AO");
}
void DFHF::AO2USO()
{
    if (Da_->nirrep() == 1) return;

    timer_on("AO2USO");
    int nao = Ja_ao_->colspi()[0];

    double** J_ao = Ja_ao_->pointer();
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nso = Da_->colspi()[h];
        if (nso == 0) continue;

        double** J_so = Ja_->pointer(h);
        double** U = AO2USO_->pointer(h);
        double** Temp = block_matrix(nso,nao);

        C_DGEMM('T','N',nso,nao,nao,1.0,U[0],nso,J_ao[0],nao,0.0,Temp[0],nao);
        C_DGEMM('N','N',nso,nso,nao,1.0,Temp[0],nao,U[0],nso,0.0,J_so[0],nso);

        free_block(Temp);
    }

    if (!is_jk_) {
        timer_off("AO2USO");
        return;
    }

    double** K_ao = Ka_ao_->pointer();
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nso = Da_->colspi()[h];
        if (nso == 0) continue;

        double** K_so = Ka_->pointer(h);
        double** U = AO2USO_->pointer(h);
        double** Temp = block_matrix(nso,nao);

        C_DGEMM('T','N',nso,nao,nao,1.0,U[0],nso,K_ao[0],nao,0.0,Temp[0],nao);
        C_DGEMM('N','N',nso,nso,nao,1.0,Temp[0],nao,U[0],nso,0.0,K_so[0],nso);

        free_block(Temp);
    }

    if (restricted_) {
        timer_off("AO2USO");
        return;
    }

    K_ao = Kb_ao_->pointer();
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nso = Da_->colspi()[h];
        if (nso == 0) continue;

        double** K_so = Kb_->pointer(h);
        double** U = AO2USO_->pointer(h);
        double** Temp = block_matrix(nso,nao);

        C_DGEMM('T','N',nso,nao,nao,1.0,U[0],nso,K_ao[0],nao,0.0,Temp[0],nao);
        C_DGEMM('N','N',nso,nso,nao,1.0,Temp[0],nao,U[0],nso,0.0,K_so[0],nso);

        free_block(Temp);
    }
    timer_off("AO2USO");
}

// ==> Block computations of J and K, JK algorithms  <== //

void DFHF::compute_JK_block_J(double** Qmnp, int nrows, int max_rows)
{
    timer_on("Form J");
    int ntri = schwarz_->get_nfun_pairs();
    C_DGEMV('N', nrows, ntri, 1.0, Qmnp[0], ntri, Dtri_, 1, 0.0, dQ_, 1);
    C_DGEMV('T', nrows, ntri, 1.0, Qmnp[0], ntri, dQ_, 1, 1.0, Jtri_, 1);
    timer_off("Form J");
}
void DFHF::compute_JK_block_K(double** Qmnp, int nrows, int max_rows, bool is_alpha)
{
    timer_on("Form K");
    int nbf = primary_->nbf();
    int nalpha = nalpha_;
    int nelec;
    if (is_alpha)
        nelec = nalpha_;
    else
        nelec = nbeta_;
    int ntri = schwarz_->get_nfun_pairs();

    double** Cp;
    double** Kp;
    if (is_alpha) {
        Cp = Ca_ao_->pointer();
        Kp = Ka_ao_->pointer();
    } else {
        Cp = Cb_ao_->pointer();
        Kp = Kb_ao_->pointer();
    }

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    int rank = 0;

    #ifdef HAVE_MKL
        int mkl_nthread = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    int m, n , ij, index;
    timer_on("(Q|mi)");
    #pragma omp parallel for private (m, n , ij, index, rank) schedule (dynamic)
    for (m = 0; m<nbf; m++) {

        rank = 0;
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        int n, ij;
        for (index = 0; index<index_sizes_[m]; index++) {
            ij = m_ij_indices_[m][index];
            n = n_indices_[m][index];

            C_DCOPY(nrows,&Qmnp[0][ij],ntri,&QS_[rank][0][index],nbf);
            C_DCOPY(nelec,Cp[n],1,&Ctemp_[rank][0][index],nbf);
        }

        C_DGEMM('N','T',nelec,nrows,index_sizes_[m],1.0,Ctemp_[rank][0],nbf,QS_[rank][0],nbf, 0.0, Ep_[m], nrows);
    }
    timer_off("(Q|mi)");

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_nthread);
    #endif

    timer_on("K");
    C_DGEMM('N','T',nbf,nbf,nrows*nelec,1.0,Ep_[0],max_rows*nalpha,Ep_[0],max_rows*nalpha,1.0,Kp[0], nbf);
    timer_off("K");

    timer_off("Form K");
}

// ==> One shot computation of J, J core algorithm <== //

void DFHF::compute_J_core()
{
    timer_on("Form J");
    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();
    int ntri = schwarz_->get_nfun_pairs();

    double** Qmnp = Qmn_->pointer();

    SharedMatrix Dt(new Matrix("D Total",nbf,nbf));
    double** Dtp = Dt->pointer();
    Dt->copy(Da_ao_);
    Dt->add(Db_ao_);
    double* Dtri = new double[ntri];
    double* Jtri = new double[ntri];
    int* schwarz_funs = schwarz_->get_schwarz_funs();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        double perm = (mu == nu ? 1.0 : 2.0);
        Dtri[munu] = perm * Dtp[mu][nu];
    }
    Dt.reset();

    double *dQ = new double[naux];
    C_DGEMV('N', naux, ntri, 1.0, Qmnp[0], ntri, Dtri, 1, 0.0, dQ, 1);
    // TODO pivoting
    C_DPOTRS('L', naux, 1, Jinv_->get_metric()->pointer()[0], naux, dQ, naux);
    C_DGEMV('T', naux, ntri, 1.0, Qmnp[0], ntri, dQ, 1, 0.0, Jtri, 1);

    double** Jp = Ja_ao_->pointer();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        Jp[mu][nu] += Jtri[munu];
        if (mu != nu)
            Jp[nu][mu] += Jtri[munu];
    }

    delete[] dQ;
    delete[] Jtri;
    delete[] Dtri;
    timer_off("Form J");
}

// ==> Drivers <== //

void DFHF::form_J_DF()
{
    initialize();

    USO2AO();

    if (is_disk_) {
        // TODO
        //do {
        //    SharedMatrix Q = disk_iter_->next_block();
        //    int current_rows = disk_iter_->current_rows();
        //    double** Qmnp = Q->pointer();
        //    compute_JK_block_J(Qmnp, current_rows, max_rows);
        //} while (!disk_iter_->finished());
    } else {
        compute_J_core();
    }

    AO2USO();
}
void DFHF::form_JK_DF()
{
    // Initialize if needed (build Qmn)
    initialize();

    // Back D and C from USO -> SO
    USO2AO();

    // Standard indexing
    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();
    int nalpha = nalpha_;
    int ntri = schwarz_->get_nfun_pairs();

    // Number of threads (needed for allocation/memory)
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    // Determine max rows, nblocks, possible disk iterator
    ULI max_rows;
    int nblocks;
    if (is_disk_) {
        if (!is_disk_initialized_) {
            is_disk_initialized_ = true;

            // Static disk overhead, CTemp * nthread, Jtri, Dtri, indexing
            ULI overhead = nthread*(ULI)nalpha*nbf + nbf*(ULI)nbf + nbf + 2L* (ULI)ntri;
            // Dynamic disk overhead, (A|mn), (A|mi), Q_S * nthread
            ULI per_row = nalpha*(ULI)nbf + nthread*(ULI)nbf + (ULI)ntri;

            // How many rows at a time, on disk?
            max_rows = (memory_ - overhead) / per_row;
            if (max_rows > naux)
                max_rows = naux;
            if (max_rows < 1)
                max_rows = 1;

            // Now split things in two
            max_rows = max_rows / 2 + (((max_rows % 2) == 0) ? 0 : 1);

            // Determine number of blocks (gimp shouldn't hurt much, K is fine-grained)
             nblocks = naux / max_rows;
            if (nblocks * max_rows != naux)
                nblocks++;

            // Build the disk iterator
            disk_iter_ = boost::shared_ptr<DFHFDiskIterator>(new DFHFDiskIterator(psio_, ntri, naux, max_rows));
            disk_iter_->set_unit(unit_);
         } else {
            max_rows = disk_iter_->max_rows();
            nblocks = disk_iter_->nblock();
         }
    } else {
        // Static core overhead, (A|mn), CTemp * nthread, Jtri, Dtri, indexing
        ULI overhead = nthread*(ULI)nalpha*nbf + ntri*(ULI)naux + nbf*(ULI)nbf + nbf + 2L * (ULI)ntri;
        // Dynamic core overhead, (A|mi),  Q_S * nthread
        ULI per_row = nalpha*(ULI)nbf + nthread*(ULI)nbf;

        // How many rows at a time, on core?
        max_rows = (memory_ - overhead) / per_row;
        if (max_rows > naux)
            max_rows = naux;
        if (max_rows < 1)
            max_rows = 1;

        // Determine number of blocks (gimp shouldn't hurt much, K is fine-grained)
        nblocks = naux / max_rows;
        if (nblocks * max_rows != naux)
            nblocks++;
    }

    // JK J allocation
    Dtri_ = new double[ntri];
    Jtri_ = new double[ntri];
    dQ_ = new double[naux];
    memset(static_cast<void*>(Jtri_), '\0', ntri*sizeof(double));

    // Build D triangular matrix
    double** Dap = Da_ao_->pointer();
    double** Dbp = restricted_ ? Da_ao_->pointer() : Db_ao_->pointer();
    int* schwarz_funs = schwarz_->get_schwarz_funs();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        double perm = (mu == nu ? 1.0 : 2.0);
        Dtri_[munu] = perm * (Dap[mu][nu] + Dbp[mu][nu]);
    }

    // JK K allocation
    Ep_ = block_matrix(nbf,nalpha*(ULI)max_rows);
    //QS temp matrix for DGEMM
    QS_ = new double**[nthread];
    for (int T = 0; T < nthread; T++)
        QS_[T] = block_matrix(max_rows,nbf);
    // Temp matrix for sparse DGEMM if sieve exists
    Ctemp_ = new double**[nthread];
    for (int T = 0; T < nthread; T++)
        Ctemp_[T] = block_matrix(nalpha,nbf);
    // Index array for non-canonical ordering of mn
    m_ij_indices_ = init_int_matrix(nbf,nbf);
    // Index array of n for given m (in order of above)
    n_indices_ = init_int_matrix(nbf,nbf);
    // sizes of above for schwarz sieve
    index_sizes_ = init_int_array(nbf);

    for (int ij = 0; ij<ntri; ij++) {
        int m = schwarz_funs[2*ij];
        int n = schwarz_funs[2*ij+1];

        m_ij_indices_[m][index_sizes_[m]] = ij;
        n_indices_[m][index_sizes_[m]] = n;
        index_sizes_[m]++;
        if (m != n){
            m_ij_indices_[n][index_sizes_[n]] = ij;
            n_indices_[n][index_sizes_[n]] = m;
            index_sizes_[n]++;
        }
    }

    // Block formation of J/K
    if (is_disk_) {
        int counter = 0;
        do {
            double** Qmnp = disk_iter_->next_block();
            int current_rows = disk_iter_->current_rows();
            compute_JK_block_J(Qmnp, current_rows, max_rows);
            compute_JK_block_K(Qmnp, current_rows, max_rows, true);
            if (!restricted_)
                compute_JK_block_K(Qmnp, current_rows, max_rows, false);
        } while (!disk_iter_->finished());
    } else {
        double** Qmnp = Qmn_->pointer();
        for (int block = 0; block < nblocks; block++) {
            int current_rows = max_rows;
            if (block == nblocks - 1)
                current_rows = naux - (block*max_rows);
            compute_JK_block_J(&Qmnp[block*(ULI)max_rows], current_rows, max_rows);
            compute_JK_block_K(&Qmnp[block*(ULI)max_rows], current_rows, max_rows, true);
            if (!restricted_)
                compute_JK_block_K(&Qmnp[block*(ULI)max_rows], current_rows, max_rows, false);
        }
    }

    // JK J copy out
    double** Jp = Ja_ao_->pointer();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        Jp[mu][nu] += Jtri_[munu];
        if (mu != nu)
            Jp[nu][mu] += Jtri_[munu];
    }

    // JK J frees
    delete[] Dtri_;
    delete[] Jtri_;
    delete[] dQ_;

    // JK K frees
    free_block(Ep_);
    for (int thread = 0; thread < nthread; thread++) {
        free_block(QS_[thread]);
        free_block(Ctemp_[thread]);
    }
    delete[] QS_;
    delete[] Ctemp_;
    free(m_ij_indices_[0]);
    free(m_ij_indices_);
    free(n_indices_[0]);
    free(n_indices_);
    free(index_sizes_);

    // Move J and K from AO -> USO
    AO2USO();
}

DFHFDiskIterator::DFHFDiskIterator(boost::shared_ptr<PSIO> psio, int ntri, int naux, int max_rows) :
    psio_(psio), ntri_(ntri), naux_(naux), max_rows_(max_rows), unit_(PSIF_DFSCF_BJ)
{
    common_init();
}
DFHFDiskIterator::~DFHFDiskIterator()
{
    // Close and keep the disk
    psio_->close(unit_,1);
}
void DFHFDiskIterator::common_init()
{
    // Open the disk
    psio_->open(unit_,PSIO_OPEN_OLD);

    // Sizing
    nblocks_ =  naux_ / max_rows_;
    int ngimp = naux_ % max_rows_;
    if (ngimp != 0) nblocks_++;

    if (nblocks_ < 2) throw SanityCheckError("SCF::DFHFDiskIterator: Nblocks must be > 1, or you'd be on core",__FILE__,__LINE__);

    block_sizes_.resize(nblocks_);
    block_starts_.resize(nblocks_);

    // TODO: ungimp later
    for (int i = 0; i < nblocks_; i++) {
        block_starts_[i] = i * max_rows_;
        block_sizes_[i] = max_rows_;
        if (i == nblocks_ - 1) {
            block_sizes_[i] = naux_ - i * max_rows_;
        }
    }

    // Initial order in the cyclic structure
    for (int i = 0; i < nblocks_; i++)
        blocks_.push_back(i);

    //fprintf(outfile,"  DFHFDiskIterator::common_init Debug:\n\n");
    //fprintf(outfile,"  Max rows = %d, Naux = %d, nblocks = %d\n\n", max_rows_, naux_, nblocks_);
    //fprintf(outfile,"  Block Starts/Block Sizes\n");
    //for (int i = 0; i < nblocks_; i++) {
    //    fprintf(outfile,"   Block %3d, %3d/%3d\n", i, block_starts_[i], block_sizes_[i]);
    //}
    //fprintf(outfile,"  Block Order\n");
    //for (int i = 0; i < nblocks_; i++) {
    //    fprintf(outfile,"  Local Index = %3d, Global Index = %3d\n", i, blocks_[i]);
    //}
    //fflush(outfile);

    // Build the blocks
    A_ = SharedMatrix(new Matrix("(Q|mn) Block A", max_rows_, ntri_));
    B_ = SharedMatrix(new Matrix("(Q|mn) Block B", max_rows_, ntri_));

    // Build the AIO object
    aio_ = boost::shared_ptr<AIOHandler>(new AIOHandler(psio_));

    // Read the first two bits in, so they are ready
    read(A_, block_starts_[0], block_sizes_[0]);
    synchronize();
    read(B_, block_starts_[1], block_sizes_[1]);
    synchronize();

    iteration_ = 0;
}
double** DFHFDiskIterator::next_block()
{
    synchronize();
    //fprintf(outfile, "  DFHFDiskIterator::next_block Debug, Block %d\n", iteration_);

    // Synchronize the AIOHandler, if needed
    // Not needed for the first two iterations of the cycle
    if (iteration_ > 1) {
        //synchronize();
        //fprintf(outfile, "  Flushing AIO\n");
    }

    // If not first or last block, there's reading to be done
    if (iteration_ > 0 && iteration_ < nblocks_ - 1) {
        // The address in the cycle is two greater than the iteration number
        int index = blocks_[iteration_ + 1];
        if (iteration_ % 2 == 1) {
            // Odd iteration, read into A
            read(A_, block_starts_[index], block_sizes_[index]);
            //fprintf(outfile,"  Reading index %d, start = %3d, size = %3d into A\n", index, \
                block_starts_[index], block_sizes_[index]);
        } else {
            // Even iteration, read into B
            read(B_, block_starts_[index], block_sizes_[index]);
            //fprintf(outfile,"  Reading index %d, start = %3d, size = %3d into B\n", index, \
                block_starts_[index], block_sizes_[index]);
        }
    }

    // Determine current rows
    int index = blocks_[iteration_];
    current_rows_ = block_sizes_[index];

    SharedMatrix val;

    //fprintf(outfile, "  Current rows = %d\n", current_rows_);
    fflush(outfile);
    // Return a block
    // Odd iteration, return B
    if ((iteration_++) % 2 == 1) {
        return B_->pointer();
        //fprintf(outfile, "  Pointing to B, which is global block %d\n", index);
    // Even iteration, return A
    } else {
        return A_->pointer();
        //fprintf(outfile, "  Pointing to A, which is global block %d\n", index);
    }
}
void DFHFDiskIterator::read(SharedMatrix A, int start, int rows)
{
    // Get the block pointer
    double** Ap = A->pointer();

    // Get an address
    psio_address addr = psio_get_address(PSIO_ZERO, sizeof(double)*start*ntri_);

    // Post the read
    psio_->read(unit_, "(Q|mn) Integrals", (char*) Ap[0], sizeof(double)*rows*ntri_, addr, &addr);
    //aio_->read(unit_, "(Q|mn) Integrals", (char*) Ap[0], sizeof(double)*rows*ntri_, addr, &addr);
}
void DFHFDiskIterator::synchronize()
{
    //aio_->synchronize();
}
bool DFHFDiskIterator::finished()
{
    if (iteration_ == nblocks_) {
        reset();
        return true;
    }
    return false;
}
void DFHFDiskIterator::reset()
{
    iteration_ = 0;

    // Reoder the cyclic structure, with the last two blocks first
    int k = blocks_[nblocks_ - 2];
    int counter = 0;
    for (int j = k; j < nblocks_; j++) {
        blocks_[counter++] = j;
    }
    for (int j = 0; j < k; j++) {
        blocks_[counter++] = j;
    }

    // If odd number of blocks, must switch A and B to have A first in next iteration
    if (nblocks_ % 2 == 1)
        boost::swap(A_,B_);

    //fprintf(outfile,"  DFHFDiskIterator::reset Debug: Block Order\n");
    //for (int i = 0; i < nblocks_; i++) {
        //fprintf(outfile,"  Local Index = %3d, Global Index = %3d\n", i, blocks_[i]);
    //}
    //fflush(outfile);
}

}}
