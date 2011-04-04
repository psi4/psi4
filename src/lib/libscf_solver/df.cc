#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
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

DFHF::DFHF(shared_ptr<BasisSet> basis, shared_ptr<PSIO> psio, Options& opt) :
    primary_(basis), psio_(psio), options_(opt)
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
        fprintf(outfile, "   by Rob Parrish\n\n");
    }    

    is_initialized_ = false;
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

    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, primary_->molecule(), "RI_BASIS_SCF");
    parser.reset();

    if (print_) {
        fprintf(outfile, "  ==> Auxiliary Basis: %s <==\n\n", options_.get_str("RI_BASIS_SCF").c_str());
        auxiliary_->print_by_level(outfile, print_);
    }

    timer_on("Schwarz");
    schwarz_ = shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, options_.get_double("SCHWARZ_CUTOFF")));
    timer_off("Schwarz");

    Jinv_ = shared_ptr<FittingMetric>(new FittingMetric(auxiliary_, true));

    // Make a memory decision here
    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = ((ULI)primary_->nbf())*ntri;
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

    if (print_) {
        if (options_.get_str("RI_INTS_IO") == "LOAD") {
            fprintf(outfile, "  Will attempt to load (Q|mn) integrals from File 97.\n\n");    
        } else if (options_.get_str("RI_INTS_IO") == "SAVE") {
            fprintf(outfile, "  Will save (Q|mn) integrals to File 97.\n\n");    
        }
    }

    shared_ptr<IntegralFactory> integral(new IntegralFactory(primary_,primary_,primary_,primary_));
    shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
    AO2USO_ = shared_ptr<Matrix>(pet->aotoso());
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
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        Jinv_.reset();
        return;
    } 

    // Indexing
    int naux = auxiliary_->nbf();
    int ntri = schwarz_->get_nfun_pairs();
    int nshellpairs = schwarz_->get_nshell_pairs();
    int* schwarz_shell_pairs = schwarz_->get_schwarz_shells();
    int* schwarz_fun_pairs = schwarz_->get_schwarz_funs();
    long int* schwarz_shell_pairs_r = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs_r = schwarz_->get_schwarz_funs_reverse();

    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();
    ULI buffer_memory = memory_ - 2*two_memory; // Two is for buffer space in fitting

    std::vector<ULI> memory_per_MU; // Ceiling memory of all functions for given MU in canonical order
    std::vector<int> MU_starts; // First MU shell index in a block
    std::vector<int> NU_starts; // First NU shell index for the given MU (as [MU, 0] is not always cool)
    std::vector<int> schwarz_pair_starts; // MUNU starting indices in the abbreviated schwarz_shell_pairs array
    std::vector<int> schwarz_pair_sizes; // Number of MUNU indices in the abbreviated schwarz_shell_pairs array
    std::vector<int> munu_offsets; // Offset to first munu compound index (reduced) in each block

    // Figure out how much each MU row costs (in triangular, shell sieved form, ceiling value)
    memory_per_MU.resize(primary_->nshell());
    for (int MU = 0; MU < primary_->nshell(); MU++) memory_per_MU[MU] = 0L;
    for (int shellpair = 0; shellpair < nshellpairs; shellpair++) {
        memory_per_MU[schwarz_shell_pairs[2*shellpair]] += (ULI)primary_->shell(schwarz_shell_pairs[2*shellpair + 1])->nfunction();
    }
    for (int MU = 0; MU < primary_->nshell(); MU++) memory_per_MU[MU] *= (ULI)primary_->shell(MU)->nfunction() * naux;

    // Figure out where appropriate blocks start in MU
    MU_starts.push_back(0);
    ULI memory_counter = 0L;
    for (int MU = 0; MU < primary_->nshell(); MU++) {
        if (memory_counter + memory_per_MU[MU] >= buffer_memory) {
            MU_starts.push_back(MU);
            memory_counter = 0L;
        }
        memory_counter += memory_per_MU[MU];
    }

    // How many blocks are there?
    int nblock = MU_starts.size();
    schwarz_pair_starts.resize(nblock);
    schwarz_pair_sizes.resize(nblock);
    munu_offsets.resize(nblock);
    NU_starts.resize(nblock);

    // Determine block metadata
    schwarz_pair_starts[0] = 0;
    NU_starts[0] = schwarz_shell_pairs[1];
    int block_index = 1;
    int munu_size = 0;
    int MUNU_size = 0;
    int max_cols = 0;
    for (int shellpair = 0; shellpair < nshellpairs; shellpair++) {
        int MU = schwarz_shell_pairs[2*shellpair];
        int NU = schwarz_shell_pairs[2*shellpair + 1];

        if (block_index < nblock && MU == MU_starts[block_index]) {
            schwarz_pair_starts[block_index] = shellpair;
            schwarz_pair_sizes[block_index - 1] = MUNU_size;
            NU_starts[block_index] = NU;
            if (munu_size > max_cols);
                max_cols = munu_size;
            block_index++;
            MUNU_size = 0;
            munu_size = 0;
        }

        if (MU == NU)
            munu_size += primary_->shell(MU)->nfunction() * (primary_->shell(NU)->nfunction() + 1) / 2;
        else
            munu_size += primary_->shell(MU)->nfunction() * primary_->shell(NU)->nfunction();

        MUNU_size++;
    }
    schwarz_pair_sizes[nblock - 1] = MUNU_size;

    // Determine munu offsets by block (use the schwarz backmap)
    // This is perhaps a bit excessive, but will always work
    for (int block = 1; block < nblock; block++) {
        int MU = MU_starts[block];
        int NU = NU_starts[block];
        int mu = primary_->shell(MU)->function_index();
        int nu = primary_->shell(NU)->function_index();
        for (int m = 0; m < primary_->shell(MU)->nfunction(); m++) {
            for (int n = 0; n < primary_->shell(NU)->nfunction(); n++) {
                int om = mu + m;
                int on = nu + n;
                if (schwarz_fun_pairs_r[om*(om+1)/2 + on] >= 0){
                    munu_offsets[block] = schwarz_fun_pairs_r[om*(om+1)/2 + on];
                    break;
                }
            }
        }
    }

    // Thread setup
    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif

    //Get a TEI for each thread
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    // Get an eig fitting metric
    if (options_.get_str("FITTING_TYPE") == "EIG") {
        Jinv_->form_eig_inverse();
    } else {
        throw PSIEXCEPTION("Fitting Metric type is not implemented.");
    }
    double** Jinvp = Jinv_->get_metric()->pointer();

    // The (A|mn)/(Q|mn) tensor
    Qmn_ = shared_ptr<Matrix>(new Matrix("(Q|mn) (Disk Chunk)", naux, max_cols));
    // A buffer for fitting purposes
    shared_ptr<Matrix>Amn (new Matrix("(Q|mn) (Buffer)",naux,naux));

    double** Qmnp = Qmn_->pointer();
    double** Amnp = Amn->pointer();

    // Prestripe
    double* prestripe = new double[ntri];
    memset(static_cast<void*>(prestripe), '\0', ntri*sizeof(double));
    psio_address pre_addr = PSIO_ZERO;
    for (int Q = 0; Q < naux; Q++) {
            psio_->write(PSIF_DFSCF_BJ,"(Q|mn) Integrals",(char *) &(prestripe[0]),sizeof(double)*ntri,pre_addr,&pre_addr);
    }
    delete[] prestripe;

    // Loop over blocks
    for (int block = 0; block < nblock; block++) {
        int schwarz_start = schwarz_pair_starts[block];
        int schwarz_sizes = schwarz_pair_sizes[block];
        int offset = munu_offsets[block];
        int current_cols = 0;

        // COmpute current cols (prevents required reduction later)
        for (int MUNU = 0; MUNU < schwarz_sizes; MUNU++) {

            int MU = schwarz_shell_pairs[2*(MUNU + schwarz_start) + 0];
            int NU = schwarz_shell_pairs[2*(MUNU + schwarz_start) + 1];
            int nummu = primary_->shell(MU)->nfunction();
            int numnu = primary_->shell(NU)->nfunction();
            int mu = primary_->shell(MU)->function_index();
            int nu = primary_->shell(NU)->function_index();
            for (int dm = 0; dm < nummu; dm++) {
                int omu = mu + dm;
                for (int dn = 0; dn < numnu;  dn++) {
                    int onu = nu + dm;
                    if (omu >= onu && schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] >= 0) {
                            current_cols++;
                    }
                }
            }
        }

        // Form (A|mn) integrals
        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = 0; MUNU < schwarz_sizes; MUNU++) {
            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[2*(MUNU + schwarz_start) + 0];
            int NU = schwarz_shell_pairs[2*(MUNU + schwarz_start) + 1];
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
                        int onu = nu + dm;
                        if (omu >= onu && schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] >= 0) {
                            for (int dp = 0; dp < nump; dp ++) {
                                int op = p + dp;
                                Qmnp[op][schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] - offset] = buffer[rank][op*nummu*numnu + omu*numnu + onu];
                            }
                        }
                    }
                }
            }
        }

        // Apply Fitting
        for (int mn = 0; mn < current_cols; mn+=naux) {
            int cols = naux;
            if (mn + naux >= current_cols)
                cols = current_cols - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,Qmnp[Q],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        // Stripe to disk
        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + offset)*sizeof(double));
            psio_->write(PSIF_DFSCF_BJ,"(Q|mn) Integrals",(char*)Qmnp[Q],current_cols*sizeof(double),addr,&addr);
        }
    }
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

    Qmn_ = shared_ptr<Matrix>(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Try to load
    if (options_.get_str("RI_INTS_IO") == "LOAD") {
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(PSIF_DFSCF_BJ, 1);
        Jinv_.reset();
        return;
    } 
    
    //Get a TEI for each thread
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
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
    shared_ptr<Matrix> temp(new Matrix("Qmn buffer", auxiliary_->nbf(), max_cols));
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

    // Save if needed
    if (options_.get_str("RI_INTS_IO") == "SAVE") {
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        psio_->write_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(PSIF_DFSCF_BJ, 1);
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

    Qmn_ = shared_ptr<Matrix>(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Save if needed
    if (options_.get_str("RI_INTS_IO") == "LOAD") {
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFSCF_BJ, "(A|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(PSIF_DFSCF_BJ, 1);
        return;
    } 
    
    //Get a TEI for each thread
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
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
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        psio_->write_entry(PSIF_DFSCF_BJ, "(A|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(PSIF_DFSCF_BJ, 1);
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
    Da_ao_ = shared_ptr<Matrix>(new Matrix("D Alpha (AO Basis)", nao, nao));
    Ja_ao_ = shared_ptr<Matrix>(new Matrix("J (AO Basis)", nao, nao));
    if (!restricted_) 
        Db_ao_ = shared_ptr<Matrix>(new Matrix("D Beta (AO Basis)", nao, nao));

    if (is_jk_) {
        Ca_ao_ = shared_ptr<Matrix>(new Matrix("C Alpha (AO Basis)", nao, nmo));
        Ka_ao_ = shared_ptr<Matrix>(new Matrix("K Alpha (AO Basis)", nao, nao));
        if (!restricted_) {
            Cb_ao_ = shared_ptr<Matrix>(new Matrix("C Beta (AO Basis)", nao, nmo));
            Kb_ao_ = shared_ptr<Matrix>(new Matrix("K Beta (AO Basis)", nao, nao));
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

    shared_ptr<Matrix> Dt(new Matrix("D Total",nbf,nbf));
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
        //do {
        //    shared_ptr<Matrix> Q = disk_iter_->next_block();
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

    // Determine max rows
    ULI overhead;
    ULI per_row;
    if (is_disk_) {
        overhead = nthread*(ULI)nalpha*nbf + nbf*(ULI)nbf + nbf + 2L* (ULI)ntri;
        per_row = nalpha*(ULI)nbf + nthread*(ULI)nbf + (ULI)ntri;
    } else {
        overhead = nthread*(ULI)nalpha*nbf + ntri*(ULI)naux + nbf*(ULI)nbf + nbf + 2L * (ULI)ntri;
        per_row = nalpha*(ULI)nbf + nthread*(ULI)nbf;
    }
    int max_rows = (memory_ - overhead) / per_row;
    if (max_rows > naux)
        max_rows = naux;
    if (max_rows < 1)
        max_rows = 1;

    // Get's a little weird, for disk, we divide into two blocks within the big block for AIO
    if (is_disk_) {
        max_rows /= 2;
        if (max_rows < 1)
            max_rows = 1;
        
    }

    // Determine number of blocks (gimp shouldn't hurt much, K is fine-grained)
    int nblocks = naux / max_rows;
    if (nblocks * max_rows != naux)
        nblocks++;

    printf("Doing DF-SCF iteration with %d blocks\n", nblocks);

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
        do {
            shared_ptr<Matrix> Q = disk_iter_->next_block();
            int current_rows = disk_iter_->current_rows();
            double** Qmnp = Q->pointer();
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

DFHFDiskIterator::DFHFDiskIterator(shared_ptr<PSIO> psio, int ntri, int naux, int max_rows) :
    psio_(psio), ntri_(ntri), naux_(naux), max_rows_(max_rows)
{
    common_init();
}
DFHFDiskIterator::~DFHFDiskIterator()
{
    // Close and keep the disk
    psio_->close(PSIF_DFSCF_BJ,1);
}
void DFHFDiskIterator::common_init()
{
    // Open the disk
    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);

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
            block_sizes_[i] = naux_ - i* max_rows_;
        }
    }

    // Initial order in the cyclic structure
    for (int i = 0; i < nblocks_; i++)
        blocks_.push_back(i);  


    // Build the blocks
    A_ = shared_ptr<Matrix>(new Matrix("(Q|mn) Block A", max_rows_, ntri_));
    B_ = shared_ptr<Matrix>(new Matrix("(Q|mn) Block B", max_rows_, ntri_));
        
    // Read the first two bits in, so they are ready
    read(A_, block_starts_[0], block_sizes_[0]);
    aio_->synchronize();
    read(B_, block_starts_[1], block_sizes_[1]);
    aio_->synchronize();
}
boost::shared_ptr<Matrix> DFHFDiskIterator::next_block()
{
    iteration_++;    

    // Synchronize the AIOHandler, if needed
    // Not needed for the first two iterations of the cycle
    if (iteration_ > 1) {
        aio_->synchronize();
    }

    // If not first or last block, there's reading to be done
    if (iteration_ > 0 && iteration_ < nblocks_ - 1) {
        // The address in the cycle is one greater than the one's based iteration
        int index = blocks_[iteration_ + 1];
        if (iteration_ % 2 == 1) {
            // Odd iteration, read into B
            read(B_, block_starts_[index], block_sizes_[index]); 
        } else {
            // Even iteration, read into A
            read(A_, block_starts_[index], block_sizes_[index]); 
        }
    }

    // Determine current rows
    // Return index is in the cycle is one less than the one's based iteration 
    int index = blocks_[iteration_ - 1];
    current_rows_ = block_sizes_[index];
    
    // Return a block
    // Odd iteration, return A
    if (iteration_ % 2 == 1) { 
        return A_;
    // Even iteration, return B
    } else {
        return B_;
    }
}
void DFHFDiskIterator::read(shared_ptr<Matrix> A, int start, int rows)
{
    // Get the block pointer
    double** Ap = A->pointer();
    
    // Get an address
    psio_address addr = psio_get_address(PSIO_ZERO, sizeof(double)*start*ntri_);

    // Post the read
    aio_->read(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Ap[0], sizeof(double)*rows*ntri_, addr, &addr);
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
}

}}
