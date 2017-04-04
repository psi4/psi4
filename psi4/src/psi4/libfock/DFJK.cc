/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "jk_independent.h"
#include "link.h"
#include "direct_screening.h"
#include "cubature.h"
#include "points.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/cholesky.h"

#include <sstream>
#include "psi4/libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

DFJK::DFJK(std::shared_ptr<BasisSet> primary,
   std::shared_ptr<BasisSet> auxiliary) :
   JK(primary), auxiliary_(auxiliary)
{
    common_init();
}
DFJK::~DFJK()
{
}
void DFJK::common_init()
{
    df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = Process::environment.get_n_threads();
    #endif
    df_ints_io_ = "NONE";
    condition_ = 1.0E-12;
    unit_ = PSIF_DFSCF_BJ;
    is_core_ = true;
    psio_ = PSIO::shared_object();
}
SharedVector DFJK::iaia(SharedMatrix Ci, SharedMatrix Ca)
{
    // Target quantity
    Dimension dim(Ci->nirrep());
    for (int symm = 0; symm < Ci->nirrep(); symm++) {
        int rank = 0;
        for (int h = 0; h < Ci->nirrep(); h++) {
            rank += Ci->colspi()[h] * Ca->colspi()[h^symm];
        }
        dim[symm] = rank;
    }

    SharedVector Iia(new Vector("(ia|ia)", dim));

    // AO-basis quantities
    int nirrep = Ci->nirrep();
    int nocc = Ci->ncol();
    int nvir = Ca->ncol();
    int nso  = AO2USO_->rowspi()[0];

    SharedMatrix Ci_ao(new Matrix("Ci AO", nso, nocc));
    SharedMatrix Ca_ao(new Matrix("Ca AO", nso, nvir));
    SharedVector Iia_ao(new Vector("(ia|ia) AO", nocc*(ULI)nvir));

    int offset = 0;
    for (int h = 0; h < nirrep; h++) {
        int ni = Ci->colspi()[h];
        int nm = Ci->rowspi()[h];
        if (!ni || !nm) continue;
        double** Cip = Ci->pointer(h);
        double** Cp = Ci_ao->pointer();
        double** Up = AO2USO_->pointer(h);
        C_DGEMM('N','N',nso,ni,nm,1.0,Up[0],nm,Cip[0],ni,0.0,&Cp[0][offset],nocc);
        offset += ni;
    }

    offset = 0;
    for (int h = 0; h < nirrep; h++) {
        int ni = Ca->colspi()[h];
        int nm = Ca->rowspi()[h];
        if (!ni || !nm) continue;
        double** Cip = Ca->pointer(h);
        double** Cp = Ca_ao->pointer();
        double** Up = AO2USO_->pointer(h);
        C_DGEMM('N','N',nso,ni,nm,1.0,Up[0],nm,Cip[0],ni,0.0,&Cp[0][offset],nvir);
        offset += ni;
    }

    // Memory size
    int naux = auxiliary_->nbf();
    int maxrows = max_rows();

    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    // Temps
    #ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
    C_temp_.resize(omp_nthread_);
    Q_temp_.resize(omp_nthread_);
    #pragma omp parallel
    {
        C_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Ctemp", nocc, nso));
        Q_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Qtemp", maxrows, nso));
    }
    omp_set_num_threads(temp_nthread);
    #else
    for (int thread = 0; thread < omp_nthread_; thread++) {
        C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", nocc, nso)));
        Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", maxrows, nso)));
    }
    #endif

    E_left_ = SharedMatrix(new Matrix("E_left", nso, maxrows * nocc));
    E_right_ = SharedMatrix(new Matrix("E_right", nvir, maxrows * nocc));

    // Disk overhead
    psio_address addr = PSIO_ZERO;
    if (!is_core()) {
        Qmn_ = SharedMatrix(new Matrix("(Q|mn) Block", maxrows, num_nm));
        psio_->open(unit_,PSIO_OPEN_OLD);
    }

    // Blocks of Q
    double** Qmnp;
    double** Clp  = Ci_ao->pointer();
    double** Crp  = Ca_ao->pointer();
    double** Elp  = E_left_->pointer();
    double** Erp  = E_right_->pointer();
    double*  Iiap = Iia_ao->pointer();
    for (int Q = 0; Q < naux; Q += maxrows) {

        // Read block of (Q|mn) in
        int rows = (naux - Q <= maxrows ? naux - Q : maxrows);
        if (is_core()) {
            Qmnp = &Qmn_->pointer()[Q];
        } else {
            Qmnp = Qmn_->pointer();
            psio_->read(unit_,"(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]),sizeof(double)*naux*num_nm,addr,&addr);
        }

        // (mi|Q)
        #pragma omp parallel for schedule (dynamic)
        for (int m = 0; m < nso; m++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            double** Ctp = C_temp_[thread]->pointer();
            double** QSp = Q_temp_[thread]->pointer();

            const std::vector<int>& pairs = sieve_->function_to_function()[m];
            int mrows = pairs.size();

            for (int i = 0; i < mrows; i++) {
                int n = pairs[i];
                long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                C_DCOPY(rows,&Qmnp[0][ij],num_nm,&QSp[0][i],nso);
                C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nso);
            }

            C_DGEMM('N','T',nocc,rows,mrows,1.0,Ctp[0],nso,QSp[0],nso,0.0,&Elp[0][m*(ULI)nocc*rows],rows);
        }

        // (ai|Q)
        C_DGEMM('T','N',nvir,nocc*(ULI)rows,nso,1.0,Crp[0],nvir,Elp[0],nocc*(ULI)rows,0.0,Erp[0],nocc*(ULI)rows);

        // (ia|Q)(Q|ia)
        for (int i = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++) {
                double* Ep = &Erp[0][a * (ULI) nocc * rows + i * rows];
                Iiap[i * nvir + a] += C_DDOT(rows, Ep, 1, Ep, 1);
            }
        }

    }

    // Free disk overhead
    if (!is_core()) {
        Qmn_.reset();
        psio_->close(unit_,1);
    }

    // Free Temps
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();

    // SO-basis (ia|ia)
    Dimension i_offsets(Ci->nirrep());
    Dimension a_offsets(Ci->nirrep());
    for (int h = 1; h < Ci->nirrep(); h++) {
        i_offsets[h] = i_offsets[h-1] + Ci->colspi()[h-1];
        a_offsets[h] = a_offsets[h-1] + Ca->colspi()[h-1];
    }

    for (int symm = 0; symm < Ci->nirrep(); symm++) {
        offset = 0;
        for (int h = 0; h < Ci->nirrep(); h++) {
            int ni = Ci->colspi()[h];
            int na = Ca->colspi()[h^symm];
            int ioff = i_offsets[h];
            int aoff = a_offsets[h^symm];
            for (int i = 0; i < ni; i++) {
                for (int a = 0; a < na; a++) {
                    Iia->set(symm, i * na + a + offset, Iiap[(ioff + i) * nvir + (aoff + a)]);
                }
            }
            offset += ni * na;
        }
    }

    return Iia;
}
void DFJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf( "    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Algorithm:         %11s\n",  (is_core_ ? "Core" : "Disk"));
        outfile->Printf( "    Integral Cache:    %11s\n",  df_ints_io_.c_str());
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf( "    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf( "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
bool DFJK::is_core() const
{
    size_t ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    size_t mem = memory_;
    mem -= memory_overhead();
    mem -= memory_temp();

    // Two is for buffer space in fitting
    if (do_wK_)
        return (3L*three_memory + 2L*two_memory < memory_);
    else
        return (three_memory + 2L*two_memory < memory_);
}
unsigned long int DFJK::memory_temp() const
{
    unsigned long int mem = 0L;

    // J Overhead (Jtri, Dtri, d)
    mem += 2L * sieve_->function_pairs().size() + auxiliary_->nbf();
    // K Overhead (C_temp, Q_temp)
    mem += omp_nthread_ * (unsigned long int) primary_->nbf() * (auxiliary_->nbf() + max_nocc());

    return mem;
}
int DFJK::max_rows() const
{
    // Start with all memory
    unsigned long int mem = memory_;
    // Subtract J/K/wK/C/D overhead
    mem -= memory_overhead();
    // Subtract threading temp overhead
    mem -= memory_temp();

    // How much will each row cost?
    unsigned long int row_cost = 0L;
    // Copies of E tensor
    row_cost += (lr_symmetric_ ? 1L : 2L) * max_nocc() * primary_->nbf();
    // Slices of Qmn tensor, including AIO buffer (NOTE: AIO not implemented yet)
    row_cost += (is_core_ ? 1L : 1L) * sieve_->function_pairs().size();

    unsigned long int max_rows = mem / row_cost;

    if (max_rows > (unsigned long int) auxiliary_->nbf())
        max_rows = (unsigned long int) auxiliary_->nbf();
    if (max_rows < 1L)
        max_rows = 1L;

    return (int) max_rows;
}
int DFJK::max_nocc() const
{
    int max_nocc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc);
    }
    return max_nocc;
}
void DFJK::initialize_temps()
{
    J_temp_ = std::shared_ptr<Vector>(new Vector("Jtemp", sieve_->function_pairs().size()));
    D_temp_ = std::shared_ptr<Vector>(new Vector("Dtemp", sieve_->function_pairs().size()));
    d_temp_ = std::shared_ptr<Vector>(new Vector("dtemp", max_rows_));


    #ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
    C_temp_.resize(omp_nthread_);
    Q_temp_.resize(omp_nthread_);
    #pragma omp parallel
    {
        C_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf()));
        Q_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Qtemp", max_rows_, primary_->nbf()));
    }
    omp_set_num_threads(temp_nthread);
    #else
        for (int thread = 0; thread < omp_nthread_; thread++) {
            C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf())));
            Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", max_rows_, primary_->nbf())));
        }
    #endif

    E_left_ = SharedMatrix(new Matrix("E_left", primary_->nbf(), max_rows_ * max_nocc_));
    if (lr_symmetric_)
        E_right_ = E_left_;
    else
        E_right_ = std::shared_ptr<Matrix>(new Matrix("E_right", primary_->nbf(), max_rows_ * max_nocc_));

}
void DFJK::initialize_w_temps()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);

    #ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
        C_temp_.resize(omp_nthread_);
        Q_temp_.resize(omp_nthread_);
        #pragma omp parallel
        {
            C_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf()));
            Q_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Qtemp", max_rows_w, primary_->nbf()));
        }
    omp_set_num_threads(temp_nthread);
    #else
        for (int thread = 0; thread < omp_nthread_; thread++) {
            C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf())));
            Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", max_rows_w, primary_->nbf())));
        }
    #endif

    E_left_  = SharedMatrix(new Matrix("E_left", primary_->nbf(), max_rows_w * max_nocc_));
    E_right_ = SharedMatrix(new Matrix("E_right", primary_->nbf(), max_rows_w * max_nocc_));

}
void DFJK::free_temps()
{
    J_temp_.reset();
    D_temp_.reset();
    d_temp_.reset();
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();
}
void DFJK::free_w_temps()
{
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();
}
void DFJK::preiterations()
{

    // DF requires constant sieve, must be static throughout object life
    if (!sieve_) {
        sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    }

    // Core or disk?
    is_core_ =  is_core();


    if (is_core_)
        initialize_JK_core();
    else
        initialize_JK_disk();

    if (do_wK_) {
        if (is_core_)
            initialize_wK_core();
        else
            initialize_wK_disk();
    }
}

void DFJK::compute_JK()
{
    max_nocc_ = max_nocc();
    max_rows_ = max_rows();

    if (do_J_ || do_K_) {
        initialize_temps();
        if (is_core_)
            manage_JK_core();
        else
            manage_JK_disk();
        free_temps();
    }

    if (do_wK_) {
        initialize_w_temps();
        if (is_core_)
            manage_wK_core();
        else
            manage_wK_disk();
        free_w_temps();
        // Bring the wK matrices back to Hermitian
        if (lr_symmetric_) {
            for (size_t N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }
}
void DFJK::postiterations()
{
    Qmn_.reset();
    Qlmn_.reset();
    Qrmn_.reset();
}
void DFJK::initialize_JK_core()
{
    size_t ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    Qmn_ = SharedMatrix(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
        return;
    }

    const int primary_maxam = primary_->max_am();
    const int aux_maxam = auxiliary_->max_am();
    const int primary_nshell = primary_->nshell();
    const int aux_nshell = auxiliary_->nshell();

    const std::vector<long int>& schwarz_shell_pairs = sieve_->shell_pairs_reverse();
    const std::vector<long int>& schwarz_fun_pairs = sieve_->function_pairs_reverse();


    //Get a TEI for each thread
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri = new std::shared_ptr<TwoBodyAOInt>[nthread];
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    buffer[0] = eri[0]->buffer();

    for (int Q = 1; Q<nthread; Q++) {
        if(eri[0]->cloneable())
            eri[Q] = std::shared_ptr<TwoBodyAOInt>(eri[0]->clone());
        else
            eri[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());

        buffer[Q] = eri[Q]->buffer();
    }

    // shell pair blocks 
    std::vector<ShellPairBlock> p_blocks = eri[0]->get_blocks12();
    std::vector<ShellPairBlock> mn_blocks = eri[0]->get_blocks34();


    timer_on("JK: (A|mn)");


    // We are calculting integrals of (p|mn). Note that mn is on the ket side.
    // Integrals may be vectorized over the ket shells, so we want more of the
    // work to be on that side.
    //
    // Variable description:
    // mn = mu nu (in old terminology)
    // p_block_idx = index of what block of p we are doing
    // mn_block_idx = index of what block of mn we are doing
    // p_block = The block of p we are doing (a vector of std::pair<int, int>)
    // mn_block = The block of mn we are doing
    // mn_pair = a single pair of indicies, corresponding to the indicies of the
    //           m and n shells in the primary basis set
    // p_pair = a single pair of indicies, corresponding to the indicies of the
    //          p shell in the auxiliary basis set
    // num_m, num_n, num_p = number of basis functions in shell m, n, and p
    // m_start, n_start, p_start = starting index of the basis functions of this shell
    //                             relative to the beginning of their corresponding basis
    //                             sets.
    // im, in, ip = indices for looping over all the basis functions of the shells
    // im_idx, in_idx, ip_idx = indices of a particular basis function in a shell
    //                          relative to the start of the basis set
    //                          (ie, m_start + im)
    // sfp = screen function pair. The actual output index, taking into account
    //                             whether or not functions have been screened.

    #pragma omp parallel for schedule (dynamic) num_threads(nthread)
    for(size_t mn_block_idx = 0; mn_block_idx < mn_blocks.size(); mn_block_idx++)
    {
        #ifdef _OPENMP
            const int rank = omp_get_thread_num();
        #else
            const int rank = 0;
        #endif
    
        const auto & mn_block = mn_blocks[mn_block_idx];

        // loop over all the blocks of P
        for(int p_block_idx = 0; p_block_idx < p_blocks.size(); ++p_block_idx)
        {

            // compute the
            eri[rank]->compute_shell_blocks(p_block_idx, mn_block_idx);
            double const * my_buffer = buffer[rank];
 
            const auto & p_block = p_blocks[p_block_idx];

            for(const auto & mn_pair : mn_block)
            {
                const int m = mn_pair.first;
                const int n = mn_pair.second;
                const int num_m = primary_->shell(m).nfunction();
                const int num_n = primary_->shell(n).nfunction();
                const int m_start = primary_->shell(m).function_index();
                const int n_start = primary_->shell(n).function_index();
                const int num_mn = num_m * num_n;

                for(const auto & p_pair : p_block)
                {
                    // remember that this vector will only contain one shell pair
                    const int p = p_pair.first;
                    const int num_p = auxiliary_->shell(p).nfunction();
                    const int p_start = auxiliary_->shell(p).function_index();

                    for(int im = 0; im < num_m; ++im)
                    {
                        const int im_idx = m_start + im;
                        
                        for(int in = 0; in < num_n; ++in)
                        {
                            const int in_idx = n_start + in;
                            const int imn_idx = im*num_n + in;
                            const int sfp_idx = (im_idx*(im_idx+1))/2+in_idx; 

                            int sfp;

                            // note the assignment in the following conditional
                            if(im_idx >= in_idx && (sfp = schwarz_fun_pairs[sfp_idx]) > -1)
                            {
                                for(int ip = 0; ip < num_p; ++ip)
                                {
                                    const int ip_idx = p_start + ip;
                                    Qmnp[ip_idx][sfp] = my_buffer[ip*num_mn + imn_idx];
                                }
                            }
                        }
                    }

                    my_buffer += num_mn * num_p;
                }
            }
        }
    }

    timer_off("JK: (A|mn)");

    delete []buffer;
    delete []eri;

    timer_on("JK: (A|Q)^-1/2");

    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1/2");

    ULI max_cols = (memory_-three_memory-two_memory) / auxiliary_->nbf();
    if (max_cols < 1)
        max_cols = 1;
    if (max_cols > ntri)
        max_cols = ntri;
    SharedMatrix temp(new Matrix("Qmn buffer", auxiliary_->nbf(), max_cols));
    double** tempp = temp->pointer();

    size_t nblocks = ntri / max_cols;
    if ((ULI)nblocks*max_cols != ntri) nblocks++;

    size_t ncol = 0;
    size_t col = 0;

    timer_on("JK: (Q|mn)");

    for (size_t block = 0; block < nblocks; block++) {

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

    timer_off("JK: (Q|mn)");
    //Qmn_->print();

    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
    }
}
void DFJK::initialize_JK_disk()
{
    // Try to load
    if (df_ints_io_ == "LOAD") {
        return;
    }

    int nshell = primary_->nshell();
    int naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int,int> >& schwarz_shell_pairs = sieve_->shell_pairs();
    const std::vector<std::pair<int,int> >& schwarz_fun_pairs = sieve_->function_pairs();
    int nshellpairs = schwarz_shell_pairs.size();
    int ntri = schwarz_fun_pairs.size();
    const std::vector<long int>&  schwarz_shell_pairs_r = sieve_->shell_pairs_reverse();
    const std::vector<long int>&  schwarz_fun_pairs_r = sieve_->function_pairs_reverse();

    // ==> Memory Sizing <== //
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();
    ULI buffer_memory = memory_ - 2*two_memory; // Two is for buffer space in fitting

    //outfile->Printf( "Buffer memory = %ld words\n", buffer_memory);

    //outfile->Printf("Schwarz Shell Pairs:\n");
    //for (int MN = 0; MN < nshellpairs; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);
    //}

    //outfile->Printf("Schwarz Function Pairs:\n");
    //for (int MN = 0; MN < ntri; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);
    //}

    //outfile->Printf("Schwarz Reverse Shell Pairs:\n");
    //for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    //outfile->Printf("Schwarz Reverse Function Pairs:\n");
    //for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell
    std::shared_ptr<IntVector> MN_mem(new IntVector("Memory per MN pair", nshell * (nshell + 1) / 2));
    int *MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < ntri; mn++) {
        int m = schwarz_fun_pairs[mn].first;
        int n = schwarz_fun_pairs[mn].second;

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

    //outfile->Printf("  # Memory per M row #\n\n");
    //for (int M = 0; M < nshell; M++)
    //    outfile->Printf("   %3d: %10ld\n", M+1,M_memp[M]);
    //outfile->Printf("\n");

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
    std::shared_ptr<IntVector> MN_start(new IntVector("MUNU start per M row", nshell));
    int* MN_startp = MN_start->pointer();

    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (schwarz_shell_pairs[MN].first == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    std::shared_ptr<IntVector> mn_start(new IntVector("munu start per M row", nshell));
    int* mn_startp = mn_start->pointer();

    mn_startp[0] = schwarz_fun_pairs[0].first;
    int m_index = 1;
    for (int mn = 0; mn < ntri; mn++) {
        if (primary_->function_to_shell(schwarz_fun_pairs[mn].first) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row
    std::shared_ptr<IntVector> MN_col(new IntVector("MUNU cols per M row", nshell));
    int* MN_colp = MN_col->pointer();

    for (int M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1];
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];

    // Figure out the mn columns per M row
    std::shared_ptr<IntVector> mn_col(new IntVector("munu cols per M row", nshell));
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

    //outfile->Printf("Block, MN start, MN cols, mn start, mn cols\n");
    //for (int block = 0; block < nblock; block++) {
    //    outfile->Printf("  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block], mn_col_b[block]);
    //}
    //

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

    timer_on("JK: (A|Q)^-1");

    psio_->open(unit_,PSIO_OPEN_NEW);
    std::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Dispatch the prestripe
    aio->zero_disk(unit_,"(Q|mn) Integrals",naux,ntri);

    // Form the J symmetric inverse
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    timer_off("JK: (A|Q)^-1");

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    // ==> ERI initialization <== //
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block];
        int mn_start_val = mn_start_b[block];
        int MN_col_val = MN_col_b[block];
        int mn_col_val = mn_col_b[block];

        // ==> (A|mn) integrals <== //

        timer_on("JK: (A|mn)");

        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
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

        timer_off("JK: (A|mn)");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)");

        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val)
                cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        timer_off("JK: (Q|mn)");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn) Write");

        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"(Q|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }

        timer_off("JK: (Q|mn) Write");
    }

    // ==> Close out <== //
    Qmn_.reset();
    delete[] eri;
    delete[] buffer;

    psio_->close(unit_,1);
}
void DFJK::initialize_wK_core()
{
    int naux = auxiliary_->nbf();
    int ntri = sieve_->function_pairs().size();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif
    int rank = 0;

    // Check that the right integrals are using the correct omega
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_,1);
    }

    Qlmn_ = SharedMatrix(new Matrix("Qlmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qlmn_->pointer();

    Qrmn_ = SharedMatrix(new Matrix("Qrmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmn2p = Qrmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "Left (Q|w|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->read_entry(unit_, "Right (Q|w|mn) Integrals", (char*) Qmn2p[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
        return;
    }

    // => Left Integrals <= //

    //Get a TEI for each thread
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    const std::vector<long int>& schwarz_shell_pairs = sieve_->shell_pairs_reverse();
    const std::vector<long int>& schwarz_fun_pairs = sieve_->function_pairs_reverse();

    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    //The integrals (A|mn)

    timer_on("JK: (A|mn)^L");

    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^L");

    delete []buffer;
    delete []eri;

    // => Fitting <= //

    timer_on("JK: (A|Q)^-1");

    // Fitting metric
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_full_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1");

    timer_on("JK: (Q|mn)^L");

    // Fitting in one GEMM (being a clever bastard)
    C_DGEMM('N','N',naux,ntri,naux,1.0,Jinvp[0],naux,Qmn2p[0],ntri,0.0,Qmnp[0],ntri);

    timer_off("JK: (Q|mn)^L");

    // => Right Integrals <= //

    const double **buffer2 = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri2 = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    //The integrals (A|w|mn)

    timer_on("JK: (A|mn)^R");

    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri2[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer2[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^R");

    delete []buffer2;
    delete []eri2;

    // Try to save
    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->write_entry(unit_, "Left (Q|w|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->write_entry(unit_, "Right (Q|w|mn) Integrals", (char*) Qmn2p[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
        psio_->close(unit_,1);
    }
}
void DFJK::initialize_wK_disk()
{
    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_,1);
    }

    size_t nshell = primary_->nshell();
    size_t naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int,int> >& schwarz_shell_pairs = sieve_->shell_pairs();
    const std::vector<std::pair<int,int> >& schwarz_fun_pairs = sieve_->function_pairs();
    int nshellpairs = schwarz_shell_pairs.size();
    int ntri = schwarz_fun_pairs.size();
    const std::vector<long int>&  schwarz_shell_pairs_r = sieve_->shell_pairs_reverse();
    const std::vector<long int>&  schwarz_fun_pairs_r = sieve_->function_pairs_reverse();

    // ==> Memory Sizing <== //
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();
    ULI buffer_memory = memory_ - 2L*two_memory; // Two is for buffer space in fitting

    //outfile->Printf( "Buffer memory = %ld words\n", buffer_memory);

    //outfile->Printf("Schwarz Shell Pairs:\n");
    //for (int MN = 0; MN < nshellpairs; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);
    //}

    //outfile->Printf("Schwarz Function Pairs:\n");
    //for (int MN = 0; MN < ntri; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);
    //}

    //outfile->Printf("Schwarz Reverse Shell Pairs:\n");
    //for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    //outfile->Printf("Schwarz Reverse Function Pairs:\n");
    //for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell
    std::shared_ptr<IntVector> MN_mem(new IntVector("Memory per MN pair", nshell * (nshell + 1) / 2));
    int *MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < ntri; mn++) {
        int m = schwarz_fun_pairs[mn].first;
        int n = schwarz_fun_pairs[mn].second;

        int M = primary_->function_to_shell(m);
        int N = primary_->function_to_shell(n);

        MN_memp[M * (M + 1) / 2 + N] += naux;
    }

    //MN_mem->print(outfile);

    // Figure out exactly how much memory per M row
    ULI* M_memp = new ULI[nshell];
    memset(static_cast<void*>(M_memp), '\0', nshell*sizeof(ULI));

    for (size_t M = 0; M < nshell; M++) {
        for (size_t N = 0; N <= M; N++) {
            M_memp[M] += MN_memp[M * (M + 1) / 2 + N];
        }
    }

    //outfile->Printf("  # Memory per M row #\n\n");
    //for (int M = 0; M < nshell; M++)
    //    outfile->Printf("   %3d: %10ld\n", M+1,M_memp[M]);
    //outfile->Printf("\n");

    // Find and check the minimum required memory for this problem
    ULI min_mem = naux*(ULI) ntri;
    for (size_t M = 0; M < nshell; M++) {
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
    std::shared_ptr<IntVector> MN_start(new IntVector("MUNU start per M row", nshell));
    int* MN_startp = MN_start->pointer();

    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (schwarz_shell_pairs[MN].first == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    std::shared_ptr<IntVector> mn_start(new IntVector("munu start per M row", nshell));
    int* mn_startp = mn_start->pointer();

    mn_startp[0] = schwarz_fun_pairs[0].first;
    int m_index = 1;
    for (int mn = 0; mn < ntri; mn++) {
        if (primary_->function_to_shell(schwarz_fun_pairs[mn].first) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row
    std::shared_ptr<IntVector> MN_col(new IntVector("MUNU cols per M row", nshell));
    int* MN_colp = MN_col->pointer();

    for (size_t M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1];
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];

    // Figure out the mn columns per M row
    std::shared_ptr<IntVector> mn_col(new IntVector("munu cols per M row", nshell));
    int* mn_colp = mn_col->pointer();

    for (size_t M = 1; M < nshell; M++) {
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
    for (size_t M = 0; M < nshell; M++) {
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

    //outfile->Printf("Block, MN start, MN cols, mn start, mn cols\n");
    //for (int block = 0; block < nblock; block++) {
    //    outfile->Printf("  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block], mn_col_b[block]);
    //}
    //

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
    psio_->open(unit_,PSIO_OPEN_OLD);
    std::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Dispatch the prestripe
    aio->zero_disk(unit_,"Left (Q|w|mn) Integrals",naux,ntri);

    // Form the J full inverse
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_full_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    // ==> ERI initialization <== //
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block];
        int mn_start_val = mn_start_b[block];
        int MN_col_val = MN_col_b[block];
        int mn_col_val = mn_col_b[block];

        // ==> (A|mn) integrals <== //

        timer_on("JK: (A|mn)^L");

        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
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

        timer_off("JK: (A|mn)^L");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)^L");

        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= (size_t)mn_col_val)
                cols = mn_col_val - mn;

            for (size_t Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        timer_off("JK: (Q|mn)^L");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn)^L Write");

        psio_address addr;
        for (size_t Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"Left (Q|w|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }

        timer_off("JK: (Q|mn)^L Write");

    }

    Qmn_.reset();
    delete[] eri;

    // => Right Integrals <= //

    const double **buffer2 = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri2 = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    ULI maxP = auxiliary_->max_function_per_shell();
    ULI max_rows = memory_ / ntri;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    size_t counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    SharedMatrix Amn2(new Matrix("(A|mn) Block", max_rows, ntri));
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (size_t block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = auxiliary_->shell(Qstart).function_index();
        int nrows  = (Qstop == auxiliary_->nshell() ?
                     auxiliary_->nbf() -
                     auxiliary_->shell(Qstart).function_index() :
                     auxiliary_->shell(Qstop).function_index() -
                     auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (size_t QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq =  auxiliary_->shell(Q).function_index();
            int sm =  primary_->shell(M).function_index();
            int sn =  primary_->shell(N).function_index();

            eri2[thread]->compute_shell(Q,0,M,N);

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    if (m >= n && schwarz_fun_pairs_r[m*(m+1)/2 + n] >= 0) {
                        long int delta = schwarz_fun_pairs_r[m*(m+1)/2 + n];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] =
                            buffer2[thread][oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_,"Right (Q|w|mn) Integrals",(char*)Amn2p[0],sizeof(double)*nrows*ntri,next_AIA,&next_AIA);

        timer_off("JK: (Q|mn)^R Write");

    }
    Amn2.reset();
    delete[] eri2;
    delete[] buffer;
    delete[] buffer2;

    psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
    psio_->close(unit_,1);
}
void DFJK::rebuild_wK_disk()
{
    // Already open
    outfile->Printf( "    Rebuilding (Q|w|mn) Integrals (new omega)\n\n");

    size_t naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int,int> >& schwarz_fun_pairs = sieve_->function_pairs();
    int ntri = schwarz_fun_pairs.size();
    const std::vector<long int>&  schwarz_fun_pairs_r = sieve_->function_pairs_reverse();

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer2 = new const double*[nthread];
    std::shared_ptr<TwoBodyAOInt> *eri2 = new std::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    ULI maxP = auxiliary_->max_function_per_shell();
    ULI max_rows = memory_ / ntri;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        size_t nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    SharedMatrix Amn2(new Matrix("(A|mn) Block", max_rows, ntri));
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (size_t block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = auxiliary_->shell(Qstart).function_index();
        int nrows  = (Qstop == auxiliary_->nshell() ?
                     auxiliary_->nbf() -
                     auxiliary_->shell(Qstart).function_index() :
                     auxiliary_->shell(Qstop).function_index() -
                     auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (size_t QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq =  auxiliary_->shell(Q).function_index();
            int sm =  primary_->shell(M).function_index();
            int sn =  primary_->shell(N).function_index();

            eri2[thread]->compute_shell(Q,0,M,N);

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    if (m >= n && schwarz_fun_pairs_r[m*(m+1)/2 + n] >= 0) {
                        long int delta = schwarz_fun_pairs_r[m*(m+1)/2 + n];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] =
                            buffer2[thread][oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_,"Right (Q|w|mn) Integrals",(char*)Amn2p[0],sizeof(double)*nrows*ntri,next_AIA,&next_AIA);

        timer_off("JK: (Q|mn)^R Write");

    }
    Amn2.reset();
    delete[] eri2;
    delete[] buffer2;

    psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
    // No need to close
}
void DFJK::manage_JK_core()
{
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[Q],naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[Q],naux);
            timer_off("JK: K");
        }
    }
}
void DFJK::manage_JK_disk()
{
    int ntri = sieve_->function_pairs().size();
    Qmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_, ntri));
    psio_->open(unit_,PSIO_OPEN_OLD);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn) Read");
        psio_->read(unit_,"(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn) Read");

        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[0],naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[0],naux);
            timer_off("JK: K");
        }
    }
    psio_->close(unit_,1);
    Qmn_.reset();
}
void DFJK::manage_wK_core()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[Q],&Qrmn_->pointer()[Q],naux);
        timer_off("JK: wK");
    }
}
void DFJK::manage_wK_disk()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    int ntri = sieve_->function_pairs().size();
    Qlmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_w, ntri));
    Qrmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_w, ntri));
    psio_->open(unit_,PSIO_OPEN_OLD);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn)^L Read");
        psio_->read(unit_,"Left (Q|w|mn) Integrals", (char*)(Qlmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn)^L Read");

        addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn)^R Read");
        psio_->read(unit_,"Right (Q|w|mn) Integrals", (char*)(Qrmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn)^R Read");

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[0],&Qrmn_->pointer()[0],naux);
        timer_off("JK: wK");
    }
    psio_->close(unit_,1);
    Qlmn_.reset();
    Qrmn_.reset();
}
void DFJK::block_J(double** Qmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    unsigned long int num_nm = function_pairs.size();

    for (size_t N = 0; N < J_ao_.size(); N++) {

        double** Dp   = D_ao_[N]->pointer();
        double** Jp   = J_ao_[N]->pointer();
        double*  J2p  = J_temp_->pointer();
        double*  D2p  = D_temp_->pointer();
        double*  dp   = d_temp_->pointer();
        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            D2p[mn] = (m == n ? Dp[m][n] : Dp[m][n] + Dp[n][m]);
        }

        timer_on("JK: J1");
        C_DGEMV('N',naux,num_nm,1.0,Qmnp[0],num_nm,D2p,1,0.0,dp,1);
        timer_off("JK: J1");

        timer_on("JK: J2");
        C_DGEMV('T',naux,num_nm,1.0,Qmnp[0],num_nm,dp,1,0.0,J2p,1);
        timer_off("JK: J2");
        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void DFJK::block_K(double** Qmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    for (size_t N = 0; N < K_ao_.size(); N++) {

        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp  = C_left_ao_[N]->pointer();
        double** Crp  = C_right_ao_[N]->pointer();
        double** Elp  = E_left_->pointer();
        double** Erp  = E_right_->pointer();
        double** Kp   = K_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {

            timer_on("JK: K1");

            #pragma omp parallel for schedule (dynamic)
            for (int m = 0; m < nbf; m++) {

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = sieve_->function_to_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                    C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nbf);
                }

                C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux);
            }

            timer_off("JK: K1");

        }

        if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N-1].get())) {

            if (C_right_[N].get() == C_left_[N].get()) {
                ::memcpy((void*) Erp[0], (void*) Elp[0], sizeof(double) * naux * nocc * nbf);
            } else {

                timer_on("JK: K1");

                #pragma omp parallel for schedule (dynamic)
                for (int m = 0; m < nbf; m++) {

                    int thread = 0;
                    #ifdef _OPENMP
                        thread = omp_get_thread_num();
                    #endif

                    double** Ctp = C_temp_[thread]->pointer();
                    double** QSp = Q_temp_[thread]->pointer();

                    const std::vector<int>& pairs = sieve_->function_to_function()[m];
                    int rows = pairs.size();

                    for (int i = 0; i < rows; i++) {
                        int n = pairs[i];
                        long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                        C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                        C_DCOPY(nocc,Crp[n],1,&Ctp[0][i],nbf);
                    }

                    C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux);
                }

                timer_off("JK: K1");

            }

        }

        timer_on("JK: K2");
        C_DGEMM('N','T',nbf,nbf,naux*nocc,1.0,Elp[0],naux*nocc,Erp[0],naux*nocc,1.0,Kp[0],nbf);
        timer_off("JK: K2");
    }

}
void DFJK::block_wK(double** Qlmnp, double** Qrmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    for (size_t N = 0; N < wK_ao_.size(); N++) {

        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp  = C_left_ao_[N]->pointer();
        double** Crp  = C_right_ao_[N]->pointer();
        double** Elp  = E_left_->pointer();
        double** Erp  = E_right_->pointer();
        double** wKp   = wK_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {

            timer_on("JK: wK1");

            #pragma omp parallel for schedule (dynamic)
            for (int m = 0; m < nbf; m++) {

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = sieve_->function_to_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux,&Qlmnp[0][ij],num_nm,&QSp[0][i],nbf);
                    C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nbf);
                }

                C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux);
            }

            timer_off("JK: wK1");

        }

        timer_on("JK: wK1");

        #pragma omp parallel for schedule (dynamic)
        for (int m = 0; m < nbf; m++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            double** Ctp = C_temp_[thread]->pointer();
            double** QSp = Q_temp_[thread]->pointer();

            const std::vector<int>& pairs = sieve_->function_to_function()[m];
            int rows = pairs.size();

            for (int i = 0; i < rows; i++) {
                int n = pairs[i];
                long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                C_DCOPY(naux,&Qrmnp[0][ij],num_nm,&QSp[0][i],nbf);
                C_DCOPY(nocc,Crp[n],1,&Ctp[0][i],nbf);
            }

            C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux);
        }

        timer_off("JK: wK1");

        timer_on("JK: wK2");
        C_DGEMM('N','T',nbf,nbf,naux*nocc,1.0,Elp[0],naux*nocc,Erp[0],naux*nocc,1.0,wKp[0],nbf);
        timer_off("JK: wK2");
    }
}
}
