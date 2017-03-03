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

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "thce.h"
#include "lreri.h"
#include "psi4/psi4-dec.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;


namespace psi {

LRERI::LRERI(std::shared_ptr<BasisSet> primary) :
    primary_(primary)
{
    common_init();
}
LRERI::~LRERI()
{
}
void LRERI::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    memory_ = 32000000L;
}
void LRERI::load_wavefunction(std::shared_ptr<Wavefunction> ref)
{
    if (ref->Ca() == ref->Cb()) {
        set_C(ref->Ca_subset("AO", "ALL"));

        int nfocc = ref->frzcpi().sum();
        int nfvir = ref->frzvpi().sum();
        int naocc = ref->nalphapi().sum() - nfocc;
        int navir = ref->nmopi().sum() - naocc - nfocc - nfvir;

        add_space("FROZEN_OCC", 0, nfocc);
        add_space("ACTIVE_OCC", nfocc, nfocc + naocc);
        add_space("ACTIVE_VIR", nfocc + naocc, nfocc + naocc + navir);
        add_space("ACTIVE_ALL", nfocc, nfocc + naocc + navir);
        add_space("FROZEN_VIR", nfocc + naocc + navir, nfocc + naocc + navir + nfvir);
        add_space("OCC", 0, nfocc + naocc);
        add_space("VIR", nfocc + naocc, nfocc + naocc + navir + nfvir);
        add_space("ALL", 0, nfocc + naocc + navir + nfvir);
    } else {
        std::vector<std::shared_ptr<Matrix> > list;
        list.push_back(ref->Ca_subset("AO", "ALL"));
        list.push_back(ref->Cb_subset("AO", "ALL"));
        set_C(Matrix::horzcat(list));

        int offset = ref->nmopi().sum();
        int nfocc_a = ref->frzcpi().sum();
        int nfvir_a = ref->frzvpi().sum();
        int naocc_a = ref->nalphapi().sum() - nfocc_a;
        int navir_a = ref->nmopi().sum() - naocc_a - nfocc_a - nfvir_a;
        int nfocc_b = ref->frzcpi().sum();
        int nfvir_b = ref->frzvpi().sum();
        int naocc_b = ref->nbetapi().sum() - nfocc_b;
        int navir_b = ref->nmopi().sum() - naocc_b - nfocc_b - nfvir_b;

        add_space("FROZEN_OCC_A", 0, nfocc_a);
        add_space("ACTIVE_OCC_A", nfocc_a, nfocc_a + naocc_a);
        add_space("ACTIVE_VIR_A", nfocc_a + naocc_a, nfocc_a + naocc_a + navir_a);
        add_space("ACTIVE_ALL_A", nfocc_a, nfocc_a + naocc_a + navir_a);
        add_space("FROZEN_VIR_A", nfocc_a + naocc_a + navir_a, nfocc_a + naocc_a + navir_a + nfvir_a);
        add_space("OCC_A", 0, nfocc_a + naocc_a);
        add_space("VIR_A", nfocc_a + naocc_a, nfocc_a + naocc_a + navir_a + nfvir_a);
        add_space("ALL_A", 0, nfocc_a + naocc_a + navir_a + nfvir_a);
        add_space("FROZEN_OCC_B", offset, offset + nfocc_b);
        add_space("ACTIVE_OCC_B", offset + nfocc_b, offset + nfocc_b + naocc_b);
        add_space("ACTIVE_VIR_B", offset + nfocc_b + naocc_b, offset + nfocc_b + naocc_b + navir_b);
        add_space("ACTIVE_ALL_B", offset + nfocc_b, offset + nfocc_b + naocc_b + navir_b);
        add_space("FROZEN_VIR_B", offset + nfocc_b + naocc_b + navir_b, offset + nfocc_b + naocc_b + navir_b + nfvir_b);
        add_space("OCC_B", offset, offset + nfocc_b + naocc_b);
        add_space("VIR_B", offset + nfocc_b + naocc_b, offset + nfocc_b + naocc_b + navir_b + nfvir_b);
        add_space("ALL_B", offset, offset + nfocc_b + naocc_b + navir_b + nfvir_b);
    }
}
void LRERI::load_options(Options& options)
{
    print_ = options.get_int("PRINT");
    debug_ = options.get_int("DEBUG");
    bench_ = options.get_int("BENCH");
    memory_ = (0.9 * Process::environment.get_memory() / 8L);
}
void LRERI::set_C(std::shared_ptr<Matrix> C)
{
    clear();
    C_ = C;
}
void LRERI::add_space(const std::string& key, int start, int end)
{
    spaces_[key] = std::pair<int,int>(start,end);
    spaces_order_.push_back(key);
}
void LRERI::clear()
{
    C_.reset();
    spaces_.clear();
    spaces_order_.clear();
}
std::shared_ptr<Matrix> LRERI::Jm12(std::shared_ptr<BasisSet> auxiliary, double J_cutoff)
{
    // Everybody likes them some inverse square root metric, eh?

    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif

    int naux = auxiliary->nbf();

    std::shared_ptr<Matrix> J(new Matrix("J", naux, naux));
    double** Jp = J->pointer();

    std::shared_ptr<IntegralFactory> Jfactory(new IntegralFactory(auxiliary, BasisSet::zero_ao_basis_set(), auxiliary, BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jeri;
    for (int thread = 0; thread < nthread; thread++) {
        Jeri.push_back(std::shared_ptr<TwoBodyAOInt>(Jfactory->eri()));
    }

    std::vector<std::pair<int, int> > Jpairs;
    for (int M = 0; M < auxiliary->nshell(); M++) {
        for (int N = 0; N <= M; N++) {
            Jpairs.push_back(std::pair<int,int>(M,N));
        }
    }
    long int num_Jpairs = Jpairs.size();

    #pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (long int PQ = 0L; PQ < num_Jpairs; PQ++) {

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        std::pair<int,int> pair = Jpairs[PQ];
        int P = pair.first;
        int Q = pair.second;

        Jeri[thread]->compute_shell(P,0,Q,0);

        int np = auxiliary->shell(P).nfunction();
        int op = auxiliary->shell(P).function_index();
        int nq = auxiliary->shell(Q).nfunction();
        int oq = auxiliary->shell(Q).function_index();

        const double* buffer = Jeri[thread]->buffer();

        for (int p = 0; p < np; p++) {
        for (int q = 0; q < nq; q++) {
            Jp[p + op][q + oq] =
            Jp[q + oq][p + op] =
                (*buffer++);
        }}
    }
    Jfactory.reset();
    Jeri.clear();

    // > Invert J < //

    J->power(-1.0/2.0, J_cutoff);

    return J;
}


DFERI::DFERI(std::shared_ptr<BasisSet> primary,
             std::shared_ptr<BasisSet> auxiliary) :
    LRERI(primary), auxiliary_(auxiliary)
{
    common_init();
}
DFERI::~DFERI()
{
}
void DFERI::common_init()
{
    schwarz_cutoff_ = 0.0;
    J_cutoff_ = 1.0E-10;
    keep_raw_integrals_ = false;
    omega_ = 0.0;
}
std::shared_ptr<DFERI> DFERI::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
{
    std::shared_ptr<DFERI> df(new DFERI(primary,auxiliary));
    df->load_options(options);
    return df;
}
std::shared_ptr<DFERI> DFERI::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options, std::shared_ptr<Wavefunction> ref)
{
    std::shared_ptr<DFERI> df = DFERI::build(primary,auxiliary,options);
    df->load_wavefunction(ref);
    return df;
}
void DFERI::load_options(Options& options)
{
    LRERI::load_options(options);
    J_cutoff_ = options.get_double("DF_FITTING_CONDITION");
    schwarz_cutoff_ = options.get_double("INTS_TOLERANCE");
}
void DFERI::add_pair_space(const std::string& name, const std::string& space1, const std::string& space2, double power, bool transpose12)
{
    pair_spaces_order_.push_back(name);
    pair_spaces_[name] = std::pair<std::string,std::string>(space1,space2);
    pair_powers_[name] = power;
    pair_transposes_[name] = transpose12;
}
void DFERI::clear_pair_spaces()
{
    pair_spaces_.clear();
    pair_spaces_order_.clear();
    pair_powers_.clear();
    pair_transposes_.clear();
    ints_.clear();
}
void DFERI::clear()
{
    clear_pair_spaces();
    LRERI::clear();
}
void DFERI::print_header(int level)
{
    outfile->Printf( "  ==> DFERI: Density Fitted 3-Index Tensors <==\n\n");
    if (omega_ != 0.0) {
        outfile->Printf( "    LRC Omega      = %11.3E\n", omega_);
    }
    outfile->Printf( "    Schwarz cutoff = %11.3E\n", schwarz_cutoff_);
    outfile->Printf( "    J cutoff       = %11.3E\n", J_cutoff_);
    outfile->Printf( "    Mem (GB)       = %11zu\n", (memory_ * 8L / 1073741824L));
    outfile->Printf( "\n");

    if (level > 1) {
        outfile->Printf( "   => Primary Basis <=\n\n");
        primary_->print_by_level("outfile", print_);
    }

    outfile->Printf( "   => Auxiliary Basis <=\n\n");
    auxiliary_->print_by_level("outfile", print_);

    if (level > 1) {
        outfile->Printf( "   => Orbital Spaces: <=\n\n");
        outfile->Printf( "    %12s %12s %12s\n", "Space", "Start", "End");
        for (int i = 0; i < spaces_order_.size(); i++) {
            outfile->Printf( "    %12s %12d %12d\n", spaces_order_[i].c_str(), spaces_[spaces_order_[i]].first, spaces_[spaces_order_[i]].second);
        }
        outfile->Printf( "\n");
    }

    if (level > 1) {
        outfile->Printf( "   => Required Orbital Pair Spaces: <=\n\n");
        outfile->Printf( "    %12s %12s %12s %11s %11s\n", "Tensor", "Space 1", "Space 2", "J Power", "Transpose12");
        for (int i = 0; i < pair_spaces_order_.size(); i++) {
            outfile->Printf( "    %12s %12s %12s %11.3E %11s\n", pair_spaces_order_[i].c_str(), pair_spaces_[pair_spaces_order_[i]].first.c_str(), pair_spaces_[pair_spaces_order_[i]].second.c_str(),
                pair_powers_[pair_spaces_order_[i]],
                pair_transposes_[pair_spaces_order_[i]]);
        }
        outfile->Printf( "\n");
    }
}
int DFERI::size_Q(void)
{
    return auxiliary_->nbf();
}
std::shared_ptr<Matrix> DFERI::Jpow(double power)
{
    // Everybody likes them some inverse square root metric, eh?

    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif

    int naux = auxiliary_->nbf();

    std::shared_ptr<Matrix> J(new Matrix("J", naux, naux));
    double** Jp = J->pointer();

    std::shared_ptr<IntegralFactory> Jfactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jeri;
    for (int thread = 0; thread < nthread; thread++) {
        if (omega_ == 0.0) {
            Jeri.push_back(std::shared_ptr<TwoBodyAOInt>(Jfactory->eri()));
        } else {
            Jeri.push_back(std::shared_ptr<TwoBodyAOInt>(Jfactory->erf_eri(omega_)));
        }
    }

    std::vector<std::pair<int, int> > Jpairs;
    for (int M = 0; M < auxiliary_->nshell(); M++) {
        for (int N = 0; N <= M; N++) {
            Jpairs.push_back(std::pair<int,int>(M,N));
        }
    }
    long int num_Jpairs = Jpairs.size();

    #pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (long int PQ = 0L; PQ < num_Jpairs; PQ++) {

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        std::pair<int,int> pair = Jpairs[PQ];
        int P = pair.first;
        int Q = pair.second;

        Jeri[thread]->compute_shell(P,0,Q,0);

        int np = auxiliary_->shell(P).nfunction();
        int op = auxiliary_->shell(P).function_index();
        int nq = auxiliary_->shell(Q).nfunction();
        int oq = auxiliary_->shell(Q).function_index();

        const double* buffer = Jeri[thread]->buffer();

        for (int p = 0; p < np; p++) {
        for (int q = 0; q < nq; q++) {
            Jp[p + op][q + oq] =
            Jp[q + oq][p + op] =
                (*buffer++);
        }}
    }
    Jfactory.reset();
    Jeri.clear();

    // > Invert J < //
    if (power != 1.0) {
        J->power(power, J_cutoff_);
    }

    return J;
}
void DFERI::compute()
{
    // => Allocation <= //

    allocate();

    // => (A|mn) C_mp C_nq => (A|pq) => (pq|A) <= //

    transform();

    // => (pq|A) J_{AQ}^-1/2 => (pq|Q) <= //

    fit();
}
void DFERI::allocate()
{
    ints_.clear();
    for (int i = 0; i < pair_spaces_order_.size(); i++) {
        std::string name = pair_spaces_order_[i];
        std::string space1 = pair_spaces_[name].first;
        std::string space2 = pair_spaces_[name].second;

        int size1 = spaces_[space1].second - spaces_[space1].first;
        int size2 = spaces_[space2].second - spaces_[space2].first;

        std::string naux = "NAUX";

        if (pair_transposes_[name]) {
            std::swap(size1,size2);
            std::swap(space1,space2);
        }

        ints_[name + "_temp"] = DiskTensor::build(name + "_temp",
            naux, auxiliary_->nbf(),
            space1, size1,
            space2, size2);
        ints_[name] = DiskTensor::build(name,
            space1, size1,
            space2, size2,
            naux, auxiliary_->nbf());
    }
}
void DFERI::transform()
{
    // > Sizing < //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // > Threading < //

    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif

    // > Task merging < //

    std::vector<std::vector<std::string> > tasks;
    std::vector<std::string> rem = pair_spaces_order_;

    // Lump all common first half-transforms together
    while (rem.size()) {
        std::vector<std::string> task;
        std::vector<std::string> rem2;
        std::string space1 = pair_spaces_[rem[0]].first;
        for (int i = 0; i < rem.size(); i++) {
            if (space1 == pair_spaces_[rem[i]].first)
                task.push_back(rem[i]);
            else
                rem2.push_back(rem[i]);
        }
        tasks.push_back(task);
        rem = rem2;
    }

    // > Maximum orbital sizes < //

    size_t max1 = 0L;
    size_t max12 = 0L;

    for (int i = 0; i < pair_spaces_order_.size(); i++) {
        std::string name = pair_spaces_order_[i];
        std::string space1 = pair_spaces_[name].first;
        std::string space2 = pair_spaces_[name].second;

        int size1 = spaces_[space1].second - spaces_[space1].first;
        int size2 = spaces_[space2].second - spaces_[space2].first;

        size_t size12 = size1 * (size_t) size2;

        max1 = (max1 < size1 ? size1 : max1);
        max12 = (max12 < size12 ? size12 : max12);
    }

    // > Row requirements < //

    unsigned long int per_row = 0L;
    // (Q|mn)
    per_row += nso * (unsigned long int) nso;
    // (Q|mi)
    per_row += max1 * (unsigned long int) nso;
    // (Q|ia)
    per_row += max12;

    // > Maximum number of rows < //

    unsigned long int max_rows = (memory_ / per_row);
    //max_rows = 3L * auxiliary_->max_function_per_shell(); // Debug
    if (max_rows < auxiliary_->max_function_per_shell()) {
        throw PSIEXCEPTION("Out of memory in DFERI.");
    }
    max_rows = (max_rows > auxiliary_->nbf() ? auxiliary_->nbf() : max_rows);

    // > Shell block assignments < //

    std::vector<int> shell_starts;
    shell_starts.push_back(0);
    int fcount = auxiliary_->shell(0).nfunction();
    for (int Q = 1; Q < auxiliary_->nshell(); Q++) {
        if (fcount + auxiliary_->shell(Q).nfunction() > max_rows) {
            shell_starts.push_back(Q);
            fcount = auxiliary_->shell(Q).nfunction();
        } else {
            fcount += auxiliary_->shell(Q).nfunction();
        }
    }
    shell_starts.push_back(auxiliary_->nshell());

    // > Task printing (Debug) < //

    //outfile->Printf("Auxiliary Composition:\n\n");
    //for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
    //    outfile->Printf("%3d: %2d\n", Q, auxiliary_->shell(Q).nfunction());
    //}
    //outfile->Printf("\n");

    //outfile->Printf("Max Rows: %zu\n\n", max_rows);

    //outfile->Printf("Task Starts:\n\n");
    //for (int task = 0; task < shell_starts.size() - 1; task++) {
    //    outfile->Printf("%3d: %3d\n", task, shell_starts[task]);
    //}
    //outfile->Printf("\n");

    // => ERI Objects <= //

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int thread = 0; thread < nthread; thread++) {
        if (omega_ == 0.0) {
            eri.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        } else {
            eri.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
        }
    }

    // => ERI Sieve <= //

    std::shared_ptr<ERISieve> sieve(new ERISieve(primary_, schwarz_cutoff_));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    long int nshell_pairs = (long int) shell_pairs.size();

    // => Temporary Tensors <= //

    // > Three-index buffers < //
    std::shared_ptr<Matrix> Amn(new Matrix("(A|mn)", max_rows, nso * (unsigned long int) nso));
    std::shared_ptr<Matrix> Ami(new Matrix("(A|mi)", max_rows, nso * (unsigned long int) max1));
    std::shared_ptr<Matrix> Aia(new Matrix("(A|ia)", max_rows, max12));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    // > C-matrix weirdness < //

    double** Cp = C_->pointer();
    int lda = C_->colspi()[0];

    // ==> Master Loop <== //

    for (int block = 0; block < shell_starts.size() - 1; block++) {

        // > Block characteristics < //

        int Pstart = shell_starts[block];
        int Pstop  = shell_starts[block+1];
        int nPshell = Pstop - Pstart;
        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? auxiliary_->nbf() : auxiliary_->shell(Pstop).function_index());
        int rows = pstop - pstart;

        // > (Q|mn) ERIs < //

        ::memset((void*) Amnp[0], '\0', sizeof(double) * rows * nso * nso);

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int PMN = 0L; PMN < nPshell * nshell_pairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P  = PMN / nshell_pairs + Pstart;
            int MN = PMN % nshell_pairs;
            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            eri[thread]->compute_shell(P,0,M,N);

            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();
            int np = auxiliary_->shell(P).nfunction();
            int om = primary_->shell(M).function_index();
            int on = primary_->shell(N).function_index();
            int op = auxiliary_->shell(P).function_index();

            const double* buffer = eri[thread]->buffer();

            for (int p = 0; p < np; p++) {
            for (int m = 0; m < nm; m++) {
            for (int n = 0; n < nn; n++) {
                Amnp[p + op - pstart][(m + om) * nso + (n + on)] =
                Amnp[p + op - pstart][(n + on) * nso + (m + om)] =
                (*buffer++);
            }}}
        }

        for (int ind1 = 0; ind1 < tasks.size(); ind1++) {

            std::string space1 = pair_spaces_[tasks[ind1][0]].first;
            int start1 = spaces_[space1].first;
            int end1   = spaces_[space1].second;
            int n1      = end1 - start1;
            double* C1p = &Cp[0][start1];

            C_DGEMM('N','N',rows*nso,n1,nso,1.0,Amnp[0],nso,C1p,lda,0.0,Amip[0],n1);

            for (int ind2 = 0; ind2 < tasks[ind1].size(); ind2++) {
                std::string space2 = pair_spaces_[tasks[ind1][ind2]].second;
                int start2 = spaces_[space2].first;
                int end2   = spaces_[space2].second;
                int n2      = end2 - start2;
                double* C2p = &Cp[0][start2];

                size_t n12 = n1 * (size_t) n2;
                size_t no1 = nso * (size_t) n1;

                std::string name = tasks[ind1][ind2];
                bool transpose12 = pair_transposes_[name];

                if (transpose12) {
                    #pragma omp parallel for num_threads(nthread)
                    for (int Q = 0; Q < rows; Q++) {
                        C_DGEMM('T','N',n2,n1,nso,1.0,C2p,lda,Amip[0] + Q*no1,n1,0.0,Aiap[0] + Q*n12,n1);
                    }
                } else {
                    #pragma omp parallel for num_threads(nthread)
                    for (int Q = 0; Q < rows; Q++) {
                        C_DGEMM('T','N',n1,n2,nso,1.0,Amip[0] + Q*no1,n1,C2p,lda,0.0,Aiap[0] + Q*n12,n2);
                    }
                }

                //Amn->print();
                //Ami->print();
                //Aia->print();

                std::shared_ptr<Tensor> A = ints_[name + "_temp"];
                FILE* fh = A->file_pointer();
                fwrite(Aiap[0],sizeof(double),rows*n12,fh);
            }
        }
    }
}
void DFERI::fit()
{
    int naux = auxiliary_->nbf();

    size_t max_pairs = 0L;
    for (int i = 0; i < pair_spaces_order_.size(); i++) {
        std::string name = pair_spaces_order_[i];
        std::shared_ptr<Tensor> A = ints_[name];
        size_t pairs = A->sizes()[0] * (size_t) A->sizes()[1];
        max_pairs = (pairs > max_pairs ? pairs : max_pairs);
    }

    size_t mem = memory_;
    mem -= naux * (size_t) naux;
    size_t max_rows = mem / (2L * naux);
    max_rows = (max_rows > max_pairs ? max_pairs : max_rows);

    std::shared_ptr<Matrix> T1(new Matrix("T1", naux, max_rows));
    std::shared_ptr<Matrix> T2(new Matrix("T2", max_rows, naux));
    double** T1p = T1->pointer();
    double** T2p = T2->pointer();

    std::set<double> unique_pows;
    for (int i = 0; i < pair_spaces_order_.size(); i++) {
        std::string name = pair_spaces_order_[i];
        double power = pair_powers_[name];
        unique_pows.insert(power);
    }

    for (std::set<double>::iterator it = unique_pows.begin();
        it != unique_pows.end(); ++it) {

        double power = (*it);

        std::shared_ptr<Matrix> J = Jpow(power);
        double** Jp = J->pointer();

        for (int i = 0; i < pair_spaces_order_.size(); i++) {
            std::string name = pair_spaces_order_[i];
            double powi = pair_powers_[name];

            if (powi != power) continue;

            std::shared_ptr<Tensor> A = ints_[name];
            size_t pairs = A->sizes()[0] * (size_t) A->sizes()[1];

            std::shared_ptr<Tensor> AT = ints_[name + "_temp"];

            FILE* fh = A->file_pointer();
            FILE* fhT = AT->file_pointer();

            for (size_t pair = 0L; pair < pairs; pair += max_rows) {
                size_t npairs = (pair + max_rows >= pairs? pairs - pair : max_rows);

                double* Ttp = T1p[0];
                for (int Q = 0; Q < naux; Q++) {
                    fseek(fhT,Q*pairs*sizeof(double) + pair*sizeof(double),SEEK_SET);
                    fread(Ttp,sizeof(double),npairs,fhT);
                    Ttp += npairs;
                }

                //T1->print();

                C_DGEMM('T','N',npairs,naux,naux,1.0,T1p[0],npairs,Jp[0],naux,0.0,T2p[0],naux);

                //T2->print();

                fwrite(T2p[0],sizeof(double),npairs*(size_t)naux,fh);
            }

            if (!keep_raw_integrals_) {
                ints_.erase(name + "_temp");
            }
        }
    }
}

LSTHCERI::LSTHCERI(std::shared_ptr<BasisSet> primary,
    std::shared_ptr<BasisSet> auxiliary,
    std::shared_ptr<Matrix> X) :
    LRERI(primary),
    auxiliary_(auxiliary),
    X_(X)
{
    common_init();
}
LSTHCERI::~LSTHCERI()
{
}
void LSTHCERI::common_init()
{
    S_cutoff_ = 1.0E-10;
    J_cutoff_ = 1.0E-10;
    schwarz_cutoff_ = 0.0;
    balance_ = false;
}
std::shared_ptr<LSTHCERI> LSTHCERI::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X, Options& options)
{
    std::shared_ptr<LSTHCERI> df(new LSTHCERI(primary,auxiliary,X));
    df->load_options(options);
    return df;
}
std::shared_ptr<LSTHCERI> LSTHCERI::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X, Options& options, std::shared_ptr<Wavefunction> ref)
{
    std::shared_ptr<LSTHCERI> df = LSTHCERI::build(primary,auxiliary,X,options);
    df->load_wavefunction(ref);
    return df;
}
void LSTHCERI::load_options(Options& options)
{
    LRERI::load_options(options);
    J_cutoff_ = options.get_double("THC_J_CUTOFF");
    S_cutoff_ = options.get_double("THC_S_CUTOFF");
    schwarz_cutoff_ = options.get_double("INTS_TOLERANCE");
    balance_ = options.get_bool("THC_BALANCE");
}
void LSTHCERI::add_eri_space(const std::string& name, const std::string& space1, const std::string& space2, const std::string& space3, const std::string& space4)
{
    eri_spaces_order_.push_back(name);
    std::vector<std::string> task;
    task.push_back(space1);
    task.push_back(space2);
    task.push_back(space3);
    task.push_back(space4);
    eri_spaces_[name] = task;
}
void LSTHCERI::clear()
{
    LRERI::clear();
    eri_spaces_.clear();
    eri_spaces_order_.clear();
    ints_.clear();
}
void LSTHCERI::print_header(int level)
{
    outfile->Printf( "  ==> LSTHCERI: LS-THC 2-Index Tensors <==\n\n");

    outfile->Printf( "    Schwarz cutoff = %11.3E\n", schwarz_cutoff_);
    outfile->Printf( "    J cutoff       = %11.3E\n", J_cutoff_);
    outfile->Printf( "    S cutoff       = %11.3E\n", S_cutoff_);
    outfile->Printf( "    Balance        = %11s\n", (balance_ ? "Yes" : "No"));
    outfile->Printf( "    Mem (GB)       = %11zu\n", (memory_ * 8L / 1073741824L));
    outfile->Printf( "\n");

    if (level > 1) {
        outfile->Printf( "   => Primary Basis <=\n\n");
        primary_->print_by_level("outfile", print_);
    }

    if (auxiliary_) {
        outfile->Printf( "   => Auxiliary Basis <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }

    if (level > 1) {
        outfile->Printf( "   => Orbital Spaces: <=\n\n");
        outfile->Printf( "    %12s %12s %12s\n", "Space", "Start", "End");
        for (int i = 0; i < spaces_order_.size(); i++) {
            outfile->Printf( "    %12s %12d %12d\n", spaces_order_[i].c_str(), spaces_[spaces_order_[i]].first, spaces_[spaces_order_[i]].second);
        }
        outfile->Printf( "\n");
    }

    if (level > 1) {
        outfile->Printf( "   => Required ERI Spaces: <=\n\n");
        outfile->Printf( "    %12s %12s %12s %12s %12s\n", "Tensor", "Space 1", "Space 2", "Space 3", "Space 4");
        for (int i = 0; i < eri_spaces_order_.size(); i++) {
            std::string tensor = eri_spaces_order_[i];
            std::vector<std::string> task = eri_spaces_[tensor];
            outfile->Printf( "    %12s %12s %12s %12s %12s\n", tensor.c_str(), task[0].c_str(), task[1].c_str(), task[2].c_str(), task[3].c_str());
        }
        outfile->Printf( "\n");
    }
}
void LSTHCERI::compute()
{
    //print_header();
    ints_.clear();

    // => Roll some X matrices <= //

    std::map<std::string, std::shared_ptr<Tensor> > Xs = build_X();

    // => Grow some L matrices <= //

    std::map<std::string, std::shared_ptr<Tensor> > Es = build_E(Xs);
    std::map<std::string, std::shared_ptr<Tensor> > Ss = build_S(Xs);
    std::map<std::string, std::shared_ptr<Tensor> > Ls = build_L(Es,Ss);
    Es.clear();

    // => Slam some Z matrices together <= //

    std::map<std::string, std::shared_ptr<Tensor> > Zs = build_Z(Ls);

    // => Pack this guy up <= //

    pack(Xs,Zs,Ls,Ss);
}
void LSTHCERI::compute_meth()
{
    print_header();
    meths_.clear();

    // => Roll some X matrices <= //

    std::map<std::string, std::shared_ptr<Tensor> > Xs = build_X(true);
    std::map<std::string, std::shared_ptr<Tensor> > Ss = build_S(Xs, true);
    pack_meth(Xs,Ss);
}
std::map<std::string, std::shared_ptr<Tensor> > LSTHCERI::build_X(bool meth)
{
    std::map<std::string, std::shared_ptr<Tensor> > Xs;
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        for (int l = 0; l < 3; l++) {
            std::string space = spaces[l];
            if (!Xs.count(space)) {
                int start = spaces_[space].first;
                int end   = spaces_[space].second;

                int na = end - start;
                int nn = C_->colspi()[0];
                int nm = X_->rowspi()[0];
                int nP = X_->colspi()[0];

                std::shared_ptr<Tensor> X =
                    (meth ? CoreTensor::build("T_" + space, space, na,"NP", nP) : CoreTensor::build("X_" + space, space, na,"NP", nP));

                double** X1p = X_->pointer();
                double* X2p = X->pointer();
                double** Cp  = C_->pointer();

                C_DGEMM('T','N',na,nP,nm,1.0,&Cp[0][start],nn,X1p[0],nP,0.0,X2p,nP);

                if (balance_) {
                    for (int P = 0; P < nP; P++) {
                        double w = C_DDOT(na, X2p + P, nP, X2p + P, nP);
                        C_DSCAL(na, pow(w, -1.0/2.0), X2p + P, nP);
                    }
                }

                Xs[space] = X;
            }
        }
    }
    return Xs;
}
std::map<std::string, std::shared_ptr<Tensor> > LSTHCERI::build_E(std::map<std::string, std::shared_ptr<Tensor> >& Xs)
{
    // > Unique Tasks < //

    std::map<std::string, std::pair<std::string, std::string> > pair_spaces;
    std::vector<std::string> pair_spaces_order;
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        for (int l = 0; l < 4; l+=2) {
            std::string space1 = spaces[l];
            std::string space2 = spaces[l+1];
            if (!pair_spaces.count(space1 + "_" + space2)) {
                pair_spaces_order.push_back(space1 + "_" + space2);
                pair_spaces[space1 + "_" + space2] = std::pair<std::string,std::string>(space1,space2);
            }
        }
    }

    // > Task merging < //

    std::vector<std::vector<std::string> > tasks;
    std::vector<std::string> rem = pair_spaces_order;

    // Lump all common first half-transforms together
    while (rem.size()) {
        std::vector<std::string> task;
        std::vector<std::string> rem2;
        std::string space1 = pair_spaces[rem[0]].first;
        for (int i = 0; i < rem.size(); i++) {
            if (space1 == pair_spaces[rem[i]].first)
                task.push_back(rem[i]);
            else
                rem2.push_back(rem[i]);
        }
        tasks.push_back(task);
        rem = rem2;
    }

    // > Sizing < //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int nP = X_->colspi()[0];

    // > Threading < //

    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif

    // > Maximum orbital sizes < //

    size_t max1 = 0L;
    size_t max2 = 0L;

    for (int i = 0; i < pair_spaces_order.size(); i++) {
        std::string name = pair_spaces_order[i];
        std::string space1 = pair_spaces[name].first;
        std::string space2 = pair_spaces[name].second;

        int size1 = spaces_[space1].second - spaces_[space1].first;
        int size2 = spaces_[space2].second - spaces_[space2].first;

        max1 = (max1 < size1 ? size1 : max1);
        max2 = (max2 < size2 ? size2 : max2);
    }

    // > Row requirements < //

    size_t per_row = 0L;
    // (Q|mn)
    per_row += nso * (size_t) nso;
    // (Q|mi)
    per_row += max1 * (size_t) nso;
    // E_Q^P
    per_row += nP;

    // > Overhead < //

    size_t mem = memory_;
    mem -= nthread * max2 * (size_t) nP;
    for (std::map<std::string, std::shared_ptr<Tensor> >::iterator it = Xs.begin();
        it != Xs.end(); ++it) {
        std::shared_ptr<Tensor> X = (*it).second;
        mem -= X->sizes()[0] * (size_t) X->sizes()[1];
    }

    // > C-matrix weirdness < //

    double** Cp = C_->pointer();
    int lda = C_->colspi()[0];

    // > R pre-transformed objects < //

    std::map<std::string, std::shared_ptr<Tensor> > Rs;
    for (int ind1 = 0; ind1 < tasks.size(); ind1++) {
        for (int ind2 = 0; ind2 < tasks[ind1].size(); ind2++) {
            std::string space = pair_spaces[tasks[ind1][ind2]].second;
            if (!Rs.count(space)) {

                std::shared_ptr<Tensor> X = Xs[space];

                int start1 = spaces_[space].first;
                int end1   = spaces_[space].second;
                int n1      = end1 - start1;

                std::shared_ptr<Tensor> R = CoreTensor::build("R_" + space, "NSO", nso, "NP", nP);
                double* Xp = X->pointer();
                double* Rp = R->pointer();

                C_DGEMM('N','N',nso,nP,n1,1.0,Cp[0] + start1,lda,Xp,nP,0.0,Rp,nP);
                Rs[space] = R;
                mem -= nso * (size_t) nP;
            }
        }
    }

    // > Maximum number of rows < //

    unsigned long int max_rows = (mem / per_row);
    max_rows = (max_rows < auxiliary_->max_function_per_shell() ? auxiliary_->max_function_per_shell() : max_rows);
    max_rows = (max_rows > auxiliary_->nbf() ? auxiliary_->nbf() : max_rows);

    // > Shell block assignments < //

    std::vector<int> shell_starts;
    shell_starts.push_back(0);
    int index = 0;
    for (int Qshell = 1; Qshell < auxiliary_->nshell()-1; Qshell++) {
        if (auxiliary_->shell(Qshell+1).function_index() - auxiliary_->shell(shell_starts[index]).function_index() > max_rows) {
            shell_starts.push_back(Qshell);
            index++;
        }
    }
    shell_starts.push_back(auxiliary_->nshell());

    for (int i=0; i<shell_starts.size()-1; i++) {
        if (i == shell_starts.size() - 2) {
            if (max_rows < auxiliary_->nbf() - auxiliary_->shell(shell_starts[i]).function_index()) {
              throw PSIEXCEPTION("Out of memory in DFERI.");
            }
        } else {
            if (max_rows < auxiliary_->shell(shell_starts[i+1]).function_index() - auxiliary_->shell(shell_starts[i]).function_index()) {
              throw PSIEXCEPTION("Out of memory in DFERI.");
            }
        }
    }

    // => ERI Objects <= //

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    }

    // => ERI Sieve <= //

    std::shared_ptr<ERISieve> sieve(new ERISieve(primary_, schwarz_cutoff_));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    long int nshell_pairs = (long int) shell_pairs.size();

    // => Temporary Tensors <= //

    // > Three-index buffers < //
    std::shared_ptr<Matrix> Amn(new Matrix("(A|mn)", max_rows, nso * (unsigned long int) nso));
    std::shared_ptr<Matrix> Ami(new Matrix("(A|mi)", max_rows, nso * (unsigned long int) max1));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();

    std::vector<std::shared_ptr<Matrix> > QiP;
    for (int thread = 0; thread < nthread; thread++) {
        QiP.push_back(std::shared_ptr<Matrix>(new Matrix("QiP", max1, nP)));
    }

    std::shared_ptr<Matrix> E(new Matrix("E", max_rows, nP));
    double** Ep = E->pointer();

    // > E Targets < //

    std::map<std::string, std::shared_ptr<Tensor> > Es;
    for (int k = 0; k < pair_spaces_order.size(); k++) {
        std::string name = pair_spaces_order[k];
        std::string space1 = pair_spaces[name].first;
        std::string space2 = pair_spaces[name].second;

        std::shared_ptr<Tensor> ET = DiskTensor::build("E_" + space1 + "_" + space2,
            "NAUX", naux, "NP", nP, false, false);

        Es[space1 + "_" + space2] = ET;
    }

    // ==> Master Loop <== //

    for (int block = 0; block < shell_starts.size() - 1; block++) {

        // > Block characteristics < //

        int Pstart = shell_starts[block];
        int Pstop  = shell_starts[block+1];
        int nPshell = Pstop - Pstart;
        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? auxiliary_->nbf() : auxiliary_->shell(Pstop).function_index());
        int rows = pstop - pstart;

        // > (Q|mn) ERIs < //

        ::memset((void*) Amnp[0], '\0', sizeof(double) * rows * nso * nso);

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int PMN = 0L; PMN < nPshell * nshell_pairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P  = PMN / nshell_pairs + Pstart;
            int MN = PMN % nshell_pairs;
            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            eri[thread]->compute_shell(P,0,M,N);

            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();
            int np = auxiliary_->shell(P).nfunction();
            int om = primary_->shell(M).function_index();
            int on = primary_->shell(N).function_index();
            int op = auxiliary_->shell(P).function_index();

            const double* buffer = eri[thread]->buffer();

            for (int p = 0; p < np; p++) {
            for (int m = 0; m < nm; m++) {
            for (int n = 0; n < nn; n++) {
                Amnp[p + op - pstart][(m + om) * nso + (n + on)] =
                Amnp[p + op - pstart][(n + on) * nso + (m + om)] =
                (*buffer++);
            }}}
        }

        for (int ind1 = 0; ind1 < tasks.size(); ind1++) {

            std::string space1 = pair_spaces[tasks[ind1][0]].first;
            int start1 = spaces_[space1].first;
            int end1   = spaces_[space1].second;
            int n1      = end1 - start1;
            double* C1p = &Cp[0][start1];

            std::shared_ptr<Tensor> X1 = Xs[space1];
            double* X1p = X1->pointer();

            C_DGEMM('N','N',rows*nso,n1,nso,1.0,Amnp[0],nso,C1p,lda,0.0,Amip[0],n1);

            for (int ind2 = 0; ind2 < tasks[ind1].size(); ind2++) {
                std::string space2 = pair_spaces[tasks[ind1][ind2]].second;
                int start2 = spaces_[space2].first;
                int end2   = spaces_[space2].second;
                int n2      = end2 - start2;

                std::shared_ptr<Tensor> R1 = Rs[space2];
                double* R1p = R1->pointer();

                size_t no1 = nso * (size_t) n1;

                #pragma omp parallel for num_threads(nthread)
                for (int Q = 0; Q < rows; Q++) {

                    int thread = 0;
                    #ifdef _OPENMP
                        thread = omp_get_thread_num();
                    #endif

                    double* QiPp = QiP[thread]->pointer()[0];

                    C_DGEMM('T','N',n1,nP,nso,1.0,Amip[0] + Q*no1,n1,R1p,nP,0.0,QiPp,nP);
                    for (int P = 0; P < nP; P++) {
                        Ep[Q][P] = C_DDOT(n1,X1p + P,nP,QiPp + P,nP);
                    }
                }

                std::string name = tasks[ind1][ind2];
                std::shared_ptr<Tensor> A = Es[name];
                FILE* fh = A->file_pointer();
                fwrite(Ep[0],sizeof(double),rows*nP,fh);
            }
        }
    }

    return Es;
}
std::map<std::string, std::shared_ptr<Tensor> > LSTHCERI::build_S(std::map<std::string, std::shared_ptr<Tensor> >& Xs, bool meth)
{
    std::map<std::string, std::shared_ptr<Tensor> > Ss;
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        for (int l = 0; l < 4; l+=2) {
            std::string space1 = spaces[l];
            std::string space2 = spaces[l+1];
            if (!Ss.count(space1 + "_" + space2)) {
                std::shared_ptr<Tensor> X1 = Xs[space1];
                std::shared_ptr<Tensor> X2 = Xs[space2];
                int nP = X1->sizes()[1];
                int n1 = X1->sizes()[0];
                int n2 = X2->sizes()[0];

                std::shared_ptr<Matrix> S1(new Matrix("S1", nP, nP));
                std::shared_ptr<Tensor> S2 = (meth ?
                    CoreTensor::build("STinv_" + space1 + "_" + space2, "NP", nP, "NP", nP) :
                    CoreTensor::build("Sinv_" + space1 + "_" + space2, "NP", nP, "NP", nP));

                double*  X1p = X1->pointer();
                double*  X2p = X2->pointer();
                double*  S1p = S1->pointer()[0];
                double*  S2p  = S2->pointer();

                C_DGEMM('T','N',nP,nP,n1,1.0,X1p,nP,X1p,nP,0.0,S1p,nP);
                C_DGEMM('T','N',nP,nP,n2,1.0,X2p,nP,X2p,nP,0.0,S2p,nP);

                double* S1tp = S1p;
                double* S2tp = S2p;
                for (size_t ind = 0L; ind < nP * (size_t) nP; ind++) {
                    (*S1tp++) *= (*S2tp++);
                }

                S1->hermitivitize();
                S1->power(-1.0,S_cutoff_);

                ::memcpy(S2p,S1p,sizeof(double) * nP * nP);

                S2->swap_out();
                Ss[space1 + "_" + space2] = S2;
            }
        }
    }
    return Ss;
}
std::map<std::string, std::shared_ptr<Tensor> > LSTHCERI::build_L(std::map<std::string, std::shared_ptr<Tensor> >& Es,
                                                                    std::map<std::string, std::shared_ptr<Tensor> >& Ss)
{
    std::shared_ptr<Matrix> J = Jm12(auxiliary_,J_cutoff_);
    double* Jp = J->pointer()[0];

    std::map<std::string, std::shared_ptr<Tensor> > Ls;
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        for (int l = 0; l < 4; l+=2) {
            std::string space1 = spaces[l];
            std::string space2 = spaces[l+1];
            if (!Ls.count(space1 + "_" + space2)) {
                std::shared_ptr<Tensor> E = Es[space1 + "_" + space2];
                std::shared_ptr<Tensor> S = Ss[space1 + "_" + space2];
                S->swap_in();

                int nA = E->sizes()[0];
                int nP = E->sizes()[1];

                std::shared_ptr<Tensor> L = CoreTensor::build("L_" + space1 + "_" + space2,
                    "NP", nP, "NAUX", nA);
                std::shared_ptr<Matrix> T(new Matrix("LT", nA, nP));

                double* Tp = T->pointer()[0];
                double* Lp = L->pointer();
                double* Sp = S->pointer();
                FILE* fh = E->file_pointer();

                // Avert your eyes
                fseek(fh,0,SEEK_SET);
                size_t statusvalue=fread(Lp,sizeof(double),nA*(size_t)nP,fh);

                C_DGEMM('N','N',nA,nP,nA,1.0,Jp,nA,Lp,nP,0.0,Tp,nP);
                C_DGEMM('N','T',nP,nA,nP,1.0,Sp,nP,Tp,nP,0.0,Lp,nA);

                S->swap_out();
                L->swap_out();
                Ls[space1 + "_" + space2] = L;
            }
        }
    }
    return Ls;
}
std::map<std::string, std::shared_ptr<Tensor> > LSTHCERI::build_Z(std::map<std::string, std::shared_ptr<Tensor> >& Ls)
{
    std::map<std::string, std::shared_ptr<Tensor> > Zs;
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        std::shared_ptr<Tensor> L12 = Ls[spaces[0] + "_" + spaces[1]];
        std::shared_ptr<Tensor> L34 = Ls[spaces[2] + "_" + spaces[3]];
        int nP = L12->sizes()[0];
        int nA = L12->sizes()[1];
        std::shared_ptr<Tensor> Z = CoreTensor::build("Z_" + name,
            "NP", nP, "NP", nP);
        L12->swap_in();
        L34->swap_in();
        double* L12p = L12->pointer();
        double* L34p = L34->pointer();
        double* Zp = Z->pointer();
        C_DGEMM('N','T',nP,nP,nA,1.0,L12p,nA,L34p,nA,0.0,Zp,nP);
        L12->swap_out();
        L34->swap_out();
        Z->swap_out();
        Zs[name] = Z;
    }
    return Zs;
}
void LSTHCERI::pack(std::map<std::string, std::shared_ptr<Tensor> >& Xs,
                    std::map<std::string, std::shared_ptr<Tensor> >& Zs,
                    std::map<std::string, std::shared_ptr<Tensor> >& Ls,
                    std::map<std::string, std::shared_ptr<Tensor> >& Ss)
{
    ints_.clear();
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        std::vector<std::shared_ptr<Tensor> > task;
        task.push_back(Xs[spaces[0]]);
        task.push_back(Xs[spaces[1]]);
        task.push_back(Zs[name]);
        task.push_back(Xs[spaces[2]]);
        task.push_back(Xs[spaces[3]]);
        task.push_back(Ls[spaces[0] + "_" + spaces[1]]);
        task.push_back(Ls[spaces[2] + "_" + spaces[3]]);
        task.push_back(Ss[spaces[0] + "_" + spaces[1]]);
        task.push_back(Ss[spaces[2] + "_" + spaces[3]]);
        ints_[name] = task;
    }
}
void LSTHCERI::pack_meth(std::map<std::string, std::shared_ptr<Tensor> >& Xs,
                         std::map<std::string, std::shared_ptr<Tensor> >& Ss)
{
    meths_.clear();
    for (int k = 0; k < eri_spaces_order_.size(); k++) {
        std::string name = eri_spaces_order_[k];
        std::vector<std::string> spaces = eri_spaces_[name];
        std::vector<std::shared_ptr<Tensor> > task;
        task.push_back(Xs[spaces[0]]);
        task.push_back(Xs[spaces[1]]);
        task.push_back(Ss[spaces[0] + "_" + spaces[1]]);
        meths_[name] = task;
    }
}

} // Namespace
