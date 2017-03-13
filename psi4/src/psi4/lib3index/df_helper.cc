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

#include "df_helper.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/sieve.h"
#include "psi4/lib3index/dftensor.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace psi{ namespace df_helper{

DF_Helper::DF_Helper(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux)
{
    primary_ = primary;
    aux_     = aux;

    nao_  = primary_->nbf();
    naux_ = aux_    ->nbf();

    prepare_blocking();
}
DF_Helper::~DF_Helper()
{
    // close streams
    for(auto& kv : file_status_)
        fclose(kv.second.fp);

    // destroy all files!
    std::string file1, file2;
    for(auto& kv : files_){
        remove(std::get<0>(kv.second).c_str());
        remove(std::get<1>(kv.second).c_str());
    }
    for(auto& kv : AO_files_)
        remove(kv.second.c_str());
}
void DF_Helper::prepare_blocking()
{
    Qshells_ = aux_->nshell();
    pshells_ = primary_->nshell();

    Qshell_aggs_.reserve(Qshells_+1);
    pshell_aggs_.reserve(pshells_+1);

    // Aux shell blocking
    Qshell_aggs_[0] = 0;
    for(size_t i=0; i<Qshells_; i++)
        Qshell_aggs_[i+1] = Qshell_aggs_[i]+aux_->shell(i).nfunction();

    // AO shell blocking
    pshell_aggs_[0] = 0;
    for(size_t i=0; i<pshells_; i++)
        pshell_aggs_[i+1] = pshell_aggs_[i]+primary_->shell(i).nfunction();
}
void DF_Helper::AO_filename_maker(size_t i)
{
    std::string name = PSIOManager::shared_object()->get_default_path();
    name.append("dfh.AO");
    name.append(std::to_string(i));
    AO_names_.push_back(name);
    std::string file = name;
    file.append(".");
    file.append(std::to_string(getpid()));
    file.append(".");
    file.append(primary_->molecule()->name());
    file.append(".dat");
    AO_files_[name] = file;
}
void DF_Helper::filename_maker(std::string name, size_t a0, size_t a1, size_t a2)
{
    std::string pfile = "dfh.p";
    pfile.append(name);
    pfile.append(".");
    pfile.append(std::to_string(getpid()));
    pfile.append(".");
    pfile.append(primary_->molecule()->name());
    pfile.append(".dat");
    std::string file = pfile;
    file.erase(4, 1);

    std::string scratch = PSIOManager::shared_object()->get_default_path();
    std::string pfilename = scratch;
    pfilename.append(pfile);
    std::string  filename = scratch;
    filename.append(file);

    std::tuple<std::string, std::string> files(pfilename.c_str(), filename.c_str());
    files_[name] = files;

    std::tuple<size_t, size_t, size_t> sizes = std::make_tuple(a0, a1, a2);
    sizes_[pfilename]=sizes;
    sizes_[ filename]=sizes;
}
void DF_Helper::initialize()
{
    outfile->Printf("\n\n    ==> DF_Helper: initialize() <==\n\n");
    timer_on("DF_Helper~-- build-------   ");
    if(method_.compare("DIRECT") && method_.compare("STORE")){
        std::stringstream error;
        error << "DF_Helper:initialize: specified method (" << method_ << ") is incorrect";
        throw PSIEXCEPTION(error.str().c_str());
    }

    outfile->Printf("      nao           = %zu\n", nao_);
    outfile->Printf("      naux          = %zu\n", naux_);
    outfile->Printf("      Scwarz cutoff = %E\n", cutoff_);
    outfile->Printf("      mem (doubles) = %zu\n", memory_);
    outfile->Printf("      method        = %s\n", method_.c_str());
    outfile->Printf("      core?         = %d\n", (int)core_);
    outfile->Printf("      hold met?     = %d\n", (int)hold_met_);

    direct_ = (!method_.compare("DIRECT") ? true : false);
    prepare_sparsity();
    //print_masks();

    if(!(std::fabs(mpower_ - 0.0) < 1e-13)){
        if(core_)
            prepare_metric_core();
        else{
            if(!hold_met_ || !direct_)
                prepare_metric();
        }
    }

    if(core_){
        Ppq_.reserve(big_skips_[nao_]);
        prepare_AO_core();
    }
    else{
        if(!direct_){
            plargest_ = shell_blocks(memory_, 0, 1, pshells_, psteps_);
            Qlargest_ = shell_blocks(memory_, 0, 0, Qshells_, Qsteps_);
            prepare_AO();
        }
    }
    built = true;
    timer_off("DF_Helper~-- build-------   ");
    outfile->Printf("\n\n    ==> DF_Helper: initialize() <==\n\n");
}
void DF_Helper::prepare_sparsity()
{
    std::vector<double> fun_prints(nao_*nao_, 0.0);
    std::vector<double> shell_prints(pshells_*pshells_, 0.0);

    schwarz_fun_mask_.reserve(nao_*nao_);
    schwarz_shell_mask_.reserve(pshells_*pshells_);
    small_skips_.reserve(nao_+1);
    big_skips_.reserve(nao_+1);

    // prepare eri buffers
    int rank = 0;
    IntegralFactory schwarzfactory(primary_,primary_,primary_,primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    std::vector<const double*> buffer(nthreads_);
    #pragma omp parallel num_threads(nthreads_) private(rank) // ease NUMA
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
        buffer[rank] = eri[rank]->buffer();
    }

    double val;
    size_t MU, NU, mu, nu, omu, onu, nummu, numnu, index;
    #pragma omp parallel for private(MU, NU, mu, nu, omu, onu, nummu, numnu, index, val, rank) num_threads(nthreads_) schedule(guided)
    for (MU=0; MU < pshells_; ++MU) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            eri[rank]->compute_shell(MU,NU,MU,NU);
            for (mu=0; mu < nummu; ++mu) {
                omu = primary_->shell(MU).function_index() + mu;
                for (nu=0; nu < numnu; ++nu) {
                    onu = primary_->shell(NU).function_index() + nu;
                    if (omu>=onu) {
                        index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                        val = fabs(buffer[rank][index]);
                        if (shell_prints[MU*pshells_+NU] <= val){
                            shell_prints[MU*pshells_+NU]=val;
                            shell_prints[NU*pshells_+MU]=val;
                        }
                        if (fun_prints[omu*nao_+onu] <= val)
                        {
                            fun_prints[omu*nao_+onu]= val;
                            fun_prints[onu*nao_+omu]= val;
                        }
                    }
                }
            }
        }
    }

//    #pragma omp parallel for simd num_threads(nthreads_) schedule(static)
    for(size_t i=0; i<pshells_*pshells_; i++)
        schwarz_shell_mask_[i] = (shell_prints[i]<cutoff_ ? 0 : 1);

    size_t count;
    #pragma omp parallel for private(count) num_threads(nthreads_)
    for(size_t i=0; i<nao_; i++){
        count=0;
        for(size_t j=0; j<nao_; j++){
            if(fun_prints[i*nao_+j]>cutoff_){
                count++;
                schwarz_fun_mask_[i*nao_+j]=count;
            }
            else
                schwarz_fun_mask_[i*nao_+j]=0;
        }
        small_skips_[i] = count;
    }

    // build indexing skips
    big_skips_[0]=0;
    size_t coltots=0;
    for(size_t j=0; j<nao_; j++){
        size_t cols = small_skips_[j];
        size_t size = cols*naux_;
        coltots += cols;
        big_skips_[j+1] = size+big_skips_[j];
    }
    small_skips_[nao_]=coltots;
}
void DF_Helper::print_masks()
{
    outfile->Printf("\n----fun mask------");
    for(size_t i=0; i<nao_; i++){
        outfile->Printf("\n");
        for(size_t j=0; j<nao_; j++){
            outfile->Printf("%zu, ", schwarz_fun_mask_[i*nao_+j]);
        }
        outfile->Printf("  count: %zu", small_skips_[i]);
    }
    //outfile->Printf("\n----shell mask----");
    //for(i=0; i<pshell_; i++){
    //    outfile->Printf("\n");
    //    for(j=0; j<pshell_; j++){
    //        outfile->Printf("%d, ", schwarz_shell_mask_[i*pshell_+j]);
    //    }
    //}
    outfile->Printf("\n");
return;
}
void DF_Helper::prepare_AO()
{
    // prepare eri buffers -- check here for size?
    size_t rank = 0;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    #pragma omp parallel num_threads(nthreads_) private(rank)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // declare largest necessary according to plargest_
    std::vector<double> M; M.reserve(std::get<0>(plargest_));
    double* Mp = M.data();

    // prepare files
    AO_filename_maker(1);
    AO_filename_maker(2);
    std::string putf = AO_files_[AO_names_[0]];

    // prep stream
    stream_check(putf, "ab");

    size_t jump = 0;
    size_t steps = psteps_.size();
    timer_on("DF_Helper~AO construction   ");
    outfile->Printf("\n    ==> Begin AO Blocked Construction <==\n\n");
    for(size_t i=0; i<steps; i++){

        size_t start = std::get<0>(psteps_[i]);
        size_t stop  = std::get<1>(psteps_[i]);
        size_t begin = pshell_aggs_[start];
        size_t end = pshell_aggs_[stop+1]-1;
        size_t block_size = end-begin+1;

        // setup
        size_t size = big_skips_[end+1]-big_skips_[begin];

        // compute
        compute_AO_p(start, stop, Mp, eri);

        // put
        put_tensor_AO(putf, Mp, size, jump, "ab");
        jump += size;
    }
    outfile->Printf("\n    ==> End AO Blocked Construction <==");
    timer_off("DF_Helper~AO construction   ");
    timer_on("DF_Helper~metric contraction");
    contract_metric_AO(Mp);
    timer_off("DF_Helper~metric contraction");
}
void DF_Helper::prepare_AO_core()
{
    size_t rank = 0;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    #pragma omp parallel num_threads(nthreads_) private(rank)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    if(!direct_){
        std::vector<double> Qpq;
        Qpq.reserve(big_skips_[nao_]);
        timer_on("DF_Helper~AO construction   ");
        compute_AO_p(0, pshells_-1, &Qpq[0], eri);
        timer_off("DF_Helper~AO construction   ");
        double* metp = metric_prep_core(mpower_);
        timer_on("DF_Helper~metric contraction");
        contract_metric_AO_core(&Qpq[0], metp);
        timer_off("DF_Helper~metric contraction");
    }
    else{
        timer_on("DF_Helper~AO construction   ");
        compute_AO_p(0, pshells_-1, &Ppq_[0], eri);
        timer_off("DF_Helper~AO construction   ");
    }
}
std::pair<size_t, size_t> DF_Helper::shell_blocks(const size_t mem, size_t max,
  bool op, const size_t shell_tots, std::vector<std::pair<size_t, size_t>>& b)
{
    size_t extra = (op ? naux_*naux_ : 0);
    if(wMO_==0)
        wMO_ = nao_/2;
    if(max==0)
        max = wMO_*wMO_;

    size_t current, block_size, tmpbs=0, total=0, count=0, largest=0;
    for(size_t i=0; i<shell_tots; i++){

        count++;
        current=0;
        if(op){
            size_t begin = pshell_aggs_[i];
            size_t end   = pshell_aggs_[i+1]-1;
            tmpbs += end-begin+1;
            current = big_skips_[end+1]-big_skips_[begin];
        }
        else{
            size_t begin = Qshell_aggs_[i];
            size_t end   = Qshell_aggs_[i+1]-1;
            tmpbs += end-begin+1;
            current = (end-begin+1)*small_skips_[nao_];
        }
        total += current;
        if((op ? 2*total : total + (wMO_*nao_ +  2*max)*tmpbs) + extra > mem || i==shell_tots-1){
            if(count==1 && i!=shell_tots-1){
                std::stringstream error;
                error << "DF_Helper: not enough memory for AO integral blocking!";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if(i == shell_tots-1)
                b.push_back(std::make_pair(i-count+1, i));
            else{
                total -= current;
                b.push_back(std::make_pair(i-count+1, i-1));
                i--;
            }
            if(largest<total){
                largest=total;
                block_size=tmpbs;
            }
            count=0;
            total=0;
            tmpbs=0;
        }
    }
// returns pair(largest buffer size, largest block size)
return std::make_pair(largest, block_size);
}
FILE* DF_Helper::stream_check(std::string filename, std::string op)
{
    timer_on("stream checks               ");
    if(file_status_.find(filename) == file_status_.end()){
        //outfile->Printf("file not found: %s\n", filename.c_str());
        file_status_[filename].op = op;
        file_status_[filename].fp = fopen(filename.c_str(), op.c_str());
    }
    else{
        if(op.compare(file_status_[filename].op)){
    timer_on("stream change               ");
            //outfile->Printf("changing......: %s\n", filename.c_str());
            file_status_[filename].op = op;
            fflush(file_status_[filename].fp);
            fclose(file_status_[filename].fp);
            file_status_[filename].fp = fopen(filename.c_str(), op.c_str());
    timer_off("stream change               ");
        }
    }
    timer_off("stream checks               ");
return file_status_[filename].fp;
}
void DF_Helper::put_tensor(std::string file, double* b, std::pair<size_t, size_t> i0,
  std::pair<size_t, size_t> i1, std::pair<size_t, size_t> i2, std::string op)
{
    // collapse to 2D, assume file has form (i1 | i2 i3)
    size_t A2 = std::get<2>(sizes_[file]);

    size_t sta0 = std::get<0>(i0);
    size_t sto0 = std::get<1>(i0);
    size_t sta1 = std::get<0>(i1);
    size_t sto1 = std::get<1>(i1);
    size_t sta2 = std::get<0>(i2);
    size_t sto2 = std::get<1>(i2);

    size_t a0 = sto0-sta0+1;
    size_t a1 = sto1-sta1+1;
    size_t a2 = sto2-sta2+1;

    // check contiguity (a2)
    if(A2==a2){
        put_tensor(file, b, sta0, sto0, a2*sta1, a2*(sto1+1)-1, op);
    }
    else{ // loop (a0, a1)
        for(size_t j=0; j<a0; j++){
            for(size_t i=0; i<a1; i++){
                put_tensor(file, &b[j*(a1*a2)+i*a2], sta0+j, sta0+j,
                  (i+sta1)*A2+sta2, (i+sta1)*A2+sta2+a2-1, op);
            }
        }
    }
}
void DF_Helper::put_tensor(std::string file, double* Mp, const size_t start1,
 const size_t stop1, const size_t start2, const size_t stop2, std::string op)
{
    size_t a0 = stop1-start1+1;
    size_t a1 = stop2-start2+1;
    size_t A0 = std::get<0>(sizes_[file]);
    size_t A1 = std::get<1>(sizes_[file])*std::get<2>(sizes_[file]);
    size_t st = A1-a1;

    // begin stream
    FILE* fp = stream_check(file, op);

    // adjust position
    fseek(fp, (start1*A0+start2)*sizeof(double), SEEK_SET);

    // is everything contiguous?
    if(st==0){
        size_t s = fwrite(&Mp[0], sizeof(double), a0*a1, fp);
        if(!s){
            std::stringstream error;
            error << "DF_Helper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
    else
    {
        for(size_t i=start1; i<stop1; i++){
            // write
            size_t s = fwrite(&Mp[i*a1], sizeof(double), a1, fp);
            if(!s){
                std::stringstream error;
                error << "DF_Helper:put_tensor: write error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            fseek(fp, st*sizeof(double), SEEK_CUR);
        }
        // manual last one
        size_t s = fwrite(&Mp[(a0-1)*a1], sizeof(double), a1, fp);
        if(!s){
            std::stringstream error;
            error << "DF_Helper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}
void DF_Helper::put_tensor_AO(std::string file, double* Mp, size_t size, size_t start, std::string op)
{
    // begin stream - don't stream check because i knew beforehand!
    FILE* fp = file_status_[file].fp;

    // adjust position
    fseek(fp, start, SEEK_SET);

    // everything is contiguous
    size_t s=fwrite(&Mp[0], sizeof(double), size, fp);
    if(!s){
        std::stringstream error;
        error << "DF_Helper:put_tensor_AO: write error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::get_tensor_AO(std::string file, double* Mp, size_t size, size_t start)
{
    // begin stream
    FILE* fp = file_status_[file].fp;

    // adjust position
    fseek(fp, start*sizeof(double), SEEK_SET);

    // everything is contiguous
    size_t s = fread(&Mp[0], sizeof(double), size, fp);
    if(!s){
        std::stringstream error;
        error << "DF_Helper:get_tensor_AO: read error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::get_tensor_(std::string file, double* b, std::pair<size_t, size_t> i0,
  std::pair<size_t, size_t> i1, std::pair<size_t, size_t> i2)
{
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(file) != tsizes_.end() ? tsizes_[file] : sizes_[file]);

    // collapse to 2D, assume file has form (i1 | i2 i3)
    size_t A2 = std::get<2>(sizes);

    size_t sta0 = std::get<0>(i0);
    size_t sto0 = std::get<1>(i0);
    size_t sta1 = std::get<0>(i1);
    size_t sto1 = std::get<1>(i1);
    size_t sta2 = std::get<0>(i2);
    size_t sto2 = std::get<1>(i2);

    size_t a0 = sto0-sta0+1;
    size_t a1 = sto1-sta1+1;
    size_t a2 = sto2-sta2+1;

    // check contiguity (a2)
    if(A2==a2){
        get_tensor_(file, b, sta0, sto0, a2*sta1, a2*(sto1+1)-1);
    }
    else{ // loop (a0, a1)
        for(size_t j=0; j<a0; j++){
            for(size_t i=0; i<a1; i++){
                get_tensor_(file, &b[j*(a1*a2)+i*a2], sta0+j, sta0+j,
                  (i+sta1)*A2+sta2, (i+sta1)*A2+sta2+a2-1);
            }
        }
    }
}
void DF_Helper::get_tensor_(std::string file, double* b, const size_t start1,
 const size_t stop1, const size_t start2, const size_t stop2)
{
    size_t a0 = stop1-start1+1;
    size_t a1 = stop2-start2+1;

    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(file) != tsizes_.end() ? tsizes_[file] : sizes_[file]);

    size_t A0 = std::get<0>(sizes);
    size_t A1 = std::get<1>(sizes)*std::get<2>(sizes);
    size_t st = A1-a1;

    // check stream
    FILE* fp = stream_check(file, "rb");

    // adjust position
    fseek(fp, (start1*A1+start2)*sizeof(double), SEEK_SET);

    // is everything contiguous?
    if(st==0){
        size_t s = fread(&b[0], sizeof(double), a0*a1, fp);
        if(!s){
            std::stringstream error;
            error << "DF_Helper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
    else
    {
        for(size_t i=0; i<a0-1; i++){
            // read
            size_t s = fread(&b[i*a1], sizeof(double), a1, fp);
            if(!s){
                std::stringstream error;
                error << "DF_Helper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            s =fseek(fp, st*sizeof(double), SEEK_CUR);
            if(s){
                std::stringstream error;
                error << "DF_Helper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
        }
        // manual last one
        size_t s = fread(&b[(a0-1)*a1], sizeof(double), a1, fp);
        if(!s){
            std::stringstream error;
            error << "DF_Helper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}
void DF_Helper::compute_AO_Q(const size_t start, const size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri)
{
    size_t begin = Qshell_aggs_[start];
    size_t   end = Qshell_aggs_[stop+1]-1;
    size_t block_size = end-begin+1;

    // prepare eri buffers
    size_t nthread = nthreads_;
    if(eri.size() != nthreads_)
        nthread = eri.size();

    int rank = 0;
    std::vector<const double*> buffer(nthread);
    #pragma omp parallel private(rank) num_threads(nthread)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    size_t MU, nummu, NU, numnu, Pshell, numP, mu, omu, nu, onu, P, PHI;
    #pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule(guided) num_threads(nthreads_)
    for (MU=0; MU<pshells_; MU++){
       #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_ -> shell(MU).nfunction();
        for (NU=0; NU<pshells_; NU++){
            numnu = primary_ -> shell(NU).nfunction();
            if(!schwarz_shell_mask_[MU*pshells_+NU]){continue;}
            for (Pshell=start; Pshell<=stop; Pshell++){
                PHI = aux_->shell(Pshell).function_index();
                numP = aux_ -> shell(Pshell).nfunction();
                eri[rank] -> compute_shell(Pshell, 0, MU, NU);
                for (mu=0; mu<nummu; mu++){
                    omu = primary_ -> shell(MU).function_index() + mu;
                    for(nu=0; nu<numnu; nu++){
                        onu = primary_->shell(NU).function_index() + nu;
                        if(!schwarz_fun_mask_[omu*nao_+onu]){continue;}
                        for(P=0; P<numP; P++){
                            Mp[(big_skips_[omu]*block_size)/naux_
                              +(PHI+P-begin)*small_skips_[omu]+schwarz_fun_mask_[omu*nao_+onu]-1]
                              = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DF_Helper::compute_AO_p(const size_t start, const size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri)
{
    size_t begin = pshell_aggs_[start];
    size_t   end = pshell_aggs_[stop+1]-1;
    size_t block_size = end-begin+1;
    size_t startind = big_skips_[begin];
    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    outfile->Printf(", nao index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if(eri.size() != nthreads_)
        nthread = eri.size();

    int rank = 0;
    std::vector<const double*> buffer(nthread);
    #pragma omp parallel private(rank) num_threads(nthread)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    size_t MU, nummu, NU, numnu, Pshell, numP, mu, omu, nu, onu, P, PHI;
    #pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule(guided) num_threads(nthread)
    for (MU=start; MU<=stop; MU++){
       #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU<pshells_; NU++){
            numnu = primary_->shell(NU).nfunction();
            if(!schwarz_shell_mask_[MU*pshells_+NU]) {continue;}
            for (Pshell=0; Pshell<Qshells_; Pshell++){
                PHI = aux_->shell(Pshell).function_index();
                numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                for (mu=0; mu<nummu; mu++){
                    omu = primary_ -> shell(MU).function_index() + mu;
                    for(nu=0; nu<numnu; nu++){
                        onu = primary_->shell(NU).function_index() + nu;
                        if(!schwarz_fun_mask_[omu*nao_+onu]) {continue;}
                        for(P=0; P<numP; P++){
                            Mp[big_skips_[omu]-startind+(PHI+P)*small_skips_[omu]
                              + schwarz_fun_mask_[omu*nao_+onu]-1]
                              = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DF_Helper::grab_AO(const size_t start, const size_t stop, double* Mp){

    size_t begin = Qshell_aggs_[start];
    size_t   end = Qshell_aggs_[stop+1]-1;
    size_t block_size = end-begin+1;
    std::string getf = AO_files_[AO_names_[1]];

    size_t sta=0;
    for(size_t i=0; i<nao_; i++){
        size_t size = block_size*small_skips_[i];
        size_t jump = begin*small_skips_[i];
        get_tensor_AO(getf, &Mp[sta], size, big_skips_[i]+jump);
        sta+=size;
    }
}
void DF_Helper::prepare_metric_core()
{
    timer_on("DF_Helper~metric construct  ");
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(aux_, true));
    Jinv->form_fitting_metric();
    metrics_[1.0] = Jinv->get_metric();
    timer_off("DF_Helper~metric construct  ");
}
double* DF_Helper::metric_prep_core(double pow)
{
    bool on = false;
    double power;
    for(auto& kv : metrics_){
        if(!(std::fabs(pow-kv.first) > 1e-13)){
            on = true;
            power = kv.first;
            break;
        }
    }
    if(!on){
        power = pow;
        timer_on("DF_Helper~metric power      ");
        SharedMatrix J = metrics_[1.0];
        J->power(power, 1e-12);
        metrics_[power] = J;
        timer_off("DF_Helper~metric power      ");
    }
return metrics_[power]->pointer()[0];
}
void DF_Helper::prepare_metric()
{
    // construct metric
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(aux_, true));
    Jinv->form_fitting_metric();
    SharedMatrix metric = Jinv->get_metric();
    double* Mp = metric->pointer()[0];

    // create file
    std::string filename = "metric";
    filename.append(".");
    filename.append(std::to_string(1.0));
    filename_maker(filename, naux_, naux_, 1);
    metric_keys_.push_back(std::make_pair(1.0, filename));

    // store
    std::string putf = std::get<0>(files_[filename]);
    put_tensor(putf, Mp, 0, naux_-1, 0, naux_-1, "wb");
}
std::string DF_Helper::return_metfile(double pow)
{
    bool on = 0;
    std::string key;
    for(size_t i=0; i<metric_keys_.size() && !on; i++){
        double pos = std::get<0>(metric_keys_[i]);
        if(std::fabs(pos-pow)<1e-12){
            key = std::get<1>(metric_keys_[i]);
            on = 1;
        }
    }

    if(!on)
        key = compute_metric(pow);
return key;
}
std::string DF_Helper::compute_metric(double pow)
{
    // ensure J
    if(std::fabs(pow - 1.0) < 1e-13)
        prepare_metric();
    else{
        // get metric
        SharedMatrix metric(new Matrix("met", naux_, naux_));
        double* metp = metric->pointer()[0];
        std::string filename = return_metfile(1.0);

        // get and compute
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_-1, 0, naux_-1);
        metric->power(pow);

        // make new file
        std::string name = "metric";
        name.append(".");
        name.append(std::to_string(pow));
        filename_maker(name, naux_, naux_, 1);
        metric_keys_.push_back(std::make_pair(pow, name));

        // store
        std::string putf = std::get<0>(files_[name]);
        put_tensor(putf, metp, 0, naux_-1, 0, naux_-1, "wb");
    }
return return_metfile(pow);
}
void DF_Helper::contract_metric(std::string file, double* metp, double* Mp, double* Fp, const size_t tots)
{
    std::string getf = std::get<0>(files_[file]);
    std::string putf = std::get<1>(files_[file]);

    size_t a0 = std::get<0>(sizes_[getf]);
    size_t a1 = std::get<1>(sizes_[getf]);
    size_t a2 = std::get<2>(sizes_[getf]);

    size_t count = 0;
    size_t maxim = 0;

    // determine blocking
    std::vector<std::pair<size_t, size_t>> steps;
    for(size_t i=0; i<a1; i++){
        count++;
        if(tots<count*a0*a2 || i==a1-1){
            if(count==1 && i!=a1-1){
                std::stringstream error;
                error << "DF_Helper:contract_metric: not enough memory.";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if(maxim<count)
                maxim = count;

            if(i == a1-1)
                steps.push_back(std::make_pair(i-count+1, i));
            else{
                steps.push_back(std::make_pair(i-count+1, i-1));
                i--;
            }
            count=0;
        }
    }

    std::string op = "wb";
    // contract in steps
    for(size_t i=0; i<steps.size(); i++){

        size_t begin = std::get<0>(steps[i]);
        size_t end   = std::get<1>(steps[i]);
        size_t bs = end-begin + 1;

        // get
        get_tensor_(getf, Mp, 0, a0-1, begin*a2, (end+1)*a2-1);

        // contract
        C_DGEMM('N', 'N', a0, bs*a2, a0, 1.0, metp, a0, Mp, bs*a2, 0.0, Fp, bs*a2);

        // put
        put_tensor(putf, Fp, 0, a0-1, begin*a2, (end+1)*a2-1, op);
    }
}
void DF_Helper::contract_metric_AO(double* Mp)
{
    // setup
    std::string getf = AO_files_[AO_names_[0]];
    std::string putf = AO_files_[AO_names_[1]];

    // prep stream
    stream_check(getf, "rb");
    stream_check(putf, "ab");

    // declare largest necessary
    std::vector<double> F; F.reserve(std::get<0>(plargest_));
    double* Fp = F.data();

    // grab metric
    std::string filename = return_metfile(mpower_);
    std::vector<double> metric; metric.reserve(naux_*naux_);
    double* metp = metric.data();
    get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_-1, 0, naux_-1);

    // Contract metric according to previously calculated scheme
    size_t count=0;
    for(size_t i=0; i<psteps_.size(); i++){

        size_t start = std::get<0>(psteps_[i]);
        size_t stop  = std::get<1>(psteps_[i]);
        size_t begin = pshell_aggs_[start];
        size_t end   = pshell_aggs_[stop+1]-1;
        size_t block_size = end-begin+1;
        size_t size = big_skips_[end+1]-big_skips_[begin];

        // grab
        get_tensor_AO(getf, Mp, size, count);

        // loop and contract
        #pragma omp parallel for num_threads(nthreads_) schedule(guided)
        for(size_t j=0; j<block_size; j++){
            size_t mi = small_skips_[begin+j];
            size_t skips = big_skips_[begin+j]-big_skips_[begin];
            C_DGEMM('N', 'N', naux_, mi, naux_, 1.0,
              metp, naux_, &Mp[skips], mi, 0.0, &Fp[skips], mi);
        }

        // put
        put_tensor_AO(putf, Fp, size, count, "ab");
        count += size;
    }
}
void DF_Helper::contract_metric_AO_core(double* Qpq, double* metp)
{
    // loop and contract
    #pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for(size_t j=0; j<nao_; j++){
        size_t mi = small_skips_[j];
        size_t skips = big_skips_[j];
        C_DGEMM('N', 'N', naux_, mi, naux_, 1.0,
          metp, naux_, &Qpq[skips], mi, 0.0, &Ppq_[skips], mi);
    }
}
void DF_Helper::add_space(std::string key, SharedMatrix M)
{
    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if(!built){
        throw PSIEXCEPTION("DF_Helper:add_space: call initialize() before adding spaces!");
    }
    else if(a0!=nao_){
        std::stringstream error;
        error << "DF_Helper:add_space: illegal space ("<< key <<"), primary axis is not nao";
        throw PSIEXCEPTION(error.str().c_str());
    }
    else if(spaces_.find(key) != spaces_.end()){
        if(a1 != std::get<1>(spaces_[key])){
            std::stringstream error;
            error << "DF_Helper:add_space: illegal space ("<< key <<"), new space has incorrect dimension!";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
    else if(wMO_ < a1 && !direct_){
            std::stringstream error;
            error << "DF_Helper:add_space: illegal space ("<< key <<"), new space is larger than the "   <<
            "worst MO size -(" << wMO_ << "<" << a1 << ")- specified at the time when initialize() was " <<
            "called, use set_MO_hint() before calling initialize()";
            throw PSIEXCEPTION(error.str().c_str());
    }
    else
        sorted_spaces_.push_back(std::make_pair(key, a1));

    spaces_[key] = std::make_tuple(M, a1);
}
void DF_Helper::add_transformation(std::string name, std::string key1, std::string key2, std::string order)
{
    if(spaces_.find(key1) == spaces_.end()){
        std::stringstream error;
        error << "DF_Helper:add_transformation: first space ("<< key1 <<"), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    else if(spaces_.find(key2) == spaces_.end()){
        std::stringstream error;
        error << "DF_Helper:add_transformation: second space ("<< key2 <<"), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    int op;
    if(!order.compare("Qpq"))
        op = 0;
    else if(!order.compare("pqQ"))
        op = 1;
    else if(!order.compare("pQq"))
        op = 2;
    else{
        throw PSIEXCEPTION("Matt doesnt do exceptions.");
    }

    size_t a1 = std::get<1>(spaces_[key1]);
    size_t a2 = std::get<1>(spaces_[key2]);
    filename_maker(name, naux_, a1, a2);
    transf_[name] = std::make_tuple(key1, key2, order);
}
void DF_Helper::clear()
{
    outfile->Printf("\n clearing everything (spaces, integrals)! \n");

    // clear spaces
    spaces_.clear();
    sorted_spaces_.clear();
    order_.clear();
    bspace_.clear();
    strides_.clear();

    // no ordering
    ordered_ = false;
    once_ = false;

    std::string left;
    std::string right;

    for(auto& kv : transf_){

        left  = std::get<0>(files_[kv.first]);
        right = std::get<1>(files_[kv.first]);

        // close streams
        if(!core_){
            if(direct_)
                fclose(file_status_[left].fp);
            fclose(file_status_[right].fp);
        }

        // delete stream holders
        if(!core_){
            if(direct_)
                file_status_.erase(left);
            file_status_.erase(right);
        }

        // delete file holders, sizes
        sizes_.erase(std::get<0>(files_[kv.first]));
        sizes_.erase(std::get<1>(files_[kv.first]));
        files_.erase(kv.first);
    }

    // tranposes too
    tsizes_.clear();

    // what transformations?
    if(core_)
        transf_core_.clear();
    transf_.clear();
}
std::pair<size_t, size_t> DF_Helper::identify_order()
{
    // Identify order of transformations to use strategic intermediates
    std::sort(sorted_spaces_.begin(), sorted_spaces_.end(),
      [](const std::pair<std::string, size_t> &left,
         const std::pair<std::string, size_t> &right){
            return left.second < right.second;});

    // copy transf_ keys into a list of needs
    std::list<std::string> needs;
    for(auto const& itr : transf_)
        needs.push_back(itr.first);

    // construct best transformation order
    size_t largest=0, maximum=0, small, large, op;
    for(size_t i=0; i<sorted_spaces_.size(); i++){
        bool on = false;
        size_t st = 0;
        std::string str = sorted_spaces_[i].first;
        std::list<std::string>::iterator itr, end;
        for(itr=needs.begin(), end=needs.end(); itr != end; ++itr){
            op = 0;
            op = (!(std::get<0>(transf_[*itr]).compare(str)) ? 1 : op);
            op = (!(std::get<1>(transf_[*itr]).compare(str)) ? 2 : op);
            if(op!=0)
            {
                if(!on){
                    bspace_.push_back(str);
                    on=true;
                }
                small = (op==1 ? std::get<1>(spaces_[std::get<0>(transf_[*itr])]) :
                    std::get<1>(spaces_[std::get<1>(transf_[*itr])]));
                large = (op==1 ? std::get<1>(spaces_[std::get<1>(transf_[*itr])]) :
                    std::get<1>(spaces_[std::get<0>(transf_[*itr])]));
                maximum = (maximum<small*large ? small*large : maximum);
                largest = (largest < small ? small : largest);
                order_.push_back(*itr);
                st++;
                needs.erase(itr);
                itr--;
            }
        }
        if(st>0) {strides_.push_back(st);}
    }
    print_order();
    ordered_ = true;
return std::make_pair(largest, maximum);
}
void DF_Helper::print_order()
{
    size_t o = order_.size();
    size_t b = bspace_.size();
    outfile->Printf("\n     ==> DF_Helper:--Begin Transformations Information <==\n\n");
    outfile->Printf("   Transformation order:\n");
    for(size_t i=0; i<o; i++){
        outfile->Printf("         %s: (%s, %s)\n", order_[i].c_str(),
          std::get<0>(transf_[order_[i]]).c_str(),std::get<1>(transf_[order_[i]]).c_str());
    }
    outfile->Printf("\n    Best Spaces:\n");
    for(size_t i=0; i<b; i++){
        outfile->Printf("         (space: %s, size: %zu)\n", bspace_[i].c_str(), std::get<1>(spaces_[bspace_[i]]));
    }
    outfile->Printf("\n    Transformation strides: ");
    for(size_t i=0; i<b; i++){
        outfile->Printf("%zu", strides_[i]);
        if(i<b-1)
            outfile->Printf(", ");
    }
    outfile->Printf("\n\n     ==> DF_Helper:--End Transformations Information <==\n\n");
}
void DF_Helper::transform()
{
    timer_on("DF_Helper~transform         ");
    core_ ? transform_core() : transform_disk();
    timer_off("DF_Helper~transform         ");
}
void DF_Helper::transform_core()
{
    if(!ordered_)
        identify_order();

    outfile->Printf("\n     ==> DF_Helper:--Begin Transformations (core)<==\n\n");

    // reset tranposes
    tsizes_.clear();

    // determine largest buffers needed
    size_t wMO  = std::get<1>(sorted_spaces_[sorted_spaces_.size()-1]);

    // enhance cache use
    size_t naux = naux_;
    size_t nao = nao_;
    int rank = 0;
    std::vector<std::vector<double>> C_buffers(nthreads_);
    #pragma omp parallel private(rank) num_threads(nthreads_)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(nao*wMO);
        C_buffers[rank] = Cp;
    }

    // stripe transformed integrals
    for(auto& kv : transf_){
        size_t size = std::get<1>(spaces_[std::get<0>(kv.second)])*std::get<1>(spaces_[std::get<1>(kv.second)]);
        transf_core_[kv.first].reserve(size*naux);
    }

    // declare bufs
    std::vector<double> T; T.reserve(naux_*nao*wMO);
    std::vector<double> F; F.reserve(naux_*wMO*wMO);
    std::vector<double> N; N.reserve(naux_*wMO*wMO);
    double *Mp = &Ppq_[0];
    double *Tp = &T[0];
    double *Fp = &F[0];
    double *Np = &N[0];

    size_t start = 0;
    size_t stop  = Qshells_-1;
    size_t block_size = naux;
    size_t begin = 0;
    size_t end = naux-1;

    size_t count=0;
    for(size_t i=0; i<bspace_.size(); i++){

        // grab best space
        std::string bspace = bspace_[i];
        size_t bsize = std::get<1>(spaces_[bspace]);
        double* Bpt = std::get<0>(spaces_[bspace])->pointer()[0];

        timer_on("DF_Helper~transform - first ");
        // form temp, thread over spM (nao)
        #pragma omp parallel for firstprivate(nao, naux, bsize, block_size) private(rank) schedule(guided) num_threads(nthreads_)
        for(size_t k=0; k<nao_; k++){

            // truncate transformation matrix according to fun_mask
            size_t jump = (big_skips_[k]*block_size)/naux;
            size_t sp_size = small_skips_[k];

#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            for(size_t m=0, sp_count= -1; m<nao_; m++){
                if(schwarz_fun_mask_[k*nao_+m]){
                    sp_count++;
                    C_DCOPY(bsize, &Bpt[m*bsize], 1, &C_buffers[rank][sp_count*bsize], 1);
                }
            }

            // (Qm)(mb)->(Qb)
            C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0,
              &Mp[jump], sp_size, &C_buffers[rank][0], bsize, 0.0, &Tp[k*block_size*bsize], bsize);
        }
        timer_off("DF_Helper~transform - first ");

        // to completion per transformation
        for(size_t k=0; k<strides_[i]; k++){
            std::string left = std::get<0>(transf_[order_[count + k]]);
            std::string right = std::get<1>(transf_[order_[count + k]]);
            std::tuple<SharedMatrix, size_t> I =
                (bspace.compare(left) == 0 ? spaces_[right] : spaces_[left]);

            size_t wsize = std::get<1>(I);
            double* Wp = std::get<0>(I)->pointer()[0];

            // (wp)(p|Qb)->(w|Qb)
            timer_on("DF_Helper~transform - second");
            C_DGEMM('T', 'N', wsize, block_size*bsize, nao_, 1.0, Wp, wsize,
              Tp, block_size*bsize, 0.0, Fp, block_size*bsize);
            timer_off("DF_Helper~transform - second");

            // Transpose in memory for convenient formats
            size_t st1 = bsize*wsize;
            size_t st2 = bsize*block_size;

            if(bspace.compare(left)==0){ // (w|Qb)->(Q|bw)
                timer_on("DF_Helper~transform - transp");
                #pragma omp parallel for num_threads(nthreads_) firstprivate(st1, st2)
                for(size_t x=0; x<block_size; x++){
                    for(size_t y=0; y<bsize; y++){
                        for(size_t z=0; z<wsize; z++){
                            Np[x*st1+y*wsize+z]=Fp[z*st2+x*bsize+y];
                        }
                    }
                }
                timer_off("DF_Helper~transform - transp");
            }
            else{ // (w|Qb)->(Q|wb)
                timer_on("DF_Helper~transform - transp");
                #pragma omp parallel for num_threads(nthreads_) firstprivate(st1, st2)
                for(size_t x=0; x<block_size; x++){
                    for(size_t z=0; z<wsize; z++){
                        #pragma omp simd
                        for(size_t y=0; y<bsize; y++){
                            Np[x * st1 + z * bsize + y] = Fp[z * st2 + x * bsize + y];
                        }
                    }
                }
                timer_off("DF_Helper~transform - transp");
            }
            double* Lp = transf_core_[order_[count+k]].data();
            if(direct_){
                double* metp = metric_prep_core(mpower_);
                timer_on("DF_Helper~direct contraction");
                C_DGEMM('N', 'N', naux, st1, naux, 1.0, metp, naux, Np, st1, 0.0, Lp, st1);
                timer_off("DF_Helper~direct contraction");
            }
            else{
                timer_on("DF_Helper~transform - write ");
                C_DCOPY(naux*st1, Np, 1, Lp, 1);
                timer_off("DF_Helper~transform - write ");
            }
        }
        count += strides_[i];
    }
    outfile->Printf("\n     ==> DF_Helper:--End Transformations (core)<==\n\n");
}

void DF_Helper::transform_disk() {

    if(!ordered_)
        info_ = identify_order();

    outfile->Printf("\n     ==> DF_Helper:--Begin Transformations (disk)<==\n\n");

    // reset transposes
    tsizes_.clear();

    // determine largest buffers needed
    size_t tots = std::get<1>(Qlargest_);
    size_t wMO = std::get<0>(info_);
    size_t max = std::get<1>(info_);
    wMO_ = wMO;

    // prep stream, blocking
    if (!direct_)
        stream_check(AO_files_[AO_names_[1]], "rb");
    else{
        Qsteps_.clear();
        Qlargest_ = shell_blocks(memory_, max, 0, Qshells_, Qsteps_);
    }

    timer_on("DF_Helper~transform - setup ");
    // enhance cache use
    size_t naux = naux_;
    size_t nao = nao_;
    int rank = 0;
    std::vector<std::vector<double>> C_buffers(nthreads_);
    // prepare eri buffers
    size_t nthread = nthreads_;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthread);
#pragma omp parallel private(rank) num_threads(nthreads_)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(nao * wMO);
        C_buffers[rank] = Cp;
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // scope these declarations
    {
        // declare bufs
        std::vector<double> M;
        M.reserve(std::get<0>(Qlargest_));
        std::vector<double> T;
        T.reserve(naux_ * nao * wMO);
        std::vector<double> F;
        F.reserve(naux_ * max);
        std::vector<double> N;
        N.reserve(naux_ * max);
        double* Mp = M.data();
        double* Tp = T.data();
        double* Fp = F.data();
        double* Np = N.data();
        timer_off("DF_Helper~transform - setup ");

        // transform in steps
        for (size_t j = 0; j < Qsteps_.size(); j++) {
            // Qshell step info
            size_t start = std::get<0>(Qsteps_[j]);
            size_t stop = std::get<1>(Qsteps_[j]);
            size_t begin = Qshell_aggs_[start];
            size_t end = Qshell_aggs_[stop + 1] - 1;
            size_t block_size = end - begin + 1;

            // print step info
            outfile->Printf("      Qshell: (%zu, %zu)", start, stop);
            outfile->Printf(", PHI: (%zu, %zu), size: %zu\n", begin, end, block_size);

            // get AO chunk according to directive
            timer_on("DF_Helper~transform - grab  ");
            (direct_ ? compute_AO_Q(start, stop, Mp, eri) : grab_AO(start, stop, Mp));
            timer_off("DF_Helper~transform - grab  ");

            size_t count = 0;
            for (size_t i = 0; i < bspace_.size(); i++) {
                // grab best space
                std::string bspace = bspace_[i];
                size_t bsize = std::get<1>(spaces_[bspace]);
                double* Bpt = std::get<0>(spaces_[bspace])->pointer()[0];

                timer_on("DF_Helper~transform - first ");
// form temp, thread over spM (nao)
#pragma omp parallel for firstprivate(nao, naux, bsize, block_size) private(rank) schedule(guided) \
                                          num_threads(nthreads_)
                for (size_t k = 0; k < nao_; k++) {
                    // truncate transformation matrix according to fun_mask
                    size_t jump = (big_skips_[k] * block_size) / naux_;
                    size_t sp_size = small_skips_[k];

#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    for (size_t m = 0, sp_count = -1; m < nao_; m++) {
                        if (schwarz_fun_mask_[k * nao_ + m]) {
                            sp_count++;
                            C_DCOPY(bsize, &Bpt[m * bsize], 1, &C_buffers[rank][sp_count * bsize],
                                    1);
                        }
                    }

                    // (Qm)(mb)->(Qb)
                    C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0, &Mp[jump], sp_size,
                            &C_buffers[rank][0], bsize, 0.0, &Tp[k * block_size * bsize], bsize);
                }
                timer_off("DF_Helper~transform - first ");

                // to completion per transformation
                for (size_t k = 0; k < strides_[i]; k++) {
                    std::string left = std::get<0>(transf_[order_[count + k]]);
                    std::string right = std::get<1>(transf_[order_[count + k]]);
                    std::tuple<SharedMatrix, size_t> I =
                        (bspace.compare(left) == 0 ? spaces_[right] : spaces_[left]);

                    size_t wsize = std::get<1>(I);
                    double* Wp = std::get<0>(I)->pointer()[0];

                    timer_on("DF_Helper~transform - second");
                    // (wp)(p|Qb)->(w|Qb)
                    C_DGEMM('T', 'N', wsize, block_size * bsize, nao_, 1.0, Wp, wsize, Tp,
                            block_size * bsize, 0.0, Fp, block_size * bsize);
                    timer_off("DF_Helper~transform - second");

                    // setup putf
                    std::string putf = (!direct_ ? std::get<1>(files_[order_[count + k]])
                                                 : std::get<0>(files_[order_[count + k]]));

                    std::string op = ( j ? "ab" : "wb");

                    // Transpose in memory for convenient formats
                    size_t st1 = bsize * wsize;
                    size_t st2 = bsize * block_size;

                    if (bspace.compare(left) == 0) { // (w|Qb)->(Q|bw)
                        timer_on("DF_Helper~transform - transp");
                        #pragma omp parallel for num_threads(nthreads_) firstprivate(st1, st2)
                        for (size_t x = 0; x < block_size; x++) {
                            for (size_t y = 0; y < bsize; y++) {
                                for (size_t z = 0; z < wsize; z++) {
                                    Np[x * st1 + y * wsize + z] = Fp[z * st2 + x * bsize + y];
                                }
                            }
                        }
                        timer_off("DF_Helper~transform - transp");
                        timer_on("DF_Helper~transform - write ");
                        put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, bsize - 1),
                                             std::make_pair(0, wsize - 1), op);
                        timer_off("DF_Helper~transform - write ");
                    } else {  // (w|Qb)->(Q|wb)
                        timer_on("DF_Helper~transform - transp");
                        #pragma omp parallel for num_threads(nthreads_) firstprivate(st1, st2)
                        for (size_t x = 0; x < block_size; x++) {
                            for (size_t z = 0; z < wsize; z++) {
                                #pragma omp simd
                                for (size_t y = 0; y < bsize; y++) {
                                    Np[x * st1 + z * bsize + y] = Fp[z * st2 + x * bsize + y];
                                }
                            }
                        }
                        timer_off("DF_Helper~transform - transp");
                        timer_on("DF_Helper~transform - write ");
                        put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, wsize - 1),
                                             std::make_pair(0, bsize - 1), op);
                        timer_off("DF_Helper~transform - write ");
                    }
                }
                count += strides_[i];
            }
        }
        outfile->Printf("\n     ==> DF_Helper:--End Transformations (disk)<==\n\n");
    }
    if (direct_) {
        // total size allowed, in doubles
        size_t total_mem = (memory_ - naux_ * naux_) / 2;
        timer_on("DF_Helper~transform - setup2");
        std::vector<double> M;
        M.reserve(total_mem);
        std::vector<double> F;
        F.reserve(total_mem);
        double* Mp = M.data();
        double* Fp = F.data();
        double* metp;
        timer_off("DF_Helper~transform - setup2");
        std::vector<std::string>::iterator itr;

        if (hold_met_) {
            timer_on("DF_Helper~transform - getmet");
            std::shared_ptr<FittingMetric> Jinv(new FittingMetric(aux_, true));
            Jinv->form_fitting_metric();
            SharedMatrix metric = Jinv->get_metric();
            metric->power(mpower_, 1e-12);
            metp = metric->pointer()[0];
            timer_off("DF_Helper~transform - getmet");
            timer_on("DF_Helper~direct contracts  ");
            for (itr = order_.begin(); itr != order_.end(); itr++)
                contract_metric(*itr, metp, Mp, Fp, total_mem);
        } else {
            std::vector<double> metric;
            metric.reserve(naux * naux);
            metp = metric.data();
            timer_on("DF_Helper~transform - getmet");
            std::string mfilename = return_metfile(mpower_);
            get_tensor_(std::get<0>(files_[mfilename]), metp, 0, naux_ - 1, 0, naux_ - 1);
            timer_off("DF_Helper~transform - getmet");
            timer_on("DF_Helper~direct contracts  ");
            for (itr = order_.begin(); itr != order_.end(); itr++)
                contract_metric(*itr, metp, Mp, Fp, total_mem);
        }
        timer_off("DF_Helper~direct contracts  ");
    }
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M) {
    fill_tensor(name, M, std::make_pair(0,0), std::make_pair(0,0), std::make_pair(0,0));
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> a1) {
    fill_tensor(name, M, a1, std::make_pair(0, 0), std::make_pair(0, 0));
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> a1,
                           std::pair<size_t, size_t> a2) {
    fill_tensor(name, M, a1, a2, std::make_pair(0, 0));
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0,
                           std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2) {

    std::string filename = std::get<1>(files_[name]);
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (std::get<1>(t0) == 0)
        t0 = std::make_pair(0, std::get<0>(sizes) - 1);
    if (std::get<1>(t1) == 0)
        t1 = std::make_pair(0, std::get<1>(sizes) - 1);
    if (std::get<1>(t2) == 0)
        t2 = std::make_pair(0, std::get<2>(sizes) - 1);

    check_transformation_name(name);
    check_transformation_tuple(name, t0, t1, t2);
    check_transformation_matrix(name, M, t0, t1, t2);

    size_t sta0 = std::get<0>(t0);
    size_t sto0 = std::get<1>(t0);
    size_t sta1 = std::get<0>(t1);
    size_t sto1 = std::get<1>(t1);
    size_t sta2 = std::get<0>(t2);
    size_t sto2 = std::get<1>(t2);

    size_t A0 = (sto0-sta0+1);
    size_t A1 = (sto1-sta1+1);
    size_t A2 = (sto2-sta2+1);

    double* Mp = M->pointer()[0];
    if (core_) {

        size_t a0 = std::get<0>(sizes);
        size_t a1 = std::get<1>(sizes);
        size_t a2 = std::get<2>(sizes);

        double* Fp = &transf_core_[name][0];
        #pragma omp parallel num_threads(nthreads_)
        for (size_t i = 0; i < A0; i++) {
            for (size_t j = 0; j < A1; j++) {
                #pragma omp simd
                for(size_t k = 0; k < A2; k++){
                    Mp[i*A1*A2 + j*A2 + k] = Fp[(sta0+i)*a1*a2 + (sta1+j)*a2 + (sta2+k)];
                }
            }
        }
    } else { get_tensor_(filename, Mp, t0, t1, t2); }
    M->set_numpy_shape({(int)A0, (int)A1, (int)A2});
}
SharedMatrix DF_Helper::get_tensor(std::string name) {
 return get_tensor(name, std::make_pair(0,0), std::make_pair(0,0), std::make_pair(0,0));
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::pair<size_t, size_t> a1) {

 return get_tensor(name, a1, std::make_pair(0, 0), std::make_pair(0, 0));
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::pair<size_t, size_t> a1,
                                   std::pair<size_t, size_t> a2) {
 return get_tensor(name, a1, a2, std::make_pair(0, 0));
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::pair<size_t, size_t> t0,
                                   std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2) {

    std::string filename = std::get<1>(files_[name]);
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (std::get<1>(t0) == 0)
        t0 = std::make_pair(0, std::get<0>(sizes) - 1);
    if (std::get<1>(t1) == 0)
        t1 = std::make_pair(0, std::get<1>(sizes) - 1);
    if (std::get<1>(t2) == 0)
        t2 = std::make_pair(0, std::get<2>(sizes) - 1);

    check_transformation_name(name);
    check_transformation_tuple(name, t0, t1, t2);

    size_t sta0 = std::get<0>(t0);
    size_t sto0 = std::get<1>(t0);
    size_t sta1 = std::get<0>(t1);
    size_t sto1 = std::get<1>(t1);
    size_t sta2 = std::get<0>(t2);
    size_t sto2 = std::get<1>(t2);

    size_t A0 = (sto0-sta0+1);
    size_t A1 = (sto1-sta1+1);
    size_t A2 = (sto2-sta2+1);

    SharedMatrix M(new Matrix("M", A0, A1*A2));
    double* Mp = M->pointer()[0];

    if (core_) {

        size_t a0 = std::get<0>(sizes);
        size_t a1 = std::get<1>(sizes);
        size_t a2 = std::get<2>(sizes);

        double* Fp = transf_core_[name].data();
        #pragma omp parallel num_threads(nthreads_)
        for (size_t i = 0; i < A0; i++) {
            for (size_t j = 0; j < A1; j++) {
                #pragma omp simd
                for(size_t k = 0; k < A2; k++){
                    Mp[i*A1*A2 + j*A2 + k] = Fp[(sta0+i)*a1*a2 + (sta1+j)*a2 + (sta2+k)];
                }
            }
        }
    } else { get_tensor_(filename, Mp, t0, t1, t2); }
    M->set_numpy_shape({(int)A0, (int)A1, (int)A2});
    return M;
}
void DF_Helper::check_transformation_name(std::string name){
    if (files_.find(name) == files_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::check_transformation_matrix(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0,
    std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2)
{
    size_t A0 = std::get<1>(t0) - std::get<0>(t1) + 1;
    size_t A1 =(std::get<1>(t1) - std::get<0>(t1) + 1) * (std::get<1>(t2) - std::get<0>(t2) + 1);

    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if (A0 != a0 && A1 != a1) {
        std::stringstream error;
        error << "DF_Helper:get_tensor: your matrix contridicts your tuple sizes when obtaining the" <<
            name << "integral.  ";
        error << "you gave me a matrix of size: (" << a0 << "," << a1 << "), but tuple sizes give:(" <<
            A0 << "," << A1 << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::check_transformation_tuple(std::string name, std::pair<size_t, size_t> t0,
    std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2)
{
    size_t sta0 = std::get<0>(t0);
    size_t sto0 = std::get<1>(t0);
    size_t sta1 = std::get<0>(t1);
    size_t sto1 = std::get<1>(t1);
    size_t sta2 = std::get<0>(t2);
    size_t sto2 = std::get<1>(t2);
    std::string filename = std::get<1>(files_[name]);

    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (sta0 > sto0) {
        std::stringstream error;
        error << "your axis 0 tuple has a larger start index: " << sta0 <<
            " than its stop index: " << sto0;
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (sta1 > sto1) {
        std::stringstream error;
        error << "your axis 1 tuple has a larger start index: " << sta1 <<
            " than its stop index: " << sto1;
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (sta2 > sto2) {
        std::stringstream error;
        error << "your axis 2 tuple has a larger start index: " << sta2 <<
            " than its stop index: " << sto2;
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M0 = std::get<0>(sizes);
    if (sto0 > M0 - 1) {
        std::stringstream error;
        error << "your axis 0 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered ("<< sto0 << "), but bounds is ("<< M0-1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M1 = std::get<1>(sizes);
    if (sto1 > M1 - 1) {
        std::stringstream error;
        error << "your axis 1 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered ("<< sto1 << "), but bounds is ("<< M1-1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M2 = std::get<2>(sizes);
    if (sto2 > M2 - 1) {
        std::stringstream error;
        error << "your axis 2 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered ("<< sto2 << "), but bounds is ("<< M2-1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::transpose(std::string name, std::tuple<size_t, size_t, size_t> order){
    (core_ ? transpose_core(name, order) : transpose_disk(name, order));
}
void DF_Helper::transpose_core(std::string name, std::tuple<size_t, size_t, size_t> order) {
    std::vector<double> M = transf_core_[name];
    double* Mp = M.data();
    double* Fp = transf_core_[name].data();

    size_t a0 = std::get<0>(order);
    size_t a1 = std::get<1>(order);
    size_t a2 = std::get<2>(order);

    std::string filename = std::get<1>(files_[name]);
    size_t M0 = std::get<0>(sizes_[filename]);
    size_t M1 = std::get<1>(sizes_[filename]);
    size_t M2 = std::get<2>(sizes_[filename]);
    std::tuple<size_t, size_t, size_t> sizes;

    bool on = false;
    if (a0 == 0) {
        if (a1 == 2) {
            sizes = std::make_tuple(M0, M2, M1);
            on = true;
        }
    } else if (a0 == 1) {
        if (a1 == 0) {
            sizes = std::make_tuple(M1, M0, M2);
            on = true;
        } else if (a1 == 2) {
            sizes = std::make_tuple(M1, M2, M0);
            on = true;
        }
    } else {
        if (a1 == 0) {
            sizes = std::make_tuple(M2, M0, M1);
            on = true;
        } else if (a1 == 1) {
            sizes = std::make_tuple(M2, M1, M0);
            on = true;
        }
    }

    if(!on)
        throw PSIEXCEPTION("you transposed all wrong!");

    if (a0 == 0) {
        if (a1 == 2) {  // (0|12) -> (0|21)
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[i * M1 * M2 + k * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    } else if (a0 == 1) {
        if (a1 == 0) {  // (0|12) -> (1|02)
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
#pragma omp simd
                    for (size_t k = 0; k < M2; k++) {
                        Fp[j * M0 * M2 + i * M2 + k] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 2) {  // (0|12) -> (1|20)
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[j * M0 * M2 + k * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    } else if (a0 == 2) {
        if (a1 == 0) {  // (0|12) -> (2|01)
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[k * M1 * M0 + i * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 1) {  // (0|12) -> (2|10)
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[k * M1 * M0 + j * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    }
    tsizes_[filename] = sizes;
}
void DF_Helper::transpose_disk(std::string name, std::tuple<size_t, size_t, size_t> order){

    size_t a0 = std::get<0>(order);
    size_t a1 = std::get<1>(order);
    size_t a2 = std::get<2>(order);

    // determine blocking
    std::string filename = std::get<1>(files_[name]);
    size_t M0 = std::get<0>(sizes_[filename]);
    size_t M1 = std::get<1>(sizes_[filename]);
    size_t M2 = std::get<2>(sizes_[filename]);

    size_t current=0, count=0, largest=0;;
    std::vector<std::pair<size_t, size_t>> steps;
    for(size_t i=0; i<M0; i++){
        current += M1*M2;
        count++;
        if(current*2 > memory_ +1000000|| i==M0-1){
            if(count==1 && i!=M0-1){
                std::stringstream error;
                error << "DF_Helper:contract_metric: not enough memory.";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if(i == M0-1)
                steps.push_back(std::make_pair(i-count+1, i));
            else{
                current -= M1*M2;
                steps.push_back(std::make_pair(i-count+1, i-1));
                i--;
            }
            if(largest < current)
                largest = current;
            count=0;
            current=0;
        }
    }

    // declare
    std::vector<double> M;
    std::vector<double> F;
    M.reserve(largest);
    F.reserve(largest);
    double* Mp = M.data();
    double* Fp = F.data();
    std::tuple<size_t, size_t, size_t> sizes;

    bool on = false;
    if (a0 == 0) {
        if (a1 == 2) {
            sizes = std::make_tuple(M0, M2, M1);
            on = true;
        }
    } else if (a0 == 1) {
        if (a1 == 0) {
            sizes = std::make_tuple(M1, M0, M2);
            on = true;
        } else if (a1 == 2) {
            sizes = std::make_tuple(M1, M2, M0);
            on = true;
        }
    } else {
        if (a1 == 0) {
            sizes = std::make_tuple(M2, M0, M1);
            on = true;
        } else if (a1 == 1) {
            sizes = std::make_tuple(M2, M1, M0);
            on = true;
        }
    }

    if(!on)
        throw PSIEXCEPTION("you transposed all wrong!");

    std::string new_file = "newfilefortransposition";
    filename_maker(new_file, std::get<0>(sizes), std::get<1>(sizes), std::get<2>(sizes));
    std::string new_filename = std::get<1>(files_[new_file]);
    std::string op = "wb";

    for(size_t m=0; m<steps.size(); m++){

        size_t start = std::get<0>(steps[m]);
        size_t stop  = std::get<1>(steps[m]);
        M0 = stop-start+1;

        // grab
        get_tensor_(filename, Mp, start, stop, 0, M1*M2 - 1);

        if (a0 == 0) {
            if (a1 == 2) {  // (0|12) -> (0|21)
#pragma omp parallel num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[i * M1 * M2 + k * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(start, stop), std::make_pair(0, M2 - 1),
                           std::make_pair(0, M1 - 1), op);
            }
        } else if (a0 == 1) {
            if (a1 == 0) {  // (0|12) -> (1|02)
#pragma omp parallel num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
#pragma omp simd
                        for (size_t k = 0; k < M2; k++) {
                            Fp[j * M0 * M2 + i * M2 + k] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M1 - 1), std::make_pair(start, stop),
                           std::make_pair(0, M2 - 1), op);
            } else if (a1 == 2) {  // (0|12) -> (1|20)
#pragma omp parallel num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[j * M0 * M2 + k * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M1 - 1), std::make_pair(0, M2 - 1),
                           std::make_pair(start, stop), op);
            }
        } else if (a0 == 2) {
            if (a1 == 0) {  // (0|12) -> (2|01)
#pragma omp parallel num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[k * M1 * M0 + i * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M2 - 1), std::make_pair(start, stop),
                           std::make_pair(0, M1 - 1), op);
            } else if (a1 == 1) {  // (0|12) -> (2|10)
#pragma omp parallel num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[k * M1 * M0 + j * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M2 - 1), std::make_pair(0, M1 - 1),
                           std::make_pair(start, stop), op);
            }
        }
    }
    // better be careful
    remove(filename.c_str());
    rename(new_filename.c_str(), filename.c_str());

    file_status_[filename] = file_status_[new_filename];
    stream_check(filename, "rb");
    file_status_.erase(new_filename);

    files_.erase(new_file);
    tsizes_[filename] = sizes;
}
size_t DF_Helper::get_space_size(std::string name){
    if (spaces_.find(name) == spaces_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_space_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return std::get<1>(spaces_[name]);
}
size_t DF_Helper::get_tensor_size(std::string name){
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    std::tuple<size_t, size_t, size_t> s = sizes_[std::get<1>(files_[name])];
    return std::get<0>(s)*std::get<1>(s)*std::get<2>(s);
}
std::tuple<size_t, size_t, size_t> DF_Helper::get_tensor_shape(std::string name){
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return sizes_[std::get<1>(files_[name])];
}
}
}  // End namespaces
