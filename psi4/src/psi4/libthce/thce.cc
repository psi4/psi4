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


#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "thce.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"

#include <unistd.h>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi {

THCE::THCE()
{
}
THCE::~THCE()
{
}
void THCE::print(std::string out, int level) const
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   if (level >= 0) {
        printer->Printf("  ==> THCE <==\n\n");

        printer->Printf("  Tensors    = %11zu [--]\n", tensors_.size());
        printer->Printf("  Dimensions = %11zu [--]\n", dimensions_.size());
        printer->Printf("  Core       = %11zu [MB]\n", (core_doubles() * 8L) / (1024L * 1024L));
        printer->Printf("  Disk       = %11zu [MB]\n", (disk_doubles() * 8L) / (1024L * 1024L));
        printer->Printf("\n");

        printer->Printf("  Dimensions:\n\n");
        printer->Printf("  %11s %11s\n", "Name", "Size");
        for(std::map<std::string, int>::const_iterator it = dimensions_.begin();
            it != dimensions_.end(); ++it) {
            printer->Printf("  %11s %11d\n", (*it).first.c_str(), (*it).second);
        }
        printer->Printf("\n");

        printer->Printf("  Tensors:\n\n");
        printer->Printf("  %11s %11s %11s %11s %11s\n", "Alias", "Name", "Order", "Storage", "Trust");
        for(std::map<std::string, std::shared_ptr<Tensor> >::const_iterator it = tensors_.begin();
            it != tensors_.end(); ++it) {
            std::string key = (*it).first;
            std::shared_ptr<Tensor> T = (*it).second;
            printer->Printf("  %11s %11s %11d %11s %11s\n", key.c_str(), T->name().c_str(), T->order(),
                (T->disk() ? "Disk" : "Core"), (T->trust() ? "Yes" : "No"));
        }
        printer->Printf("\n");
    }
    if (level >= 1) {
        printer->Printf("  Tensor Details:\n\n");
        for(std::map<std::string, std::shared_ptr<Tensor> >::const_iterator it = tensors_.begin();
            it != tensors_.end(); ++it) {
            (*it).second->print(out,level);
        }
    }
    //if (level >= 0) {
    //    printer->Printf("  \"More matter with less art\"\n");
    //    printer->Printf("\n");
    //}

}
size_t THCE::core_doubles() const
{
    // Counts only unique (non-aliased tensors)
    std::set<std::string> names;

    size_t val = 0L;
    for(std::map<std::string, std::shared_ptr<Tensor> >::const_iterator it = tensors_.begin();
        it != tensors_.end(); ++it) {
        std::shared_ptr<Tensor> T = (*it).second;
        if (!names.count(T->name())) {
            val += T->core_doubles();
            names.insert(T->name());
        }
    }
    return val;
}
size_t THCE::disk_doubles() const
{
    // Counts only unique (non-aliased tensors)
    std::set<std::string> names;

    size_t val = 0L;
    for(std::map<std::string, std::shared_ptr<Tensor> >::const_iterator it = tensors_.begin();
        it != tensors_.end(); ++it) {
        std::shared_ptr<Tensor> T = (*it).second;
        if (!names.count(T->name())) {
            val += T->disk_doubles();
            names.insert(T->name());
        }
    }
    return val;
}
void THCE::dimension_check(const std::string& name)
{
    if (!dimensions_.count(name)) {
        throw PSIEXCEPTION("Dimension " + name + " has not been declared.");
    }
}
void THCE::tensor_check(const std::string& name)
{
    if (!tensors_.count(name)) {
        throw PSIEXCEPTION("Tensor " + name + " has not been declared.");
    }
}
void THCE::new_dimension(const std::string& name, int size)
{
    dimensions_[name] = size;
}
void THCE::delete_dimension(const std::string& name)
{
    dimensions_.erase(name);
}
void THCE::new_core_tensor(const std::string& name, const std::string& dimensions, double* data, bool trust)
{
    std::vector<std::string> dims;
    std::vector<int> sizes;

    if (dimensions.length() != 0) {
        dims = split(dimensions, ",");
        for (int i = 0; i < dims.size(); i++) {
            dimension_check(dims[i]);
            sizes.push_back(dimensions_[dims[i]]);
        }
    }

    std::shared_ptr<Tensor> T(new CoreTensor(name,dims,sizes,data,trust));

    tensors_[name] = T;
}
void THCE::new_disk_tensor(const std::string& name, const std::string& dimensions, bool save, bool load)
{
    std::vector<std::string> dims;
    std::vector<int> sizes;

    if (dimensions.length() != 0) {
        dims = split(dimensions, ",");
        for (int i = 0; i < dims.size(); i++) {
            dimension_check(dims[i]);
            sizes.push_back(dimensions_[dims[i]]);
        }
    }

    std::shared_ptr<Tensor> T(new DiskTensor(name,dims,sizes,save,load));

    tensors_[name] = T;
}
void THCE::add_tensor(const std::string& name, std::shared_ptr<Tensor> tensor)
{
    tensors_[name] = tensor;
}
void THCE::delete_tensor(const std::string& name)
{
    tensors_.erase(name);
}

long int Tensor::unique_id = 0;

Tensor::Tensor(const std::string& name,
        std::vector<string>& dimensions,
        std::vector<int>& sizes) :
    name_(name), dimensions_(dimensions), sizes_(sizes)
{
    if (sizes_.size() != dimensions_.size()) {
        throw PSIEXCEPTION("Dimensions and Sizes are not the same order.");
    }
    order_ = sizes_.size();
    active_sizes_ = sizes_;
    numel_ = 1L;
    for (int k = 0; k < order_; k++) {
        numel_ *= sizes_[k];
    }

    set_filename();
}
Tensor::~Tensor()
{
}
void Tensor::set_filename()
{
    std::stringstream ss;
    ss <<  PSIOManager::shared_object()->get_default_path();
    ss <<  "/";
    ss << psi_file_prefix;
    ss << ".";
    ss << getpid();
    ss << ".";
    ss << PSIO::get_default_namespace();
    ss << ".";
    ss << Tensor::unique_id;
    ss << ".";
    ss << name_;
    ss << ".dat";
    filename_ = ss.str();

    Tensor::unique_id++;
}
void Tensor::slice(std::shared_ptr<Tensor> A, std::vector<std::tuple<bool,int,int,bool,int,int> >& topology)
{
    // => Topology Metadata/Validity <= //

    // Pointer to C
    double* Cp;
    // Pointer to A
    double* Ap;

    // File pointer to C
    FILE* FCp;
    // File pointer to A
    FILE* FAp;

    // Memory topology type (0: C = C, 1: D = C, 2: C = D)
    int mem_type = 0;

    if (disk() && A->disk()) {
        throw PSIEXCEPTION("Disk to disk slice is not permitted.");
    } else if (A->disk()) {
        swap_check();
        FAp = A->file_pointer();
        Cp = pointer();
        mem_type = 2;
    } else if (disk()) {
        A->swap_check();
        FCp = file_pointer();
        Ap = A->pointer();
        mem_type = 1;
    } else {
        A->swap_check();
        swap_check();
        Cp = pointer();
        Ap = A->pointer();
        mem_type = 0;
    }

    // Sizes in C
    std::vector<int>& sizeC = sizes_;
    // Sizes in A
    std::vector<int>& sizeA = A->sizes();

    // Active Sizes in C
    std::vector<int>& actC = active_sizes_;
    // Active Sizes in A
    std::vector<int>& actA = A->active_sizes();

    // Indices in C
    std::vector<std::string>& indC = dimensions_;
    // Indices in A
    std::vector<std::string>& indA = A->dimensions();

    // Number of indices (including united singletons)
    int nindex = topology.size();

    // Start indices in C
    std::vector<int> startC(nindex);
    // Start indices in A
    std::vector<int> startA(nindex);

    // End indices in C
    std::vector<int> endC(nindex);
    // End indices in A
    std::vector<int> endA(nindex);

    // End - Start in C
    std::vector<int> deltaC(nindex);
    // End - Start in A
    std::vector<int> deltaA(nindex);

    // Strides in C
    std::vector<size_t> strideC(nindex);
    // Strides in A
    std::vector<size_t> strideA(nindex);

    // Superblock opportunities
    std::vector<bool> full(nindex);

    // TODO: Error checking

    size_t size_data = 1L;
    size_t LDA = 1L;
    size_t LDC = 1L;
    int rA = A->order() - 1;
    int rC = order_ - 1;
    for (int ind = nindex-1; ind >= 0; ind--) {
        bool isC = std::get<0>(topology[ind]);
        bool isA = std::get<3>(topology[ind]);
        int sC = std::get<1>(topology[ind]);
        int sA = std::get<4>(topology[ind]);
        int eC = std::get<2>(topology[ind]);
        int eA = std::get<5>(topology[ind]);

        if (!isC && !isA) {
            throw PSIEXCEPTION("What exactly is F/F?");
        } else if (!isC) {
            // Singleton C
            startC[ind] = 0;
            endC[ind] = 1;
            deltaC[ind] = 1;
            strideC[ind] = LDC;

            startA[ind] = sA;
            endA[ind] = eA;
            deltaA[ind] = 1;
            strideA[ind] = LDA;
            LDA *= sizeA[rA];

            if (sizeA[rA] == deltaA[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rA--;
        } else if (!isA) {
            // Singleton A
            startC[ind] = sC;
            endC[ind] = eC;
            deltaC[ind] = 1;
            strideC[ind] = LDC;
            LDC *= sizeC[rC];

            startA[ind] = 0;
            endA[ind] = 1;
            deltaA[ind] = 1;
            strideA[ind] = LDA;

            if (sizeC[rC] == deltaC[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rC--;
        } else if (sA == -1 && sC == -1) {
            // Full/Full
            actC[rC] = actA[rA];

            startC[ind] = 0;
            endC[ind] = actC[rC];
            deltaC[ind] = endC[ind] - startC[ind];
            strideC[ind] = LDC;
            LDC *= sizeC[rC];

            startA[ind] = 0;
            endA[ind] = actA[rA];
            deltaA[ind] = endA[ind] - startA[ind];
            strideA[ind] = LDA;
            LDA *= sizeA[rA];

            if (sizeA[rA] == deltaA[ind] && sizeC[rC] == deltaC[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rC--;
            rA--;
        } else if (sA == -1) {
            // Slice/Full

            startC[ind] = sC;
            endC[ind] = eC;
            deltaC[ind] = endC[ind] - startC[ind];
            strideC[ind] = LDC;
            LDC *= sizeC[rC];

            startA[ind] = 0;
            endA[ind] = eC - sC;
            deltaA[ind] = endA[ind] - startA[ind];
            strideA[ind] = LDA;
            LDA *= sizeA[rA];

            if (sizeA[rA] == deltaA[ind] && sizeC[rC] == deltaC[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rC--;
            rA--;
        } else if (sC == -1) {
            // Full/Slice
            actC[rC] = eA - sA;

            startC[ind] = 0;
            endC[ind] = actC[rC];
            deltaC[ind] = endC[ind] - startC[ind];
            strideC[ind] = LDC;
            LDC *= sizeC[rC];

            startA[ind] = sA;
            endA[ind] = eA;
            deltaA[ind] = endA[ind] - startA[ind];
            strideA[ind] = LDA;
            LDA *= sizeA[rA];

            if (sizeA[rA] == deltaA[ind] && sizeC[rC] == deltaC[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rC--;
            rA--;
        } else {
            // Slice/Slice
            if (eC - sC != eA -sA) {
                throw PSIEXCEPTION("Slice/Slice index ranges are not the same size");
            }

            startC[ind] = sC;
            endC[ind] = eC;
            deltaC[ind] = endC[ind] - startC[ind];
            strideC[ind] = LDC;
            LDC *= sizeC[rC];

            startA[ind] = sA;
            endA[ind] = eA;
            deltaA[ind] = endA[ind] - startA[ind];
            strideA[ind] = LDA;
            LDA *= sizeA[rA];

            if (sizeA[rA] == deltaA[ind] && sizeC[rC] == deltaC[ind]) {
                full[ind] = true;
            } else {
                full[ind] = false;
            }

            rC--;
            rA--;
        }
        size_data *= deltaA[ind];
    }

    // => Index Selection <= //

    int nfast = 0;
    size_t size_fast = 1L;
    for (int ind = nindex - 1; ind >= 0; ind--) {
        if (full[ind]) {
            size_fast *= deltaA[ind];
            nfast++;
        } else {
            break;
        }
    }
    int nslow = nindex - nfast;
    size_t size_slow = size_data / size_fast;

    // => Critical index absorption into superindex <= //

    size_t delta = 1L;
    if (nslow > 0) {
        delta = deltaA[nslow-1];
    }

    /**

    // => Debug printing <= //

    outfile->Printf( "  ==> Slice <==\n\n");
    outfile->Printf( "    Total Size: %11zu\n", size_data);
    outfile->Printf( "    Slow Size:  %11zu\n", size_slow / delta);
    outfile->Printf( "    Fast Size:  %11zu\n", size_fast * delta);
    outfile->Printf( "\n");

    outfile->Printf( "    Total Dims: %11zu\n", nindex);
    outfile->Printf( "    Slow Dims:  %11zu\n", (nslow > 0 ? nslow - 1 : nslow));
    outfile->Printf( "    Fast Dims:  %11zu\n", (nslow > 0 ? nfast + 1 : nfast));
    outfile->Printf( "\n");

    outfile->Printf( "    %11s %11s %11s %11s %11s %11s %11s %11s %11s\n",
        "startC", "endC", "deltaC", "strideC",
        "startA", "endA", "deltaA", "strideA",
        "full?");
    for (int ind = 0; ind < nindex; ind++) {
        outfile->Printf( "    %11d %11d %11d %11d %11d %11d %11d %11d %11s\n",
            startC[ind], endC[ind], deltaC[ind], strideC[ind],
            startA[ind], endA[ind], deltaA[ind], strideA[ind],
            (full[ind] ? "Yes" : "No"));
    }
    outfile->Printf( "\n");

    **/

    // => Copy Operations <= //

    // Offset in tensor A
    size_t offsetA;
    // Offset in tensor C
    size_t offsetC;
    // Paging values
    size_t num, den, val;

    // => Master Loop <= //

    for (size_t ind = 0L; ind < size_slow; ind+=delta) {

        // Compute the C offset
        offsetC = 0L;
        num = ind;
        den = size_slow;
        for (int dim = 0; dim < nslow; dim++) {
            den /= deltaC[dim];
            val = num / den + startC[dim]; // The current index of the dim'th dimension
            num -= val * den;
            offsetC += val * strideC[dim];
        }

        // Compute the A offset
        offsetA = 0L;
        num = ind;
        den = size_slow;
        for (int dim = 0; dim < nslow; dim++) {
            den /= deltaA[dim];
            val = num / den + startA[dim]; // The current index of the dim'th dimension
            num -= val * den;
            offsetA += val * strideA[dim];
        }

        // Do the copy
        if (mem_type == 0) {
            ::memcpy((void*) (Cp + offsetC), (void*) (Ap + offsetA), sizeof(double) * size_fast * delta);
        } else if (mem_type == 1) {
            fseek(FCp,offsetC*sizeof(double),SEEK_SET);
            fwrite((void*) (Ap + offsetA), sizeof(double), size_fast * delta, FCp);
        } else {
            fseek(FAp,offsetA*sizeof(double),SEEK_SET);
            size_t statusvalue=fread((void*) (Cp + offsetC), sizeof(double), size_fast * delta, FAp);
        }
    }

    if (mem_type == 1) {
        fflush(FCp);
    }
    if (mem_type == 2) {
        fflush(FAp);
    }
}

CoreTensor::CoreTensor(const std::string& name,
        std::vector<string>& dimensions, std::vector<int>& sizes,
        double* data,
        bool trust) :
        Tensor(name,dimensions,sizes),
        trust_(trust)
{
    if (trust_) {
        data_ = data;
    } else {
        data_ = new double[numel_];
        if (data  == NULL) {
            ::memset((void*) data_, '\0', sizeof(double) * numel_);
        } else {
            ::memcpy((void*) data_, (void*) data, sizeof(double) * numel_);
        }
    }

    swapped_ = false;
    fh_ = NULL;
}
CoreTensor::~CoreTensor()
{
    if (!trust_) {
        if (data_ != NULL) {
            delete[] data_;
            data_ = NULL;
        }
        if (fh_ != NULL) {
            fclose(fh_);
            fh_ = NULL;
            std::string file = filename();
            remove(file.c_str());
        }
    }
}
void CoreTensor::swap_check() const
{
    if (!core() | swapped())
        throw PSIEXCEPTION("Tensor is swapped out, cannot operate on it.");
}
std::shared_ptr<Tensor> CoreTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    const std::string& dimension3, int size3,
    const std::string& dimension4, int size4,
    double* data, bool trust)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);
    dimensions.push_back(dimension3);
    sizes.push_back(size3);
    dimensions.push_back(dimension4);
    sizes.push_back(size4);

    return std::shared_ptr<Tensor>(new CoreTensor(name,dimensions,sizes,data,trust));
}
std::shared_ptr<Tensor> CoreTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    const std::string& dimension3, int size3,
    double* data, bool trust)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);
    dimensions.push_back(dimension3);
    sizes.push_back(size3);

    return std::shared_ptr<Tensor>(new CoreTensor(name,dimensions,sizes,data,trust));
}
std::shared_ptr<Tensor> CoreTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    double* data, bool trust)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);

    return std::shared_ptr<Tensor>(new CoreTensor(name,dimensions,sizes,data,trust));
}
std::shared_ptr<Tensor> CoreTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    double* data, bool trust)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);

    return std::shared_ptr<Tensor>(new CoreTensor(name,dimensions,sizes,data,trust));
}
std::shared_ptr<Tensor> CoreTensor::build(const std::string& name,
    double* data, bool trust)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    return std::shared_ptr<Tensor>(new CoreTensor(name,dimensions,sizes,data,trust));
}
void CoreTensor::print(std::string out, int level) const
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   const int print_ncol = Process::environment.options.get_int("MAT_NUM_COLUMN_PRINT");
    if (level >= 0) {
        printer->Printf( "  => CoreTensor %s <=\n\n", name_.c_str());
        printer->Printf( "    Order   = %11d\n", order_);
        printer->Printf( "    Numel   = %11zu\n", numel_);
        printer->Printf( "    Swapped = %11s\n", swapped() ? "Yes" : "No");
        printer->Printf( "    Trust   = %11s\n", trust_ ? "Yes" : "No");
        printer->Printf( "\n");

        printer->Printf( "    Dimensions:\n\n");
        printer->Printf( "    %2s %11s %11s %11s\n", "N", "Name", "Alloc Size", "Active Size");
        for (int k = 0; k < order_; k++) {
            printer->Printf( "    %2d %11s %11d %11d\n", k+1, dimensions_[k].c_str(), sizes_[k], active_sizes_[k]);
        }
        printer->Printf( "\n");
    }
    if (level >= 2) {
        if (swapped()) {
            printer->Printf( "    CoreTensor is swapped out, data is unavailable to print.\n\n");
        } else {
            size_t page_size = 1L;
            int rows = 1;
            int cols = 1;
            if (order_ >= 1) {
                page_size *= sizes_[order_ - 1];
                rows = sizes_[order_ - 1];
            }
            if (order_ >= 2) {
                page_size *= sizes_[order_ - 2];
                rows = sizes_[order_ - 2];
                cols = sizes_[order_ - 1];
            }

            printer->Printf( "    Data:\n\n");

            size_t pages = numel_ / page_size;
            for (size_t page = 0L; page < pages; page++) {

                if (order_ > 2) {
                    printer->Printf( "    Page (");
                    size_t num = page;
                    size_t den = pages;
                    size_t val;
                    for (int k = 0; k < order_ - 2; k++) {
                        den /= sizes_[k];
                        val = num / den;
                        num -= val * den;
                        printer->Printf("%zu,",val);
                    }
                    printer->Printf( "*,*):\n\n");
                }

                double* vp = data_ + page * page_size;
                if (order_ == 0) {
                    printer->Printf( "    %12.7f\n", *(vp));
                    printer->Printf("\n");
                } else if(order_ == 1) {
                    for (int i=0; i<page_size; ++i) {
                        printer->Printf( "    %5d %12.7f\n", i, *(vp + i));
                    }
                    printer->Printf("\n");
                } else {
                    int nframes = cols / print_ncol;
                    for (int j = 0; j < cols; j+= print_ncol) {
                        int ncols = (j + print_ncol >= cols ? cols - j : print_ncol);

                        // Column Header
                        printer->Printf("    %5s", "");
                        for (int jj = j; jj < j+ncols; jj++) {
                            printer->Printf(" %12d", jj);
                        }
                        printer->Printf("\n");

                        // Data
                        for (int i = 0; i < rows; i++) {
                            printer->Printf("    %5d", i);
                            for (int jj = j; jj < j+ncols; jj++) {
                                printer->Printf(" %12.7f", *(vp + i * cols + jj));
                            }
                            printer->Printf("\n");
                        }

                        // Block separator
                        printer->Printf("\n");
                    }
                }
            }
        }
    }
}
void CoreTensor::set_pointer(double* data)
{
    if (!trust_) {
        throw PSIEXCEPTION("You cannot assign a trust pointer to a non-trust CoreTensor");
    }
    data_ = data;
}
void CoreTensor::set_data(double* data)
{
    swap_check();

    ::memcpy((void*) data_, (void*) data, sizeof(double) * numel_);
}
void CoreTensor::swap_out(bool changed)
{
    if (trust_) {
        throw PSIEXCEPTION("You can't swap a trust CoreTensor.");
    }

    // First swap

    if (fh_ == NULL) {
        std::string file = filename();
        fh_ = fopen(file.c_str(), "wb+");
        fwrite((void*) data_, sizeof(double), numel_, fh_);
        fseek(fh_,0L,SEEK_SET);
        delete[] data_;
        data_ = NULL;
        swapped_ = true;
        return;
    }

    // Subsequent swaps

    if (!swapped()) {
        if (changed) {
            fseek(fh_,0L,SEEK_SET);
            fwrite((void*) data_, sizeof(double), numel_, fh_);
            fseek(fh_,0L,SEEK_SET);
        }
        delete[] data_;
        data_ = NULL;
        swapped_ = true;
    }
}
void CoreTensor::swap_in(bool read)
{
    if (trust_) {
        throw PSIEXCEPTION("You can't swap a trust CoreTensor.");
    }

    if (swapped()) {
        data_ = new double[numel_];
        if (read) {
            fseek(fh_,0L,SEEK_SET);
            size_t statusvalue=fread((void*) data_, sizeof(double), numel_, fh_);
            fseek(fh_,0L,SEEK_SET);
        } else {
            ::memset(data_,'\0', numel_*sizeof(double));
        }
        swapped_ = false;
    }
}
void CoreTensor::zero()
{
    swap_check();

    ::memset((void*) data_, '\0',sizeof(double) *  numel_);
}
void CoreTensor::scale(double val)
{
    swap_check();

    C_DSCAL(numel_, val, data_, 1);
}
void CoreTensor::add(std::shared_ptr<Tensor> A, double alpha, double beta)
{
    swap_check();
    A->swap_check();

    scale(beta);

    if (numel_ != A->numel() || order_ != A->order())
        throw PSIEXCEPTION("Unlike tensors cannot be added");

    double* Ap = A->pointer();

    C_DAXPY(numel_,alpha,Ap,1,data_,1);
}
void CoreTensor::permute(std::shared_ptr<Tensor> A, std::vector<int>& orderC)
{
    // => Swap Check <= //

    swap_check();
    A->swap_check();

    // => Topology Metadata/Validity <= //

    // Pointer to C
    double* Cp = data_;
    // Pointer to A
    double* Ap = A->pointer();

    // Sizes in C
    std::vector<int>& sizeC = sizes_;
    // Sizes in A
    std::vector<int>& sizeA = A->sizes();

    // Active Sizes in C
    std::vector<int>& actC = active_sizes_;
    // Active Sizes in A
    std::vector<int>& actA = A->active_sizes();

    // Indices in C
    std::vector<std::string>& indC = dimensions_;
    // Indices in A
    std::vector<std::string>& indA = A->dimensions();

    // Total Number of Indices
    int nindex = orderC.size();

    if (nindex != order_) {
        throw PSIEXCEPTION("Permutation Topology Error: CoreTensor C has wrong order");
    }
    if (nindex != A->order()) {
        throw PSIEXCEPTION("Permutation Topology Error: CoreTensor A has wrong order");
    }

    // => Canonical Ordering (Order in A) <= //

    for (int rA = 0; rA < nindex; rA++) {
        int rC = orderC[rA];
        if (indC[rC] != indA[rA]) {
            throw PSIEXCEPTION("Permutation Topology Error: Non-matching (name) Index");
        }
        if (sizeC[rC] != sizeA[rA]) {
            throw PSIEXCEPTION("Permutation Topology Error: Non-matching (size) Index");
        }
        // Active Inheritance
        actC[rC] = actA[rA];
    }

    // => Block superindex extents <= //

    // Number of block indices
    int nfast = 1L;
    // Size of the block superindex
    size_t fast_size = sizeA[nindex-1];
    // Last index in the block superindex in C
    int fast_indC = orderC[nindex-1];
    // Where was the C-tensor index last pass?
    int last_ind = fast_indC;
    for (int ind = nindex-2; ind >= 0; ind--) {
        int ind2 = orderC[ind];
        if (ind2 + 1 == last_ind) {
            nfast++;
            fast_size *= sizeA[ind];
        } else {
            break;
        }
    }
    // Size of the fast superindex stride in C
    size_t fast_stride = 1L;
    for (int ind = fast_indC+1; ind < nindex; ind++) {
        fast_stride *= sizeC[ind];
    }

    // Number of slow indices
    int nslow = nindex - nfast;
    // Size of slow index
    size_t slow_size = numel_ / fast_size;

    // => C stride values <= //

    // Unordered stride values in C
    std::vector<size_t> strideCU(nindex);
    strideCU[nindex-1] = 1L;
    for (int ind = nindex-2; ind >= 0; ind--) {
        strideCU[ind] = strideCU[ind+1] * sizeC[ind+1];
    }
    // Ordered stride values in C
    std::vector<size_t> strideC(nindex);
    for (int ind = 0; ind < nindex; ind++) {
        strideC[orderC[ind]] = strideCU[ind];
    }

    // => Master Loop <= //

    // Current C pointer
    double* A2p = Ap;
    // Current A pointer
    double* C2p;
    // Paging values
    size_t num, den, val;
    for (size_t ind = 0L; ind < slow_size; ind++) {

        // Place the C pointer
        C2p = Cp;
        num = ind;
        den = slow_size;
        //size_t offset = 0L;
        //printf("Ind %4d: ", ind);
        for (int dim = 0; dim < nslow; dim++) {
            den /= sizeA[dim];
            val = num / den; // The current index of the dim'th dimension
            num -= val * den;
            C2p += val * strideC[dim];
            //offset += val * strideC[dim];
            //printf("%4d ", val);
        }
        //printf("%11zu\n", offset);

        // Do the copy
        C_DCOPY(fast_size,A2p,1,C2p,fast_stride);

        // Increment the A pointer
        A2p += fast_size;
    }
}
void CoreTensor::contract(std::shared_ptr<Tensor> A, std::shared_ptr<Tensor> B, std::vector<std::tuple<std::string,int,int,int> >& topology, double alpha, double beta)
{
    // => Swap Check <= //

    swap_check();
    A->swap_check();
    B->swap_check();

    // => Topology Metadata/Validity <= //

    // Pointer to C
    double* Cp = data_;
    // Pointer to A
    double* Ap = A->pointer();
    // Pointer to B
    double* Bp = B->pointer();

    // Indices in C
    std::vector<std::string>& indC = dimensions_;
    // Indices in A
    std::vector<std::string>& indA = A->dimensions();
    // Indices in B
    std::vector<std::string>& indB = B->dimensions();

    // Sizes in C
    std::vector<int>& sizeC = sizes_;
    // Sizes in A
    std::vector<int>& sizeA = A->sizes();
    // Sizes in B
    std::vector<int>& sizeB = B->sizes();

    // Active Sizes in C
    std::vector<int>& actC = active_sizes_;
    // Active Sizes in A
    std::vector<int>& actA = A->active_sizes();
    // Active Sizes in B
    std::vector<int>& actB = B->active_sizes();

    // Number of Indices
    int nindex = topology.size();
    // Number of Hadamard Indices
    int nH = 0;
    // Number of of Left Outer Indices
    int nL = 0;
    // Number of Right Outer Indices
    int nR = 0;
    // Number of Inner Indices
    int nI = 0;

    // Name if index
    std::vector<std::string> name(nindex);
    // 0 is Hadamard, 1 is outer left, 2 is outer right, 3 is inner
    std::vector<int> type(nindex);
    // Rank of the index in A (or -1 for not present, see type)
    std::vector<int> rankA(nindex);
    // Rank of the index in B (or -1 for not present, see type)
    std::vector<int> rankB(nindex);
    // Rank of the index in C (or -1 for not present, see type)
    std::vector<int> rankC(nindex);
    // Active size of the index
    std::vector<int> act(nindex);
    // True size of the index
    std::vector<int> size(nindex);
    // Is the index gimp?
    std::vector<bool> gimp(nindex);

    // Figure everything out
    for (int ind = 0; ind < nindex; ind++) {
        name[ind] = std::get<0>(topology[ind]);
        int rC = std::get<1>(topology[ind]);
        int rA = std::get<2>(topology[ind]);
        int rB = std::get<3>(topology[ind]);
        rankA[ind] = rA;
        rankB[ind] = rB;
        rankC[ind] = rC;
        if (rA >= 0 && rB >= 0 && rC >= 0) {
            // Type Spec
            type[ind] = 0;
            nH++;
            act[ind] = actA[rA];
            size[ind] = sizeA[rA];
            gimp[ind] = (act[ind] != size[ind]);
            // Active Inheritance
            actC[rC] = actA[rA];
            // Topology Errors
            if (indC[rC] != indA[rA] || indC[rC] != indB[rB]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (name) Index " + std::get<0>(topology[ind]));
            }
            if (sizeC[rC] != sizeA[rA] || sizeC[rC] != sizeB[rB]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (size) Index " + std::get<0>(topology[ind]));
            }
            if (actB[rB] != actA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (active) Index " + std::get<0>(topology[ind]));
            }
        } else if (rA >= 0 && rC >= 0) {
            // Type Spec
            type[ind] = 1;
            nL++;
            act[ind] = actA[rA];
            size[ind] = sizeA[rA];
            gimp[ind] = (act[ind] != size[ind]);
            // Active Inheritance
            actC[rC] = actA[rA];
            // Topology Errors
            if (indC[rC] != indA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (name) Index " + std::get<0>(topology[ind]));
            }
            if (sizeC[rC] != sizeA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (size) Index " + std::get<0>(topology[ind]));
            }
        } else if (rB >= 0 && rC >= 0) {
            // Type Spec
            type[ind] = 2;
            nR++;
            act[ind] = actB[rB];
            size[ind] = sizeB[rB];
            gimp[ind] = (act[ind] != size[ind]);
            // Active Inheritance
            actC[rC] = actB[rB];
            // Topology Errors
            if (indC[rC] != indB[rB]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (name) Index " + std::get<0>(topology[ind]));
            }
            if (sizeC[rC] != sizeB[rB]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (size) Index " + std::get<0>(topology[ind]));
            }
        } else if (rA >= 0 && rB >= 0) {
            // Type Spec
            type[ind] = 3;
            nI++;
            act[ind] = actA[rA];
            size[ind] = sizeA[rA];
            gimp[ind] = (act[ind] != size[ind]);
            // Topology Errors
            if (indB[rB] != indA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (name) Index " + std::get<0>(topology[ind]));
            }
            if (sizeB[rB] != sizeA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (size) Index " + std::get<0>(topology[ind]));
            }
            if (actB[rB] != actA[rA]) {
                throw PSIEXCEPTION("Contraction Topology Error: Non-matching (active) Index " + std::get<0>(topology[ind]));
            }
        } else {
            throw PSIEXCEPTION("Contraction Topology Error: Disconnected Index " + std::get<0>(topology[ind]));
        }
    }

    if (nH + nL + nR != order_) {
        throw PSIEXCEPTION("Contraction Topology Error: CoreTensor C has wrong order");
    }
    if (nH + nL + nI != A->order()) {
        throw PSIEXCEPTION("Contraction Topology Error: CoreTensor A has wrong order");
    }
    if (nH + nI + nR != B->order()) {
        throw PSIEXCEPTION("Contraction Topology Error: CoreTensor B has wrong order");
    }

    // => Classiness Check <= //

    // Type, rank, rank, rank (if Hadamard), original index
    std::vector<std::tuple<int,int,int,int,int> > canon;
    for (int ind = 0; ind < nindex; ind++) {
        if (type[ind] == 0) {
            canon.push_back(std::tuple<int,int,int,int,int>(0,rankC[ind],rankA[ind],rankB[ind],ind));
        } else if (type[ind] == 1) {
            canon.push_back(std::tuple<int,int,int,int,int>(1,rankC[ind],rankA[ind],rankB[ind],ind));
        } else if (type[ind] == 2) {
            canon.push_back(std::tuple<int,int,int,int,int>(2,rankC[ind],rankB[ind],rankA[ind],ind));
        } else {
            canon.push_back(std::tuple<int,int,int,int,int>(3,rankA[ind],rankB[ind],rankC[ind],ind));
        }
    }
    std::sort(canon.begin(),canon.end());

    // Hadamard check
    for (int ind = 0; ind < nH; ind++) {
        if (rankA[ind] != ind || rankB[ind] != ind || rankC[ind] != ind) {
            throw PSIEXCEPTION("Nonconforming Hadamard (must be slow in all three tensors)");
        }
    }

    // Is this operation conforming for BLAS optimization?
    bool classy = true;

    int offset = 0;
    for (int ind = 1; ind < nH; ind++) {
        int index = ind + offset;
        int orig1 = std::get<4>(canon[index-1]);
        int orig2 = std::get<4>(canon[index]);
        if (rankA[orig2] - rankA[orig1] != 1 || rankB[orig2] - rankB[orig1] != 1 || rankC[orig2] - rankC[orig1] != 1) {
            classy = false;
        }
        if (gimp[orig2]) {
            classy = false;
        }
    }
    offset += nH;
    for (int ind = 1; ind < nL; ind++) {
        int index = ind + offset;
        int orig1 = std::get<4>(canon[index-1]);
        int orig2 = std::get<4>(canon[index]);
        if (rankA[orig2] - rankA[orig1] != 1 || rankC[orig2] - rankC[orig1] != 1) {
            classy = false;
        }
        if (gimp[orig2]) {
            classy = false;
        }
    }
    offset += nL;
    for (int ind = 1; ind < nR; ind++) {
        int index = ind + offset;
        int orig1 = std::get<4>(canon[index-1]);
        int orig2 = std::get<4>(canon[index]);
        if (rankB[orig2] - rankB[orig1] != 1 || rankC[orig2] - rankC[orig1] != 1) {
            classy = false;
        }
        if (gimp[orig2]) {
            classy = false;
        }
    }
    offset += nR;
    for (int ind = 1; ind < nI; ind++) {
        int index = ind + offset;
        int orig1 = std::get<4>(canon[index-1]);
        int orig2 = std::get<4>(canon[index]);
        if (rankA[orig2] - rankA[orig1] != 1 || rankB[orig2] - rankB[orig1] != 1) {
            classy = false;
        }
        if (gimp[orig2]) {
            classy = false;
        }
    }


    if (classy) {

        // Full Hadamard size
        size_t sizeH = 1L;
        // Active Hadamard size (including possible slow gimp)
        size_t actH = 1L;
        // Full Left size
        size_t sizeL = 1L;
        // Active Left size (including possible slow gimp)
        size_t actL = 1L;
        // Full Right size
        size_t sizeR = 1L;
        // Active Right size (including possible slow gimp)
        size_t actR = 1L;
        // Full Inner size
        size_t sizeI = 1L;
        // Active Inner size (including possible slow gimp)
        size_t actI = 1L;

        // Size loops
        offset = 0;
        for (int ind = 0; ind < nH; ind++) {
            int index = ind + offset;
            int orig = std::get<4>(canon[index]);
            actH *= act[orig];
            sizeH *= size[orig];
        }
        offset += nH;
        for (int ind = 0; ind < nL; ind++) {
            int index = ind + offset;
            int orig = std::get<4>(canon[index]);
            actL *= act[orig];
            sizeL *= size[orig];
        }
        offset += nL;
        for (int ind = 0; ind < nR; ind++) {
            int index = ind + offset;
            int orig = std::get<4>(canon[index]);
            actR *= act[orig];
            sizeR *= size[orig];
        }
        offset += nR;
        for (int ind = 0; ind < nI; ind++) {
            int index = ind + offset;
            int orig = std::get<4>(canon[index]);
            actI *= act[orig];
            sizeI *= size[orig];
        }

        // L before R in C?
        bool LR = false;
        // L before I in A?
        bool LI = false;
        // R before I in B?
        bool RI = false;

        if (nL > 0 && nR > 0) {
            LR = (rankC[std::get<4>(canon[nH])] < rankC[std::get<4>(canon[nH+nL])]);
        }
        if (nL > 0 && nI > 0) {
            LI = (rankA[std::get<4>(canon[nH])] < rankA[std::get<4>(canon[nH+nL+nR])]);
        }
        if (nI > 0 && nR > 0) {
            RI = (rankB[std::get<4>(canon[nH+nL])] < rankB[std::get<4>(canon[nH+nL+nR])]);
        }

        // => Contraction Operations <= //

        // > True Hadamard (Special Case) < //
        if (nI == 0 && nL == 0 && nR == 0) {
            for (size_t indH = 0L; indH < actH; indH++) {
                (*Cp) = alpha * (*Ap) * (*Bp) + beta * (*Cp);
                Ap++;
                Bp++;
                Cp++;
            }
            return;
        }

        // > Generalized Hadamard Loop < //
        for (size_t indH = 0L; indH < actH; indH++) {

            if (nI != 0 && nL == 0 && nR == 0) {
                // > DOT < //
                (*Cp) = alpha * C_DDOT(actI,Ap,1,Bp,1) + beta * (*Cp);
            } else if (nI == 0 && nL != 0 && nR == 0) {
                // > AXPY < //
                C_DSCAL(actL,beta,Cp,1);
                C_DAXPY(actL,alpha*(*Bp),Ap,1,Cp,1);
            } else if (nI == 0 && nL == 0 && nR != 0) {
                // > AXPY < //
                C_DSCAL(actR,beta,Cp,1);
                C_DAXPY(actR,alpha*(*Ap),Bp,1,Cp,1);
            } else if (nI != 0 && nL != 0 && nR == 0) {
                // > GEMV < //
                C_DGEMV((LI ? 'N' : 'T'), (LI ? actL : actI),(LI ? actI : actL),alpha,Ap,(LI ? sizeI : sizeL),Bp,1,beta,Cp,1);
            } else if (nI != 0 && nL == 0 && nR != 0) {
                // > GEMV < //
                C_DGEMV((RI ? 'N' : 'T'), (RI ? actR : actI),(RI ? actI : actR),alpha,Bp,(RI ? sizeI : sizeR),Ap,1,beta,Cp,1);
            } else if (nI == 0 && nL != 0 && nR != 0) {
                // > GER < //
                C_DSCAL(actL*actR,beta,Cp,1);
                if (LR) {
                    C_DGER(actL,actR,alpha,Ap,1,Bp,1,Cp,sizeR);
                } else {
                    C_DGER(actR,actL,alpha,Bp,1,Ap,1,Cp,sizeL);
                }
            } else if (nI != 0 && nL != 0 && nR != 0) {
                // > GEMM < //
                if (LR) {
                    C_DGEMM((LI ? 'N' : 'T'),(RI ? 'T' : 'N'),actL,actR,actI,alpha,Ap,(LI ? sizeI : sizeL),Bp,(RI ? sizeI : sizeR),beta,Cp,sizeR);
                } else {
                    C_DGEMM((RI ? 'N' : 'T'),(LI ? 'T' : 'N'),actR,actL,actI,alpha,Bp,(RI ? sizeI : sizeR),Ap,(LI ? sizeI : sizeL),beta,Cp,sizeL);
                }
            } else {
                throw PSIEXCEPTION("Contraction Topology Error: How the fuck did this happen?");
            }

            // Pointer increment
            Cp += sizeL * sizeR;
            Ap += sizeL * sizeI;
            Bp += sizeR * sizeI;
        }

    } else {
        // Non-conforming case
        throw PSIEXCEPTION("Slow-ass code is not allowed. Do you like apples? How about them apples?");
    }
}

DiskTensor::DiskTensor(const std::string& name,
        std::vector<string>& dimensions, std::vector<int>& sizes,
        bool save, bool load) : Tensor(name,dimensions,sizes), save_(save)
{
    if (load) {
        fh_ = fopen(filename().c_str(),"rb+");
    } else {
        fh_ = fopen(filename().c_str(),"wb+");
    }
}
DiskTensor::~DiskTensor()
{
    fclose(fh_);
    if (!save_) {
        remove(filename().c_str());
    }
}
void DiskTensor::swap_check() const
{
    throw PSIEXCEPTION("Tensor is DiskTensor, cannot operate on it.");
}
std::shared_ptr<Tensor> DiskTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    const std::string& dimension3, int size3,
    const std::string& dimension4, int size4,
    bool save, bool load)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);
    dimensions.push_back(dimension3);
    sizes.push_back(size3);
    dimensions.push_back(dimension4);
    sizes.push_back(size4);

    return std::shared_ptr<Tensor>(new DiskTensor(name,dimensions,sizes,save,load));
}
std::shared_ptr<Tensor> DiskTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    const std::string& dimension3, int size3,
    bool save, bool load)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);
    dimensions.push_back(dimension3);
    sizes.push_back(size3);

    return std::shared_ptr<Tensor>(new DiskTensor(name,dimensions,sizes,save,load));
}
std::shared_ptr<Tensor> DiskTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    const std::string& dimension2, int size2,
    bool save, bool load)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);
    dimensions.push_back(dimension2);
    sizes.push_back(size2);

    return std::shared_ptr<Tensor>(new DiskTensor(name,dimensions,sizes,save,load));
}
std::shared_ptr<Tensor> DiskTensor::build(const std::string& name,
    const std::string& dimension1, int size1,
    bool save, bool load)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    dimensions.push_back(dimension1);
    sizes.push_back(size1);

    return std::shared_ptr<Tensor>(new DiskTensor(name,dimensions,sizes,save,load));
}
std::shared_ptr<Tensor> DiskTensor::build(const std::string& name,
    bool save, bool load)
{
    std::vector<std::string> dimensions;
    std::vector<int> sizes;

    return std::shared_ptr<Tensor>(new DiskTensor(name,dimensions,sizes,save,load));
}
void DiskTensor::print(std::string out, int level) const
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   if (level >= 0) {
        printer->Printf( "  => DiskTensor %s <=\n\n", name_.c_str());
        printer->Printf( "    File    = %s\n", filename().c_str());
        printer->Printf( "    Save    = %11s\n", (save_ ? "Yes" : "No"));
        printer->Printf( "    Order   = %11d\n", order_);
        printer->Printf( "    Numel   = %11zu\n", numel_);
        printer->Printf( "\n");

        printer->Printf( "    Dimensions:\n\n");
        printer->Printf( "    %2s %11s %11s %11s\n", "N", "Name", "Alloc Size", "Active Size");
        for (int k = 0; k < order_; k++) {
            printer->Printf( "    %2d %11s %11d %11d\n", k+1, dimensions_[k].c_str(), sizes_[k], active_sizes_[k]);
        }
        printer->Printf( "\n");
    }

}
void DiskTensor::zero()
{
    // Executive decision: stripe up to two indices at once
    size_t fast_size = 1L;
    for (int ind = order_ - 1; ind >= 0 && ind > order_ - 3; ind--) {
        fast_size *= sizes_[ind];
    }
    size_t slow_size = numel_ / fast_size;

    double* buf = new double[fast_size];
    ::memset((void*) buf, '\0', sizeof(double) * fast_size);

    fseek(fh_,0L,SEEK_END);
    for (size_t ind = 0L; ind < slow_size; ind++) {
        fwrite((void*) buf, sizeof(double), fast_size, fh_);
    }
    fflush(fh_);

    delete[] buf;
}

}
