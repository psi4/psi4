
#include "schwarz.h"
#include "filler.h"
#include "aobasis.h"

using namespace yeti;
using namespace std;

// For now, all basis sets must be the same
CauchySchwarzValueEstimater::CauchySchwarzValueEstimater(
    TEIShellComputeFunctorPtr tbint,
    AOBasisPtr aobasis,
    TensorIndexDescr* descr
) : descr_(descr->get(0))
{
    range_ = descr_->get_range(1);
    size_ = range_->nelements();
    MultiShellMapPtr shmap = aobasis->get_multishell_map();
    //TODO Permutational symmetry
    A = new double*[size_];
    for(uli idx = 0; idx < size_; ++idx) {
        A[idx] = new double[size_];
    }

    uli start = range_->start();
    uli stop = range_->stop();
    for(uli i = start; i < stop; ++i) {
        uli nshi = shmap->nshell(i);
        uli istart = shmap->shell_start(i);
        for(uli j = start; j < stop; ++j) {
            double max_element = -1.0;
            uli nshj = shmap->nshell(j);
            uli jstart = shmap->shell_start(j);
            for(uli ish = istart; ish < istart + nshi; ++ish) {
                uli nfxn_shi = shmap->shell_size(ish);
                for(uli jsh = jstart; jsh < jstart + nshj; ++jsh) {
                    uli nfxn_shj = shmap->shell_size(jsh);
                    (*tbint)(
                        aobasis->program_shell_number(ish),
                        aobasis->program_shell_number(jsh),
                        aobasis->program_shell_number(ish),
                        aobasis->program_shell_number(jsh)
                    );
                    const double* intptr = tbint->buffer();
                    for(usi idx1 = 0; idx1 < nfxn_shi; ++idx1) {
                        for(usi idx2 = 0; idx2 < nfxn_shj; ++idx2) {
                            for(usi idx3 = 0; idx3 < nfxn_shi; ++idx3) {
                                if(idx3 != idx1) {
                                    intptr += nfxn_shj;
                                }
                                else {
                                    for(usi idx4 = 0; idx4 < nfxn_shj; ++idx4, ++intptr) {
                                        if(idx4 != idx2) {
                                            continue;
                                        }
                                        else if(*intptr > max_element) {
                                            max_element = 1.0/2.0 * fabs(*intptr);
                                        }
                                    }
                                }
                            }
                        }
                    } // end loops over functions in shell

                }
            } // end loop over shells in i, j

            A[i-start][j-start] = (max_element > 0.0 ? log10(max_element) : LOG_ZERO);
        }
    } // end loop over indices in range

}



CauchySchwarzValueEstimater::CauchySchwarzValueEstimater(
    const CauchySchwarzValueEstimater* sub_est,
    TensorIndexDescr* descr,
    usi depth
) : descr_(descr->get(0))
{
    range_ = descr_->get_range(depth);
    size_ = range_->nelements();
    //TODO Permutational symmetry
    A = new double*[size_];
    for(uli idx = 0; idx < size_; ++idx) {
        A[idx] = new double[size_];
    }

    for(uli i = 0; i < range_->nelements(); ++i) {
        IndexRange* subrangei = range_->get_subindex(i);
        for(uli j = 0; j < range_->nelements(); ++j) {
            IndexRange* subrangej = range_->get_subindex(j);
            double max_log = LOG_NONZERO;

            for(uli isub = subrangei->start(); isub < subrangei->stop(); ++isub) {
                for(uli jsub = subrangej->start(); jsub < subrangej->stop(); ++jsub) {
                    double subestval = sub_est->A[isub][jsub];
                    if(subestval > max_log) {
                        max_log = subestval;
                    }
                }
            }

        }
    }

}


CauchySchwarzValueEstimater::~CauchySchwarzValueEstimater()
{
    for(int i = 0; i < size_; ++i) {
        delete[] A[i];
    }
    delete[] A;
}

float
CauchySchwarzValueEstimater::max_log(
    const uli *indices
) const
{
    return A[indices[0]][indices[1]] + A[indices[2]][indices[3]];
}

