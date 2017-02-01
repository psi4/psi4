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

#include "psi4/psi4-dec.h"

#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"
#include "sobasis.h"
#include "matrix.h"
#include "molecule.h"

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP

namespace psi {

OneBodySOInt::OneBodySOInt(const std::shared_ptr<OneBodyAOInt> & ob,
                           const std::shared_ptr<IntegralFactory>& integral)
    : ob_(ob), integral_(integral.get()), deriv_(ob->deriv())
{
    common_init();

}

OneBodySOInt::OneBodySOInt(const std::shared_ptr<OneBodyAOInt> & ob,
                           const IntegralFactory* integral)
    : ob_(ob), integral_(integral), deriv_(ob->deriv())
{
    common_init();
}

OneBodySOInt::~OneBodySOInt()
{
}

void OneBodySOInt::common_init()
{
    b1_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis1(), integral_));

    if (ob_->basis2() == ob_->basis1())
        b2_ = b1_;
    else
        b2_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis2(), integral_));

    ob_->set_force_cartesian(b1_->petite_list()->include_pure_transform());
}

std::shared_ptr<SOBasisSet> OneBodySOInt::basis() const
{
    return b1_;
}

std::shared_ptr<SOBasisSet> OneBodySOInt::basis1() const
{
    return b1_;
}

std::shared_ptr<SOBasisSet> OneBodySOInt::basis2() const
{
    return b2_;
}

std::shared_ptr<OneBodyAOInt> OneBodySOInt::ob() const
{
    return ob_;
}

void OneBodySOInt::compute(SharedMatrix result)
{
    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            const SOTransform &t1 = b1_->sotrans(ish);
            const SOTransform &t2 = b2_->sotrans(jsh);

            int nao2 = b2_->naofunction(jsh);

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];
                    ob_->compute_shell(s1.aoshell, s2.aoshell);

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;

                            // Check the irreps to ensure symmetric quantities.
                            if (ifunc.irrep == jfunc.irrep)
                                result->add(ifunc.irrep,
                                            b1_->function_within_irrep(ish, isofunc),
                                            b2_->function_within_irrep(jsh, jsofunc),
                                            jcoef * aobuf[jaooff]);
                        }
                    }
                }
            }
        }
    }
}

void OneBodySOInt::compute(std::vector<SharedMatrix > results)
{
    // Do not worry about zeroing out result
    int nchunk = ob_->nchunk();
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            const SOTransform &t1 = b1_->sotrans(ish);
            const SOTransform &t2 = b2_->sotrans(jsh);

            int nao1 = b1_->naofunction(ish);
            int nao2 = b2_->naofunction(jsh);
            int nao = nao1*nao2;

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    ob_->compute_shell(s1.aoshell, s2.aoshell);

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;

                            // Handle chunks
                            for (int i=0; i<nchunk; ++i) {
                                double temp = jcoef * aobuf[jaooff + (i*nao)];

                                int ijirrep = ifunc.irrep ^ jfunc.irrep;
                                if (ijirrep == results[i]->symmetry()) {
                                    // Add the contribution to the matrix
                                    results[i]->add(ifunc.irrep,
                                                    b1_->function_within_irrep(ish, isofunc),
                                                    b2_->function_within_irrep(jsh, jsofunc),
                                                    temp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void OneBodySOInt::compute_deriv1(std::vector<SharedMatrix > result,
                                  const CdSalcList &cdsalcs)
{
    // Do not worry about zeroing out result.

    // Do some checks:
    if (deriv_ < 1)
        throw SanityCheckError("OneBodySOInt::compute_deriv1: integral object not created to handle derivatives.", __FILE__, __LINE__);

    if (result.size() != cdsalcs.ncd())
        throw SanityCheckError("OneBodySOInt::compute_deriv1: result vector size does not match SALC size.", __FILE__, __LINE__);

    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over unique SO shells.
    for (int ish=0; ish<ns1; ++ish) {
        const SOTransform& t1 = b1_->sotrans(ish);
        int nao1 = b1_->naofunction(ish);

        for (int jsh=0; jsh<ns2; ++jsh) {
            const SOTransform& t2= b2_->sotrans(jsh);
            int nao2 = b2_->naofunction(jsh);

            int nao12 = nao1 * nao2;

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                int center_i = ob_->basis1()->shell(s1.aoshell).ncenter();
                const CdSalcWRTAtom& cdsalc1 = cdsalcs.atom_salc(center_i);

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];
                    int center_j = ob_->basis2()->shell(s2.aoshell).ncenter();
                    const CdSalcWRTAtom& cdsalc2 = cdsalcs.atom_salc(center_j);

                    // If we're working on the same atomic center, don't even bother with the derivative
                    if (center_i == center_j)
                        continue;

                    ob_->compute_shell_deriv1(s1.aoshell, s2.aoshell);

                    // handle SO transform
                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        // SO transform coefficient
                        double icoef = ifunc.coef;
                        // AO function offset in a linear array
                        int iaofunc  = ifunc.aofunc;
                        // SO function offset in a linear array
                        int isofunc  = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        // AO function offset in a linear array
                        int iaooff   = iaofunc;
                        // Relative position of the SO function within its irrep
                        int irel     = b1_->function_within_irrep(ish, isofunc);
                        int iirrep   = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc  = jfunc.aofunc;
                            int jsofunc  = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff   = iaooff*nao2 + jaofunc;
                            int jrel     = b2_->function_within_irrep(jsh, jsofunc);
                            int jirrep   = jfunc.irrep;

                            // Need to loop over the cdsalcs

                            double jcoef_aobuf = jcoef * aobuf[jaooff + 0*nao12];
                            for (int nx=0; nx<cdsalc1.nx(); ++nx) {
                                const CdSalcWRTAtom::Component element = cdsalc1.x(nx);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }

                            jcoef_aobuf = jcoef * aobuf[jaooff + 1*nao12];
                            for (int ny=0; ny<cdsalc1.ny(); ++ny) {
                                const CdSalcWRTAtom::Component element = cdsalc1.y(ny);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }

                            jcoef_aobuf = jcoef * aobuf[jaooff + 2*nao12];
                            for (int nz=0; nz<cdsalc1.nz(); ++nz) {
                                const CdSalcWRTAtom::Component element = cdsalc1.z(nz);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }

                            jcoef_aobuf = jcoef * aobuf[jaooff + 3*nao12];
                            for (int nx=0; nx<cdsalc2.nx(); ++nx) {
                                const CdSalcWRTAtom::Component element = cdsalc2.x(nx);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }

                            jcoef_aobuf = jcoef * aobuf[jaooff + 4*nao12];
                            for (int ny=0; ny<cdsalc2.ny(); ++ny) {
                                const CdSalcWRTAtom::Component element = cdsalc2.y(ny);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }

                            jcoef_aobuf = jcoef * aobuf[jaooff + 5*nao12];
                            for (int nz=0; nz<cdsalc2.nz(); ++nz) {
                                const CdSalcWRTAtom::Component element = cdsalc2.z(nz);
                                double temp = jcoef_aobuf * element.coef;
                                if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                    result[element.salc]->add(iirrep, irel, jrel, temp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

TwoBodySOInt::TwoBodySOInt(const std::shared_ptr<TwoBodyAOInt> &tb,
                           const std::shared_ptr<IntegralFactory>& integral) :

    integral_(integral), only_totally_symmetric_(false), cdsalcs_(0)
{
    tb_.push_back(tb);
    common_init();

}

TwoBodySOInt::TwoBodySOInt(const std::vector<std::shared_ptr<TwoBodyAOInt> > &tb,
                           const std::shared_ptr<IntegralFactory>& integral) :

    tb_(tb), integral_(integral), only_totally_symmetric_(false), cdsalcs_(0)
{
    common_init();

}

TwoBodySOInt::TwoBodySOInt(const std::shared_ptr<TwoBodyAOInt>& tb,
                           const std::shared_ptr<IntegralFactory>& integral,
                           const CdSalcList& cdsalcs) :

    integral_(integral), only_totally_symmetric_(false), cdsalcs_(&cdsalcs)
{
    tb_.push_back(tb);
    common_init();

}

TwoBodySOInt::TwoBodySOInt(const std::vector<std::shared_ptr<TwoBodyAOInt> >& tb,
                           const std::shared_ptr<IntegralFactory>& integral,
                           const CdSalcList& cdsalcs) :

    tb_(tb), integral_(integral), only_totally_symmetric_(false), cdsalcs_(&cdsalcs)
{
    common_init();


}

void TwoBodySOInt::common_init()
{
    // MPI runtime settings (defaults provided for all communicators)
    nthread_ = Process::environment.get_n_threads();
    nproc_   = 1;
    me_      = 0;

    // Try to reduce some work:
    b1_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(tb_[0]->basis1(), integral_));

    if (tb_[0]->basis1() == tb_[0]->basis2())
        b2_ = b1_;
    else
        b2_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(tb_[0]->basis2(), integral_));

    if (tb_[0]->basis1() == tb_[0]->basis3())
        b3_ = b1_;
    else
        b3_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(tb_[0]->basis3(), integral_));

    if (tb_[0]->basis3() == tb_[0]->basis4())
        b4_ = b3_;
    else
        b4_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(tb_[0]->basis4(), integral_));

    for (int i=0; i<nthread_; ++i)
        tb_[i]->set_force_cartesian(b1_->petite_list()->include_pure_transform());

    size_ = b1_->max_nfunction_in_shell() *
            b2_->max_nfunction_in_shell() *
            b3_->max_nfunction_in_shell() *
            b4_->max_nfunction_in_shell();

//    outfile->Printf( "aQRS %zu, abcD %zu, abRS %zu, max_size %zu, size %zu\n", aQRS, abcD, abRS, max_size, size_);

    // Check to make sure things are consistent
    if (tb_[0]->deriv() > 0 && cdsalcs_ == 0)
        throw PSIEXCEPTION("TwoBodySOInt: AO object configured for gradients but CdSalc object not provided.");
    if (tb_[0]->deriv() == 0 && cdsalcs_ != 0)
        throw PSIEXCEPTION("TwoBodySOInt: Ao object not configured for gradients but CdSalc was provided.");

    if (tb_[0]->deriv() == 1 && cdsalcs_ != 0)
        for (int i=0; i<nthread_; ++i)
            buffer_.push_back(new double[size_ * cdsalcs_->ncd()]);
    else {
        for (int i=0; i<nthread_; ++i) {
            buffer_.push_back(new double[size_]);
//            temp_.push_back(new double[max_size]);
//            temp2_.push_back(new double[abRS]);
        }
    }

    ::memset(iirrepoff_, 0, sizeof(int) * 8);
    ::memset(jirrepoff_, 0, sizeof(int) * 8);
    ::memset(kirrepoff_, 0, sizeof(int) * 8);
    ::memset(lirrepoff_, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        iirrepoff_[h] = iirrepoff_[h-1] + b1_->nfunction_in_irrep(h-1);
        jirrepoff_[h] = jirrepoff_[h-1] + b2_->nfunction_in_irrep(h-1);
        kirrepoff_[h] = kirrepoff_[h-1] + b3_->nfunction_in_irrep(h-1);
        lirrepoff_[h] = lirrepoff_[h-1] + b4_->nfunction_in_irrep(h-1);
    }

    // We need the petite list.
    petite1_ = b1_->petite_list();
    petite2_ = b2_->petite_list();
    petite3_ = b3_->petite_list();
    petite4_ = b4_->petite_list();

    pg_ = tb_[0]->basis()->molecule()->point_group();

    dcd_ = std::shared_ptr<DCD>(new DCD(petite1_->group()));

    if (cdsalcs_) {
        int ncd = cdsalcs_->ncd();
        for (int n=0; n<nthread_; ++n) {
            deriv_.push_back(new double*[ncd]);
            for (int i=0; i<ncd; ++i)
                deriv_[n][i] = &(buffer_[n][i*size_]);
        }
    }

    cutoff_ = Process::environment.options.get_double("INTS_TOLERANCE");
}

TwoBodySOInt::~TwoBodySOInt()
{
    for (int i=0; i<nthread_; ++i) {
        delete[] buffer_[i];
//        if (tb_[0]->deriv() == 0) { // 'if' statement to be removed
//            delete[] temp_[i];
//            delete[] temp2_[i];
//        }
        if (deriv_.size())
            delete[] deriv_[i];
    }
}

std::shared_ptr<SOBasisSet> TwoBodySOInt::basis() const
{
    return b1_;
}

std::shared_ptr<SOBasisSet> TwoBodySOInt::basis1() const
{
    return b1_;
}

std::shared_ptr<SOBasisSet> TwoBodySOInt::basis2() const
{
    return b2_;
}

std::shared_ptr<SOBasisSet> TwoBodySOInt::basis3() const
{
    return b3_;
}

std::shared_ptr<SOBasisSet> TwoBodySOInt::basis4() const
{
    return b4_;
}

}
