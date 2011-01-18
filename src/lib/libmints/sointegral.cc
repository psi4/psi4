#include "sointegral.h"
#include "twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"
#include "sobasis.h"
#include "matrix.h"

#include <boost/shared_ptr.hpp>

#define DEBUG

namespace psi {

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyInt> & ob,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : ob_(ob), integral_(integral.get())
{
    b1_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis1(), integral));

    //    b1_->print();

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis2(), integral));

    only_totally_symmetric_ = 0;

    buffer_ = new double[INT_NCART(ob->basis1()->max_am())
            *INT_NCART(ob->basis2()->max_am())];
}

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyInt> & ob,
                           const IntegralFactory* integral)
    : ob_(ob), integral_(integral)
{
    b1_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis1(), integral));

    //    b1_->print();

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis2(), integral));

    only_totally_symmetric_ = 0;

    buffer_ = new double[INT_NCART(ob->basis1()->max_am())
            *INT_NCART(ob->basis2()->max_am())];
}

OneBodySOInt::~OneBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis2() const
{
    return b2_;
}

void OneBodySOInt::compute_shell(int ish, int jsh)
{
    const double *aobuf = ob_->buffer();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);

    memset(buffer_, 0, nso1*nso2*sizeof(double));

    int nao2 = b2_->naofunction(jsh);

    // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
    // and then using the transformation coefficients to obtain everything else.
    // Otherwise using the petite list doesn't save us any computational time
    // in computing the integrals, but does save us time when we use the integrals.

    // loop through the AO shells that make up this SO shell
    for (int i=0; i<t1.naoshell; ++i) {
        const SOTransformShell &s1 = t1.aoshell[i];
        for (int j=0; j<t2.naoshell; ++j) {
            const SOTransformShell &s2 = t2.aoshell[j];

            //            fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
            ob_->compute_shell(s1.aoshell, s2.aoshell);

            for (int itr=0; itr<s1.nfunc; ++itr) {
                const SOTransformFunction &ifunc = s1.func[itr];
                double icoef = ifunc.coef;
                int iaofunc = ifunc.aofunc;
                int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                int iaooff = iaofunc;
                int isooff = isofunc;

                for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                    const SOTransformFunction &jfunc = s2.func[jtr];
                    double jcoef = jfunc.coef * icoef;
                    int jaofunc = jfunc.aofunc;
                    int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                    int jaooff = iaooff*nao2 + jaofunc;
                    int jsooff = isooff*nso2 + jsofunc;

                    buffer_[jsooff] += jcoef * aobuf[jaooff];

#ifdef DEBUG
                    //                    if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
                    //                        fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                    //                                isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                    //                                ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                    //                                jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
                    //                    }
#endif
                }
            }
        }
    }
}

void OneBodySOInt::compute(boost::shared_ptr<Matrix> result)
{
    // Check symmetry.
//    if (b1_->nirrep() == 1) {
//        // Molecule has C1 symmetry. Need to use the ob_ object to get a SimpleMatrix.
//        throw FeatureNotImplemented("libmints", "OneBodySOInt with C1 symmetry. Need to direct the code to OneBodyInt.", __FILE__, __LINE__);
//    }

    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            //            fprintf(outfile, "computing ish = %d jsh = %d\n", ish, jsh);

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);

            memset(buffer_, 0, nso1*nso2*sizeof(double));

            int nao2 = b2_->naofunction(jsh);

            // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
            // and then using the transformation coefficients to obtain everything else.
            // Otherwise using the petite list doesn't save us any computational time
            // in computing the integrals, but does save us time when we use the integrals.

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    //                    fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
                    ob_->compute_shell(s1.aoshell, s2.aoshell);
                    //                    for (int z=0; z < INT_NPURE(ob_->basis1()->shell(s1.aoshell)->am()) *
                    //                         INT_NPURE(ob_->basis2()->shell(s2.aoshell)->am()); ++z) {
                    //                        fprintf(outfile, "raw: %d -> %8.5f\n", z, aobuf[z]);
                    //                    }

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        int iirrep = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            int jirrep = jfunc.irrep;

                            buffer_[jsooff] += jcoef * aobuf[jaooff];

                            //                            if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
                            //                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                            //                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                            //                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
                            //                                fprintf(outfile, "(%d|%d) += %8.5f * (%d|%d): %8.5f -> %8.5f\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
                            //                            }

                            // Check the irreps to ensure symmetric quantities.
                            if (ifunc.irrep == jfunc.irrep)
                                result->add(ifunc.irrep, b1_->function_within_irrep(ish, isofunc), b2_->function_within_irrep(jsh, jsofunc), jcoef * aobuf[jaooff]);
                        }
                    }
                }
            }
        }
    }
}

void OneBodySOInt::compute_petitelist(boost::shared_ptr<Matrix> result)
{
    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Check symmetry.
    if (b1_->nirrep() == 1) {
        // Molecule has C1 symmetry. Need to use the ob_ object to get a SimpleMatrix.
        throw FeatureNotImplemented("libmints", "OneBodySOInt with C1 symmetry. Need to direct the code to OneBodyInt.", __FILE__, __LINE__);
    }

    shared_ptr<PetiteList> petite = b1_->petitelist();

    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<=ish; ++jsh) {

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);

            memset(buffer_, 0, nso1*nso2*sizeof(double));

            int nao2 = b2_->naofunction(jsh);

            fprintf(outfile, "ish = %d, jsh = %d\n", ish, jsh);

            // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
            // and then using the transformation coefficients to obtain everything else.
            // Otherwise using the petite list doesn't save us any computational time
            // in computing the integrals, but does save us time when we use the integrals.

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];

                if (!petite->in_p1(s1.aoshell)) {
                    fprintf(outfile, "skipping i = %d\n", i);
                    continue;
                }

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    double lambda = (double)petite->lambda(s1.aoshell, s2.aoshell);

                    fprintf(outfile, "i = %d, j = %d, in_p2 = %d lambda = %8.5f\n",
                            s1.aoshell, s2.aoshell, petite->in_p2(s1.aoshell, s2.aoshell),
                            lambda);

                    ob_->compute_shell(s1.aoshell, s2.aoshell);

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        int iirrep = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            int jirrep = jfunc.irrep;

                            buffer_[jsooff] += jcoef * aobuf[jaooff] * lambda;

                            if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
//                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
//                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
//                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
//                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
                                fprintf(outfile, "(%d|%d) += %8.5f * %8.5f * (%d|%d): %8.5f -> %8.5f\n",
                                        isofunc, jsofunc, lambda, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
                            }

                            // Check the irreps to ensure symmetric quantities.
                            if (ifunc.irrep == jfunc.irrep)
                                result->add(ifunc.irrep,
                                            b1_->function_within_irrep(ish, isofunc),
                                            b2_->function_within_irrep(jsh, jsofunc),
                                            jcoef * aobuf[jaooff] * lambda);
                        }
                    }
                }
            }
        }
    }
}

#if 0
void OneBodySOInt::compute_dcr(boost::shared_ptr<Matrix> result)
{
    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    shared_ptr<PetiteList> petite = b1_->petitelist();

    // Check symmetry.
    if (b1_->nirrep() == 1) {
        // Molecule has C1 symmetry. Need to use the ob_ object to get a SimpleMatrix.
        throw FeatureNotImplemented("libmints", "OneBodySOInt with C1 symmetry. Need to direct the code to OneBodyInt.", __FILE__, __LINE__);
    }

    for (int usi=0; usi<ns1; ++usi) {

        // Not sure if this check is needed.
        for (int usj=0; usj<=usi; ++usj) {

            const SOTransform &t1 = b1_->trans(usi);
            const SOTransform &t2 = b2_->trans(usj);

            int nso1 = b1_->nfunction(usi);
            int nso2 = b2_->nfunction(usj);

            memset(buffer_, 0, nso1*nso2*sizeof(double));

            int nao2 = b2_->naofunction(usj);

            fprintf(outfile, "ish = %d, jsh = %d\n", usi, usj);

            // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
            // and then using the transformation coefficients to obtain everything else.
            // Otherwise using the petite list doesn't save us any computational time
            // in computing the integrals, but does save us time when we use the integrals.

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];

                if (!petite->in_p1(s1.aoshell)) {
                    fprintf(outfile, "skipping i = %d\n", i);
                    continue;
                }

                int atom_i = integral_->basis1()->shell(s1.aoshell)->ncenter();
                unsigned int stab_i = petite->stabilizer(atom_i);

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    int atom_j = integral_->basis2()->shell(s2.aoshell)->ncenter();
                    unsigned int stab_j = petite->stabilizer(atom_j);

                    unsigned rOperator = petite->dcr(stab_i, stab_j);
                    double lambda = petite->order() / petite->dcr_degeneracy(rOperator);

                    fprintf(outfile, "i = %d, j = %d, in_p2 = %d stab_i = %d stab_j = %d lambda = %8.5f order = %d, degen = %d\n",
                            s1.aoshell, s2.aoshell, petite->in_p2(s1.aoshell, s2.aoshell),
                            stab_i, stab_j, lambda, petite->order(), petite->dcr_degeneracy(rOperator));

                    ob_->compute_shell(s1.aoshell, s2.aoshell);

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(usi, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        int iirrep = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(usj, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            int jirrep = jfunc.irrep;

                            buffer_[jsooff] += jcoef * aobuf[jaooff] * lambda;

                            //                            if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
                            //                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                            //                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                            //                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
                            //                                fprintf(outfile, "(%d|%d) += %8.5f * (%d|%d): %8.5f -> %8.5f\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
                            //                            }

                            // Check the irreps to ensure symmetric quantities.
                            if (ifunc.irrep == jfunc.irrep)
                                result->add(ifunc.irrep,
                                            b1_->function_within_irrep(usi, isofunc),
                                            b2_->function_within_irrep(usj, jsofunc),
                                            jcoef * aobuf[jaooff] * lambda);
                        }
                    }
                }
            }
        }
    }

    result->copy_lower_to_upper();
}
#endif

///////////////////////////////////////////////////////////////////////////////

TwoBodySOInt::TwoBodySOInt(const boost::shared_ptr<TwoBodyInt> &tb,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : tb_(tb), integral_(integral)
{
    // Try to reduce some work:
    b1_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis1(), integral));

    if (tb->basis1() == tb->basis2())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis2(), integral));

    if (tb->basis1() == tb->basis3())
        b3_ = b1_;
    else
        b3_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis3(), integral));

    if (tb->basis3() == tb->basis4())
        b4_ = b3_;
    else
        b4_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis4(), integral));

    // Allocate accumulation buffer
    buffer_ = new double[16*INT_NCART(tb->basis1()->max_am())
            *INT_NCART(tb->basis2()->max_am())
            *INT_NCART(tb->basis3()->max_am())
            *INT_NCART(tb->basis4()->max_am())];
}

TwoBodySOInt::~TwoBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis2() const
{
    return b2_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis3() const
{
    return b3_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis4() const
{
    return b4_;
}

}
