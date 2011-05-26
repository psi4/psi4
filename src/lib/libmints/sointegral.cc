#include <psi4-dec.h>

#include "sointegral.h"
#include "twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"
#include "sobasis.h"
#include "matrix.h"
#include "molecule.h"

#include <boost/shared_ptr.hpp>

#define DEBUG

namespace psi {

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyAOInt> & ob,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : ob_(ob), integral_(integral.get()), deriv_(ob->deriv())
{
    common_init();

}

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyAOInt> & ob,
                           const IntegralFactory* integral)
    : ob_(ob), integral_(integral), deriv_(ob->deriv())
{
    common_init();
}

OneBodySOInt::~OneBodySOInt()
{
    delete[] buffer_;
}

void OneBodySOInt::common_init()
{

    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis1(), integral_));

    if (ob_->basis2() == ob_->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis2(), integral_));

    ob_->set_force_cartesian(b1_->petitelist()->include_pure_transform());

    int max1_nequivalent_atoms = ob_->basis1()->molecule()->max_nequivalent();
    int max2_nequivalent_atoms = ob_->basis2()->molecule()->max_nequivalent();

    size_ = max1_nequivalent_atoms*max2_nequivalent_atoms*
            INT_NCART(ob_->basis1()->max_am())*INT_NCART(ob_->basis2()->max_am());

    if (deriv_ == 1) {
        // Make sure the AOInt knows how to compute derivatives
        if (!ob_->has_deriv1())
            throw PSIEXCEPTION("OneBodySOInt::OneBodySOInt: The AO integral object doesn't provide first derivatives.");
        size_ *= 3 * ob_->basis1()->molecule()->natom();
    }
    if (deriv_ > 1)
        throw FeatureNotImplemented("libmints", "Symmetrized integral derivatives greater than first order not implemented.",
                                    __FILE__, __LINE__);

    // Grab the number of chunks from the ao object and enlarge the buffer to handle
    size_ *= ob_->nchunk();

    buffer_ = new double[size_];
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis2() const
{
    return b2_;
}

void OneBodySOInt::compute(boost::shared_ptr<Matrix> result)
{
    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);

            memset(buffer_, 0, size_*sizeof(double));

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
//                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d jaoff = %d jsooff = %d\n",
//                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
//                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
//                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc),
//                                        jaooff, jsooff);
//                                fprintf(outfile, "(%d|%d) += %8.5f * (%d|%d): %8.5f -> %8.5f\n",
//                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
//                                fflush(outfile);
//                            }

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

void OneBodySOInt::compute(std::vector<boost::shared_ptr<Matrix> > results)
{
    // Do not worry about zeroing out result
    int nchunk = ob_->nchunk();
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);
            int nso = nso1*nso2;

            memset(buffer_, 0, nchunk*nso*sizeof(double));

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

                            // Handle chunks
                            for (int i=0; i<nchunk; ++i) {
                                double temp = jcoef * aobuf[jaooff + (i*nao)];

                                // If I get the if statements working below the next line is not needed at all.
                                buffer_[jsooff + (i*nso)] += temp;

                                int ijirrep = ifunc.irrep ^ jfunc.irrep;
                                if (ijirrep == results[i]->symmetry()) {

//                                    if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
//                                        fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d jaoff = %d jsooff = %d\n",
//                                                isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff + (i*nao)], buffer_[jsooff + (i*nso)],
//                                                ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
//                                                jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc),
//                                                jaooff + (i*nao), jsooff + (i*nso));
//                                        fflush(outfile);
//                                    }

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

void OneBodySOInt::compute_deriv1(std::vector<boost::shared_ptr<Matrix> > result,
                                  const CdSalcList &cdsalcs)
{
    // Do not worry about zeroing out result.

    // Do some checks:
    if (deriv_ < 1)
        throw SanityCheckError("OneBodySOInt::compute_deriv1: integral object not created to handle derivatives.", __FILE__, __LINE__);

    if (result.size() != cdsalcs.ncd())
        throw SanityCheckError("OneBodySOInt::compute_deriv1: result vector size does not match SALC size.", __FILE__, __LINE__);

    Molecule& mol = *ob_->basis1()->molecule().get();
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over unique SO shells.
    for (int ish=0; ish<ns1; ++ish) {
        const SOTransform& t1 = b1_->trans(ish);
        int nso1 = b1_->nfunction(ish);
        int nao1 = b1_->naofunction(ish);

        for (int jsh=0; jsh<ns2; ++jsh) {
            const SOTransform& t2= b2_->trans(jsh);
            int nso2 = b2_->nfunction(jsh);
            int nao2 = b2_->naofunction(jsh);

            int nao12 = nao1 * nao2;
            int nso12 = nso1 * nso2;

            // Clear out the memory we need.
            memset(buffer_, 0, size_*sizeof(double));

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                int center_i = ob_->basis1()->shell(s1.aoshell)->ncenter();
                const CdSalcWRTAtom& cdsalc1 = cdsalcs.atom_salc(center_i);

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];
                    int center_j = ob_->basis2()->shell(s2.aoshell)->ncenter();
                    const CdSalcWRTAtom& cdsalc2 = cdsalcs.atom_salc(center_j);

                    // If we're working on the same atomic center, don't even bother with the derivative
                    if (center_i == center_j)
                        continue;

                    ob_->compute_shell_deriv1(s1.aoshell, s2.aoshell);

                    // Need to loop over the cdsalcs

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
                        // SO function offset in a lienar array
                        int isooff   = isofunc;
                        // Relative position of the SO function within its irrep
                        int irel     = b1_->function_within_irrep(ish, isofunc);
                        int iirrep   = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc  = jfunc.aofunc;
                            int jsofunc  = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff   = iaooff*nao2 + jaofunc;
                            int jsooff   = isooff*nso2 + jsofunc;
                            int jrel     = b2_->function_within_irrep(jsh, jsofunc);
                            int jirrep   = jfunc.irrep;

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

TwoBodySOInt::TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt> &tb,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : tb_(tb), integral_(integral)
{
    // Try to reduce some work:
    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis1(), integral));

    if (tb->basis1() == tb->basis2())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis2(), integral));

    if (tb->basis1() == tb->basis3())
        b3_ = b1_;
    else
        b3_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis3(), integral));

    if (tb->basis3() == tb->basis4())
        b4_ = b3_;
    else
        b4_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis4(), integral));

    tb_->set_force_cartesian(b1_->petitelist()->include_pure_transform());

    size_ = b1_->max_nfunction_in_shell() *
            b2_->max_nfunction_in_shell() *
            b3_->max_nfunction_in_shell() *
            b4_->max_nfunction_in_shell();
    buffer_ = new double[size_];

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
}

TwoBodySOInt::~TwoBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis2() const
{
    return b2_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis3() const
{
    return b3_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis4() const
{
    return b4_;
}

}
