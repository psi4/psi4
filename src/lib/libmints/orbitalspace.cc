#include "mints.h"
#include "orbitalspace.h"

namespace psi {

OrbitalSpace::OrbitalSpace(const std::string& name,
                           const SharedMatrix& full_C,   // Should this be 4 C's instead?
                           const boost::shared_ptr<Vector>& evals,
                           const boost::shared_ptr<BasisSet>& basis,
                           const boost::shared_ptr<IntegralFactory>& ints)
    : nirrep_(full_C->nirrep()), name_(name), C_(full_C),
      evals_(evals), basis_(basis), ints_(ints), dim_(full_C->colspi())
{
}

OrbitalSpace::OrbitalSpace(const std::string &name, const boost::shared_ptr<Wavefunction> &wave)
    : nirrep_(wave->nirrep()), name_(name), C_(wave->Ca()), evals_(wave->epsilon_a()),
      basis_(wave->basisset()), ints_(wave->integral()), dim_(wave->Ca()->colspi())
{
}

int OrbitalSpace::nirrep() const
{
    return nirrep_;
}

const std::string& OrbitalSpace::name() const
{
    return name_;
}

const SharedMatrix& OrbitalSpace::C() const
{
    return C_;
}

const boost::shared_ptr<Vector>& OrbitalSpace::evals() const
{
    return evals_;
}

const boost::shared_ptr<BasisSet>& OrbitalSpace::basis() const
{
    return basis_;
}

const boost::shared_ptr<IntegralFactory>& OrbitalSpace::integral() const
{
    return ints_;
}

const Dimension& OrbitalSpace::dim() const
{
    return dim_;
}

OrbitalSpace OrbitalSpace::transform(const OrbitalSpace& A, const boost::shared_ptr<BasisSet>& B)
{
    SharedMatrix SBA = overlap(B, A.basis());
    SBA->set_name("Sba");
    SharedMatrix SBB = overlap(B, B);
    SBB->set_name("SBB");

    // Follows Werner's method from Mol. Phys. 102, 21-22, 2311
    // just like HF::dualBasisProjection

    // 1. Invert SBB
    SBB->invert();
    SBB->set_name("SBB^-1");

    // 2. Form T
    SharedMatrix I = Matrix::create("I = SAB SBB SBA", SBA->colspi(), SBA->colspi());
    I->transform(SBB, SBA);

    SharedMatrix T = Matrix::create("T", A.dim(), A.dim());
    T->transform(I, A.C());
    I.reset(); // release memory

    // 3. Form T^{-1/2}
    T->power(-0.5);

    // 4. Cb = [Sbb]^-1 Sba] Ca T^{-1/2}
    // 4a. Ca T^{-1/2}
    SharedMatrix CaT = Matrix::create("Ca*T^{-1/2}", A.C()->rowspi(), A.C()->colspi());
    CaT->gemm(false, false, 1.0, A.C(), T, 0.0);

    // 4b. Sba * 4a
    SharedMatrix SbaCaT = Matrix::create("SbaCaT", SBB->rowspi(), A.C()->colspi());
    SbaCaT->gemm(false, false, 1.0, SBA, CaT, 0.0);

    // 4c. [Sbb]^-1 * 4b
    SharedMatrix Cb = Matrix::create("Cb", SBB->rowspi(), A.C()->colspi());
    Cb->gemm(false, false, 1.0, SBB, SbaCaT, 0.0);

    boost::shared_ptr<IntegralFactory> i(new IntegralFactory(B, B, B, B));

    return OrbitalSpace("Ca transformed into Cb",
                        Cb,
                        A.evals(),
                        B,
                        i);
}

SharedMatrix OrbitalSpace::overlap(const OrbitalSpace &space1, const OrbitalSpace &space2)
{
    IntegralFactory mix_ints(space1.basis(), space2.basis(), space1.basis(), space2.basis());

    SharedMatrix Smat(new Matrix("Overlap between space1 and space2",
                                 space1.C()->rowspi(), space2.C()->colspi()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    return Smat;
}

SharedMatrix OrbitalSpace::overlap(const boost::shared_ptr<BasisSet>& basis1,
                                   const boost::shared_ptr<BasisSet>& basis2)
{
    IntegralFactory mix_ints(basis1, basis2);
    SOBasisSet sobasis1(basis1, &mix_ints);
    SOBasisSet sobasis2(basis2, &mix_ints);

    SharedMatrix Smat(new Matrix("Overlap between space1 and space2",
                                 sobasis1.dimension(), sobasis2.dimension()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    return Smat;
}

} // namespace psi
