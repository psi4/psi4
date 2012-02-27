#include <boost/shared_ptr.hpp>
#include "mints.h"
#include "orbitalspace.h"
#include "orthog.h"

namespace psi {

OrbitalSpace::OrbitalSpace(const std::string& id,
                           const std::string &name,   // Should this be 4 C's instead?
                           const SharedMatrix &full_C,
                           const boost::shared_ptr<Vector> &evals,
                           const boost::shared_ptr<BasisSet> &basis,
                           const boost::shared_ptr<IntegralFactory> &ints)
    : id_(id),
      name_(name),
      C_(full_C),
      evals_(evals),
      basis_(basis),
      ints_(ints),
      dim_(full_C->colspi())
{
}

OrbitalSpace::OrbitalSpace(const std::string& id,
                           const std::string& name,
                           const SharedMatrix& full_C,
                           const boost::shared_ptr<BasisSet>& basis,
                           const boost::shared_ptr<IntegralFactory>& ints)
    : id_(id),
      name_(name),
      C_(full_C),
      basis_(basis),
      ints_(ints),
      dim_(full_C->colspi())
{
}

OrbitalSpace::OrbitalSpace(const std::string &id,
                           const std::string &name,
                           const boost::shared_ptr<Wavefunction> &wave)
    : id_(id),
      name_(name),
      C_(wave->Ca()),
      evals_(wave->epsilon_a()),
      basis_(wave->basisset()),
      ints_(wave->integral()),
      dim_(wave->Ca()->colspi())
{
}

int OrbitalSpace::nirrep() const
{
    return C_->nirrep();
}

const std::string& OrbitalSpace::id() const
{
    return id_;
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

    return OrbitalSpace("p",
                        "Ca transformed into Cb",
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

void OrbitalSpace::print() const
{
    fprintf(outfile, "    Orbital space %s (%s)\n", name_.c_str(), id_.c_str());
    fprintf(outfile, "        Basis: %s\n", basis_->name().c_str());
    fprintf(outfile, "        Dimensions: "); dim_.print();
    fprintf(outfile, "        Transformation matrix:\n");
    C_->print();
}

namespace SpaceBuilder
{
namespace { // anonymous
    OrbitalSpace orthogonalize(const std::string& id, const std::string& name,
                               const boost::shared_ptr<BasisSet>& bs,
                               const boost::shared_ptr<IntegralFactory>& ints,
                               double lindep_tol)
    {
        OneBodySOInt *o_engine = ints->so_overlap();

        SOBasisSet* so_bs = new SOBasisSet(bs, ints);
        const Dimension& SODIM = so_bs->petite_list()->SO_basisdim();
        delete so_bs;

        SharedMatrix overlap(new Matrix("Overlap", SODIM, SODIM));
        o_engine->compute(overlap);
        delete o_engine;

        fprintf(outfile, "    Orthogonalizing basis for space %s:\n", name.c_str());
        OverlapOrthog orthog(OverlapOrthog::Symmetric,
                             overlap,
                             lindep_tol,
                             0);

        SharedMatrix orthog_so = orthog.basis_to_orthog_basis();
        fprintf(outfile, "    Minimum eigenvalue in the overlap matrix is %14.10E.\n\n", orthog.min_orthog_res());

        orthog_so->transpose_this();

        PetiteList petite(bs, ints);
        SharedMatrix orthog_ao = petite.evecs_to_AO_basis(orthog_so);

        return OrbitalSpace(id, name, orthog_ao, bs, ints);
    }

    OrbitalSpace orthogonal_compliment(const OrbitalSpace& space1, const OrbitalSpace& space2, const std::string& id, const std::string& name, double lindep_tol)
    {
        fprintf(outfile, "    Projecting out '%s' from '%s' to obtain space '%s'\n",
                space1.name().c_str(), space2.name().c_str(), name.c_str());

        // If space1 is empty, return a copy of the original space.
        if (space1.dim().sum() == 0)
            return OrbitalSpace(id, name, space2.C(), space2.evals(), space2.basis(), space2.integral());

        // C12 = C1 * S12 * C2
        SharedMatrix C12 = OrbitalSpace::overlap(space1, space2);

        C12->print();
    }
} // namespace anonymous

    OrbitalSpace build_cabs_space(const OrbitalSpace &orb_space, const OrbitalSpace &ri_space, double lindep_tol)
    {
        return orthogonal_compliment(orb_space, ri_space, "p''", "CABS", lindep_tol);
    }

    OrbitalSpace build_ri_space(boost::shared_ptr<BasisSet> aux_bs, boost::shared_ptr<BasisSet> obs, boost::shared_ptr<IntegralFactory> ints, double lindep_tol)
    {
        boost::shared_ptr<BasisSet> ri_basis = obs + aux_bs;
        return orthogonalize("p'", "RIBS", ri_basis, ints, lindep_tol);
    }
} // namespace SpaceBuilder

} // namespace psi
