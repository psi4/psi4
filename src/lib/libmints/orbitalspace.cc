#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "mints.h"
#include "view.h"
#include "orbitalspace.h"
#include "orthog.h"

#include <boost/tuple/tuple.hpp>

namespace psi {

OrbitalSpace::OrbitalSpace(const std::string& id,
                           const std::string &name,
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

const boost::shared_ptr<BasisSet>& OrbitalSpace::basisset() const
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
    SharedMatrix SBA = overlap(B, A.basisset());
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
    IntegralFactory mix_ints(space1.basisset(), space2.basisset());

    PetiteList p1(space1.basisset(), space1.integral());
    PetiteList p2(space2.basisset(), space2.integral());

    SharedMatrix Smat(new Matrix("Overlap between space1 and space2",
                                 p1.SO_basisdim(), p2.SO_basisdim()));

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

namespace { // anonymous
    OrbitalSpace orthogonalize(const std::string& id, const std::string& name,
                               const boost::shared_ptr<BasisSet>& bs,
                               double lindep_tol)
    {
        boost::shared_ptr<IntegralFactory> localfactory(new IntegralFactory(bs));
        OneBodySOInt *o_engine = localfactory->so_overlap();

        SOBasisSet *so_bs = new SOBasisSet(bs, localfactory);
        const Dimension& SODIM = so_bs->petite_list()->SO_basisdim();
        delete so_bs;

        SharedMatrix overlap(new Matrix("Overlap", SODIM, SODIM));
        o_engine->compute(overlap);
        delete o_engine;

        fprintf(outfile, "    Orthogonalizing basis for space %s:\n", name.c_str());
        Dimension remaining = overlap->power(-0.5, lindep_tol);

        View Cview(overlap, overlap->rowspi(), remaining);
        SharedMatrix C = Cview();
        C->set_name("Transformation matrix");

        PetiteList petite(bs, localfactory);
        SharedMatrix orthog_ao = petite.evecs_to_AO_basis(C);

        return OrbitalSpace(id, name, orthog_ao, bs, localfactory);
    }

    OrbitalSpace orthogonal_compliment(const OrbitalSpace& space1, const OrbitalSpace& space2, const std::string& id, const std::string& name, const double& lindep_tol)
    {
        fprintf(outfile, "    Projecting out '%s' from '%s' to obtain space '%s'\n",
                space1.name().c_str(), space2.name().c_str(), name.c_str());

        // If space1 is empty, return a copy of the original space.
        if (space1.dim().sum() == 0)
            return OrbitalSpace(id, name, space2.C(), space2.evals(), space2.basisset(), space2.integral());

        // O12 = O12
        SharedMatrix O12 = OrbitalSpace::overlap(space1, space2);
        SharedMatrix C12 = Matrix::create("C12", space1.C()->colspi(), space2.C()->colspi());

        O12->print();
        space1.C()->print();
        space2.C()->print();

        // C12 = C1t * S12 * C2
        C12->transform(space1.C(), O12, space2.C());
//        Matrix temp("temp", O12->rowspi(), space2.C()->colspi());
//        temp.gemm(false, false, 1.0, O12, space2.C(), 0.0);
//        temp.print();
//        C12->gemm(true, false, 1.0, space1.C(), temp, 0.0);
        C12->print();

        // SVD C12 =
//        boost::tuple<SharedMatrix, SharedVector, SharedMatrix> svd_temps = C12->svd_temps();
//        SharedMatrix U = svd_temps.get<0>();
//        SharedVector Sigma = svd_temps.get<1>();
//        SharedMatrix V = svd_temps.get<2>();

        Dimension smallest(C12->nirrep());
        const Dimension& rowd = C12->rowspi();
        const Dimension& cold = C12->colspi();

        for (int h=0, nirrep=C12->nirrep(); h<nirrep; ++h)
            smallest[h] = rowd[h] < cold[h] ? rowd[h] : cold[h];

        SharedMatrix U = Matrix::create("U", rowd, rowd);
        SharedMatrix V = Matrix::create("V", cold, cold);
        SharedVector Sigma = Vector::create("Sigma", smallest);

        // Something is not right with our SVD call
        // transpose our C12 to become fortran like
        SharedMatrix C12fortran = C12->transpose();

        for (int h=0; h<C12fortran->nirrep(); ++h) {
            if (!C12fortran->rowspi(h) || !C12fortran->colspi(h))
                continue;

            int m = C12->rowspi(h);
            int n = C12->colspi(h);
            int k = (m < n ? m : n);

            double** Ap = block_matrix(m,n);
            ::memcpy((void*) Ap[0], (void*) C12fortran->pointer(h)[0], sizeof(double) * m * n);
            double*  Sp = Sigma->pointer(h);
            double** Up = U->pointer(h);
            double** Vp = V->pointer(h);

            int* iwork = new int[8L * k];

            // Workspace Query
            double lwork;
            int info = C_DGESDD('A',m,n,Ap[0],m,Sp,Up[0],k,Vp[0],n,&lwork,-1,iwork);

            double* work = new double[(int)lwork];

            // SVD
            info = C_DGESDD('A',m,n,Ap[0],m,Sp,Up[0],k,Vp[0],n,work,(int)lwork,iwork);

            delete[] work;
            delete[] iwork;

            if (info != 0) {
                if (info < 0) {
                    fprintf(outfile, "Matrix::svd with metric: C_DGESDD: argument %d has invalid parameter.\n", -info);
                    fflush(outfile);
                    abort();
                }
                if (info > 0) {
                    fprintf(outfile, "Matrix::svd with metric: C_DGESDD: error value: %d\n", info);
                    fflush(outfile);
                    abort();
                }
            }
            free_block(Ap);
        }
        // Computes A
//        C12->svd_a(U, Sigma, V);
        U->print();
        Sigma->print();
        V->print();

        // No need to transpose since we transposed C12
//        V = V->transpose();
        SharedMatrix Vao = Matrix::create("Vao", V->rowspi(), space2.C()->rowspi());
        Vao->gemm(false, true, 1.0, V, space2.C(), 0.0);
        Vao->print();

        int nlindep=0;
        double min_sigma = 1.0;
        double max_sigma = 0.0;

        Dimension remove(Sigma->nirrep(), "Number of orbitals to remove");

        // Walk through S and determine how many to keep
        for (int h=0, nirrep=Sigma->nirrep(); h<nirrep; ++h) {
            int nzeros=0;
            int nsigma=Sigma->dim(h);

            for (int s=0; s<nsigma; ++s) {
                double sigma = Sigma->get(h, s);

                if (sigma < lindep_tol)
                    nzeros++;
                if (sigma < min_sigma)
                  min_sigma = sigma;
                if (sigma > max_sigma)
                  max_sigma = sigma;
            }

            nlindep += nzeros;
            remove[h] = nsigma - nzeros;
        }

        remove.print();
        Dimension zero(Vao->nirrep());
        View removes(Vao, Vao->rowspi()-remove, Vao->colspi(), remove, zero);

        SharedMatrix C = removes()->transpose();
        C->set_name("Final");
        C->print();
        return OrbitalSpace(id, name, C, space2.basisset(), space2.integral());
    }
} // namespace anonymous

    OrbitalSpace OrbitalSpace::build_cabs_space(const OrbitalSpace &orb_space, const OrbitalSpace &ri_space, double lindep_tol)
    {
        return orthogonal_compliment(orb_space, ri_space, "p''", "CABS", lindep_tol);
    }

    OrbitalSpace OrbitalSpace::build_ri_space(boost::shared_ptr<BasisSet> aux_bs, boost::shared_ptr<BasisSet> obs, double lindep_tol)
    {
        boost::shared_ptr<BasisSet> ri_basis = obs + aux_bs;
        return orthogonalize("p'", "RIBS", ri_basis, lindep_tol);
    }

} // namespace psi
