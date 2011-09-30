#include "mints.h"
#include "moindexspace.h"

namespace psi {

MOIndexSpace::MOIndexSpace(const std::string& name,
                           const boost::shared_ptr<Matrix>& full_C,   // Should this be 4 C's instead?
                           const boost::shared_ptr<Vector>& evals,
                           const boost::shared_ptr<BasisSet>& basis,
                           const boost::shared_ptr<IntegralFactory>& ints)
    : nirrep_(full_C->nirrep()), name_(name), C_(full_C),
      evals_(evals), basis_(basis), ints_(ints)
{
}

MOIndexSpace::MOIndexSpace(const std::string &name, const boost::shared_ptr<Wavefunction> &wave)
    : nirrep_(wave->nirrep()), name_(name), C_(wave->Ca()), evals_(wave->epsilon_a()),
      basis_(wave->basisset()), ints_(wave->integral())
{
}

int MOIndexSpace::nirrep() const
{
    return nirrep_;
}

const std::string& MOIndexSpace::name() const
{
    return name_;
}

const boost::shared_ptr<Matrix>& MOIndexSpace::C() const
{
    return C_;
}

const boost::shared_ptr<Vector>& MOIndexSpace::evals() const
{
    return evals_;
}

const boost::shared_ptr<BasisSet>& MOIndexSpace::basis() const
{
    return basis_;
}

const boost::shared_ptr<IntegralFactory>& MOIndexSpace::integral() const
{
    return ints_;
}

const Dimension& MOIndexSpace::dim() const
{
    return dim_;
}

boost::shared_ptr<Matrix> MOIndexSpace::transform(const MOIndexSpace& space2, const MOIndexSpace& space1)
{
    boost::shared_ptr<Matrix> S21 = overlap(space2, space1);

    // Invert S21
    S21->general_invert();

    S21->set_name("Transformation Matrix between space1 and space2");

    return S21;
}

boost::shared_ptr<Matrix> MOIndexSpace::overlap(const MOIndexSpace &space2, const MOIndexSpace &space1)
{
    IntegralFactory mix_ints(space1.basis(), space2.basis(), space1.basis(), space2.basis());

    boost::shared_ptr<Matrix> Smat(new Matrix("Overlap between space1 and space2",
                                              space1.C()->rowspi(), space2.C()->colspi()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    boost::shared_ptr<Matrix> V1 = space1.C();
    boost::shared_ptr<Matrix> V2 = space2.C();

    boost::shared_ptr<Matrix> result(new Matrix("Overlap between space1 and space2", V1->colspi(), V2->colspi()));
    result->transform(V1, Smat, V2);

    return result;
}

} // namespace psi
