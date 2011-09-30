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

}

boost::shared_ptr<Matrix> MOIndexSpace::overlap(const MOIndexSpace &space2, const MOIndexSpace &space1)
{
    IntegralFactory mix_ints(space1.basis(), space2.basis(), space1.basis(), space2.basis());

    boost::shared_ptr<Matrix> Smat(new Matrix("Overlap between space1 and space2",
                                              space1.C()->rowspi(), space2.C()->colspi()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    Smat->print();

    boost::shared_ptr<Matrix> result(new Matrix("Transformation space1 to space2",
                                                space2.C()->colspi(), space1.C()->colspi()));
    result->transform(space1.C(), Smat, space2.C());

    return result;
}

} // namespace psi
