#include "filler.h"
#include "tile.h"
#include "data.h"
#include "mpqc.h"
#include "aobasis.h"

#if HAVE_MPQC

using namespace sc;
using namespace std;
using namespace yeti;
using namespace mpqc;

TwoElectronIntegralComputer::TwoElectronIntegralComputer(
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3,
    const Ref<GaussianBasisSet>& bs4
) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    bs4_(bs4),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    descr_(0)
{
    init();
}

TwoElectronIntegralComputer::TwoElectronIntegralComputer(
    const Ref<TwoBodyIntDescr>& descr,
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3,
    const Ref<GaussianBasisSet>& bs4
) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    bs4_(bs4),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    descr_(descr)
{
    init();
}

TileElementComputer*
TwoElectronIntegralComputer::copy() const
{
    return new TwoElectronIntegralComputer(IF_, bs1_, bs2_, bs3_, bs4_);
}

TwoElectronIntegralComputer::~TwoElectronIntegralComputer()
{
}

void
TwoElectronIntegralComputer::compute(Tile* tile, double* data)
{
    const size_t* indices = tile->indices();
    tbint_->compute_shell(indices[0], indices[1], indices[2], indices[3]);
    nblock_ = tile->get_data()->n();
    ::memcpy(data, buffer_, nblock_ * sizeof(double));
}

void
TwoElectronIntegralComputer::init()
{
    if (bs2_.null())
        bs2_ = bs1_;
    if (bs3_.null())
        bs3_ = bs2_;
    if (bs4_.null())
        bs4_ = bs3_;

    IF_->set_basis(bs1_, bs2_, bs3_, bs4_);

    if (descr_.nonnull())
        tbint_ = descr_->inteval();
    else
        tbint_ = IF_->electron_repulsion();

    buffer_ = tbint_->buffer();
}

DijabElementComputer::DijabElementComputer(
    const RefDiagSCMatrix& ei,
    const RefDiagSCMatrix& ej,
    const RefDiagSCMatrix& ea,
    const RefDiagSCMatrix& eb
)
    :
    ei_(ei),
    ej_(ej),
    ea_(ea),
    eb_(eb)
{
}

void
DijabElementComputer::compute(Tile* tile, double* data)
{
    IndexRangeTuple* tuple = tile->get_index_ranges();

    ni_ = tuple->get(0)->n();
    istart_ = tuple->get(0)->start();
    nj_ = tuple->get(1)->n();
    jstart_ = tuple->get(1)->start();
    na_ = tuple->get(2)->n();
    astart_ = tuple->get(2)->start();
    nb_ = tuple->get(3)->n();
    bstart_ = tuple->get(3)->start();

    uli istop = istart_ + ni_;
    uli jstop = jstart_ + nj_;
    uli astop = astart_ + na_;
    uli bstop = bstart_ + nb_;

    double* dataptr = data;
    for (size_t i=istart_; i < istop; ++i)
    {
        double ei = ei_.get_element(i);
        for (size_t j=jstart_; j <  jstop; ++j)
        {
            double ej = ej_.get_element(j);
            for (size_t a=astart_; a < astop; ++a)
            {
                double ea = ea_.get_element(a);
                for (size_t b=bstart_; b < bstop; ++b, ++dataptr)
                {
                    double eb = eb_.get_element(b);
                    (*dataptr) = ei + ej - ea - eb;
                }
            }
        }
    }
}

TileElementComputer*
DijabElementComputer::copy() const
{
    return new DijabElementComputer(ei_, ej_, ea_, eb_);
}

#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H
IntegralTransformElementComputer::IntegralTransformElementComputer(
    const Ref<IntegralTransform>& tform,
    uli offset_i,
    uli offset_j,
    uli offset_a,
    uli offset_b,
    uli nindex4
) : tform_(tform),
    nindex4_(nindex4),
    offset_i_(offset_i),
    offset_j_(offset_j),
    offset_a_(offset_a),
    offset_b_(offset_b)
{
}

void
IntegralTransformElementComputer::compute(Tile* tile, double* data)
{
    IndexRangeTuple* tuple = tile->get_index_ranges();

    ni_ = tuple->get(0)->n();
    istart_ = tuple->get(0)->start() - offset_i_;
    nj_ = tuple->get(1)->n();
    jstart_ = tuple->get(1)->start() - offset_j_;
    na_ = tuple->get(2)->n();
    astart_ = tuple->get(2)->start() - offset_a_;
    nb_ = tuple->get(3)->n();
    bstart_ = tuple->get(3)->start() - offset_b_;

    size_t istop = istart_ + ni_;
    size_t jstop = jstart_ + nj_;
    size_t astop = astart_ + na_;
    size_t bstop = bstart_ + nb_;

    double* dataptr = data;
    uli stride = nindex4_ - nb_;
    for (size_t i=istart_; i < istop; ++i)
    {
        for (size_t j=jstart_; j <  jstop; ++j)
        {
            const double* ints = tform_->retrieve_pair_block(i, j, TwoBodyInt::eri);
            const double* intptr = ints;
            intptr += astart_ * nindex4_ + bstart_;
            for (size_t a=astart_; a < astop; ++a, intptr += stride)
            {
                for (size_t b=bstart_; b < bstop; ++b, ++dataptr, ++intptr)
                {
                    (*dataptr) = (*intptr);
                }
            }
            tform_->release_pair_block(i,j,TwoBodyInt::eri);
        }
    }
}

TileElementComputer*
IntegralTransformElementComputer::copy() const
{
    return new IntegralTransformElementComputer(
        tform_,
        offset_i_,
        offset_j_,
        offset_a_,
        offset_b_,
        nindex4_
     );
}


IndexRange*
mpqc::get_ao_range(
    const Ref<GaussianBasisSet>& obs
)
{
    Ref<Molecule> mol = obs->molecule();
    //set up the yeti index range
    AOBasis aobasis(obs->label(), obs->name());
    for (int atom=0; atom < obs->ncenter(); ++atom)
    {
        AtomBasisPtr atombasis(new AtomBasis(mol->atom_symbol(atom), mol->atom_name(atom)));
        int start = obs->shell_on_center(atom,0);
        int stop = start + obs->nshell_on_center(atom);
        for (int shell=start; shell < stop; ++shell)
        {
            atombasis->add_shell(obs->shell(shell).nfunction());
        }
        aobasis.add_atom(atombasis);
    }
    IndexRange* aorange(aobasis.get_index_range());
    return aorange;
}


#endif




#endif
