#include "elementop.h"
#include "opimpl.h"

using namespace yeti;

void
ElementOp::update(Tile *tile)
{
    //do nothing by default
}

DijabOp::DijabOp(
    const double *ei,
    const double *ej,
    const double *ea,
    const double *eb
)
    :
    evals_i_(ei),
    evals_j_(ej),
    evals_a_(ea),
    evals_b_(eb)
{
}

void
DijabOp::element_op(
    IndexRangeTuple *tuple,
    Data *data
)
{
    istart_ = tuple->get(0)->start();
    istop_ = tuple->get(0)->n() + istart_;
    jstart_ = tuple->get(1)->start();
    jstop_ = tuple->get(1)->n() + jstart_;
    astart_ = tuple->get(2)->start();
    astop_ = tuple->get(2)->n() + astart_;
    bstart_ = tuple->get(3)->start();
    bstop_ = tuple->get(3)->n() + bstart_;

    data_type_switch(data->type(), scale_by_denominator, data);
}
