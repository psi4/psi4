#ifndef yeti_ELEMENTOP_H
#define yeti_ELEMENTOP_H

#include "class.h"

#include "data.hpp"
#include "index.hpp"
#include "elementop.hpp"
#include "tile.hpp"

namespace yeti {

class ElementOp : public smartptr::Countable {

    public:
        virtual void element_op(
            IndexRangeTuple* tuple,
            Data* data
        ) = 0;

        virtual void update(Tile* tile);

};

class DijabOp : public ElementOp {

    private:
        const double* evals_i_;
        const double* evals_j_;
        const double* evals_a_;
        const double* evals_b_;
        uli istart_;
        uli istop_;
        uli jstart_;
        uli jstop_;
        uli astart_;
        uli astop_;
        uli bstart_;
        uli bstop_;

    public:
        DijabOp(
            const double* ei,
            const double* ej,
            const double* ea,
            const double* eb
        );

        void
        element_op(
            IndexRangeTuple* tuple,
            Data* data
        );

        template <class data_t>
        void
        scale_by_denominator(
            Data* data
        );

};

} //namespace yeti

#endif // ELEMENTOP_H
