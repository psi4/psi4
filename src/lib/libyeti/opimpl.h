#ifndef yeti_OPIMPL_H
#define yeti_OPIMPL_H

#include "dataimpl.h"
#include "class.h"
#include "elementop.h"
#include "tile.h"

#include <libsmartptr/printstream.h>

namespace yeti {


template <typename scale_t>
class ScaleOp : public ElementOp {

    private:
        scale_t factor_;

        float maxlog_;

    public:

        ScaleOp(scale_t factor)
            : factor_(factor), maxlog_(0)
        {
            maxlog_ = log(fabs(factor));
        }

        void
        element_op(
            IndexRangeTuple* tuple,
            Data* data
        )
        {
            data_type_switch(data->type(), scale, data);
        }

        template <class data_t>
        void
        scale(
            Data* data
        )
        {
            data_t* dptr = data->get<data_t>();
            uli n = data->n();
            for (uli i=0; i < n; ++dptr, ++i)
                (*dptr) *= factor_;
        }

        void
        update(Tile *tile)
        {
            //adjust the tile's maxlog
            tile->scale_max_log(maxlog_);
        }

};

template <class data_t>
void
DijabOp::scale_by_denominator(
    Data* data
)
{
    data_t* dptr = data->get<data_t>();
    double ei, ej, ea, eb, dijab;
    for (uli i=istart_; i < istop_; ++i)
    {
        ei = evals_i_[i];
        for (uli j=jstart_; j < jstop_; ++j)
        {
            ej = evals_j_[j];
            for (uli a=astart_; a < astop_; ++a)
            {
                ea = evals_a_[a];
                for (uli b=bstart_; b < bstop_; ++b, ++dptr)
                {
                    eb = evals_b_[b];
                    dijab = ei + ej - ea - eb;
                    (*dptr) /= dijab;

                }
            }
        }
    }
}


}

#endif // OPIMPL_H
