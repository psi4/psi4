#ifndef yeti_debugimpl_h
#define yeti_debugimpl_h

#include "tile.h"
#include "data.h"
#include "index.h"
#include "filler.h"

namespace yeti {


template <typename data_t>
class EmptyFiller : public TileElementComputer {

    private:
        int callcount_;

    public:
        EmptyFiller()
         :
            callcount_(0)
        {
        }

        TileElementComputer* copy() const
        {
            return new EmptyFiller;
        }

        int callcount() {return callcount_;}

};

template <typename data_t>
class ModulusFiller :
    public TileElementComputer
{

    private:
        size_t n_;

        size_t modulus_;

        int callcount_;

        uli offset_;

        uli denom_;

    public:
        ModulusFiller(
            uli modulus = 1000,
            uli offset = 0,
            uli denom = 0
        ) :
            callcount_(0),
            modulus_(modulus),
            n_(0),
            offset_(offset),
            denom_(denom)
        {
        }

        TileElementComputer* copy() const
        {
            return new ModulusFiller(modulus_, offset_, denom_);
        }

        void compute(Tile* tile, data_t* data)
        {
            uli ntot = 1;
            IndexRangeTuple::iterator it(tile->get_index_ranges()->begin());
            IndexRangeTuple::iterator stop(tile->get_index_ranges()->end());
            //nothing to compute
            for ( ; it != stop; ++it)
            {
                IndexRange* range(*it);
                ntot *= range->n();
            }

            data_t* dptr = data;
            for (size_t i=0; i < ntot; ++i, ++dptr)
            {
                if (denom_)
                {
                    double n = i % modulus_ + offset_;
                    double d = denom_;
                    (*dptr) = n / d;
                }
                else
                {
                    (*dptr) = i % modulus_ + offset_;
                }
            }
            ++callcount_;
        }

        int callcount() {return callcount_;}

};

}

#endif
