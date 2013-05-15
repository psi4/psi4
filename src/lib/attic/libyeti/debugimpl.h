/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef yeti_debugimpl_h
#define yeti_debugimpl_h

#include "data.h"
#include "index.h"
#include "filler.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {


template <typename data_t>
class EmptyFiller :
    public TensorElementComputer
{

    private:
        int callcount_;

    public:
        EmptyFiller()
         :
            callcount_(0)
        {
        }

        TensorElementComputer* copy() const
        {
            return new EmptyFiller;
        }

        int callcount() {return callcount_;}

};

template <typename data_t>
class ModulusFiller :
    public TensorElementComputer
{

    private:
        uli n_;

        uli modulus_;

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

        TensorElementComputer* copy() const
        {
            return new ModulusFiller(modulus_, offset_, denom_);
        }

        void compute(const uli* indices, data_t* data, uli n)
        {
            data_t* dptr = data;
            for (uli i=0; i < n; ++i, ++dptr)
            {
                if (denom_)
                {
                    double numer = i % modulus_ + offset_;
                    double denom = denom_;
                    (*dptr) = numer / denom;
                }
                else
                {
                    (*dptr) = i % modulus_ + offset_;
                }
            }
            ++callcount_;
        }

        int callcount() {return callcount_;}

        TemplateInfo::type_t element_type(const uli *indices, usi depth)
        {
            return yeti::TypeInfo<data_t>::type();
        }

};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
