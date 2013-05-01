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

#ifndef SCHWARZ_H
#define SCHWARZ_H


#include "filler.h"
#include "index.h"
#include "aobasis.h"


namespace yeti {

class CauchySchwarzValueEstimater
        : public TensorValueEstimater
{
    protected:

        double** A;

        IndexDescr* descr_;

        IndexRange* range_;

        uli size_;

    public:

        CauchySchwarzValueEstimater(
            TEIShellComputeFunctorPtr tbint,
            AOBasisPtr aobasis,
            TensorIndexDescr* descr
        );

        CauchySchwarzValueEstimater(
            const CauchySchwarzValueEstimater* sub_est,
            TensorIndexDescr* descr,
            usi depth
        );

        ~CauchySchwarzValueEstimater();

        float max_log(const uli* indices) const;

};

} // end namespace yeti

#endif // SCHWARZ_H
