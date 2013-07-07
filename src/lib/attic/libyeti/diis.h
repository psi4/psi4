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

#ifndef yeti_diis_h
#define yeti_diis_h

#include "class.h"
#include "mapimpl.h"
#include <vector>
#include <list>

#include "diis.hpp"
#include "tensorparser.h"

namespace yeti {

class DiisParameters :
    public smartptr::Countable
{

    private:
        std::list<YetiTensorPtr> tensors_;

        std::list<YetiTensorPtr> residuals_;



    public:
        typedef std::list<YetiTensorPtr>::const_iterator iterator;

        DiisParameters();

        void add(
            const YetiTensorPtr& tensor,
            const YetiTensorPtr& residual
        );

        iterator tensor_begin() const;

        iterator tensor_end() const;

        iterator residual_begin() const;

        iterator residual_end() const;

        double dot_product(const DiisParametersPtr& params);



};

class DiisExtrapolation :
    public smartptr::Countable
{
    private:
        typedef std::list<DiisParametersPtr>::const_iterator iterator;

        std::list<DiisParametersPtr> params_;

        usi max_nparams_;

        double norm_;

        double* A;
        double* AP;
        int* ipiv;
        double* B;
        double* X;

        double* ferr;
        double* berr;
        double* work;
        int* iwork;

        void _extrapolate();

    public:

        DiisExtrapolation(usi max_nparams);

        ~DiisExtrapolation();

        void add(
            const YetiTensorPtr& t1,
            const YetiTensorPtr& r1,
            const YetiTensorPtr& t2 = 0,
            const YetiTensorPtr& r2 = 0,
            const YetiTensorPtr& t3 = 0,
            const YetiTensorPtr& r3 = 0,
            const YetiTensorPtr& t4 = 0,
            const YetiTensorPtr& r4 = 0
        );

        void clear();

        bool extrapolate(
            const YetiTensorPtr& t1,
            const YetiTensorPtr& t2 = 0,
            const YetiTensorPtr& t3 = 0,
            const YetiTensorPtr& t4 = 0
       );
};

}

#endif

