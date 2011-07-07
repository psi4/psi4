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

