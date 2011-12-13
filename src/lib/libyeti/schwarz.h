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
