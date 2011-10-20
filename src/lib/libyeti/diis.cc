#include "diis.h"
#include "gigmatrix.h"

#include "blas.h"
extern "C" {
extern void F_DSPSVX(const char* fact, const char* uplo, const int* n, const int* nrhs,
                       const double* AP, double* AFP, int* ipiv, const double* BB, const int* nb,
                       double* XX, const int* nx, double* rcond, double* ferr, double* berr,
                       double* work, int* iwork, int* info);
}//EndExternC

using namespace yeti;
using namespace std;


DiisParameters::DiisParameters()
{
}

void
DiisParameters::add(
    const YetiTensorPtr& tensor,
    const YetiTensorPtr& residual
)
{
    tensors_.push_back(tensor->copy());
    residuals_.push_back(residual->copy());
}

double
DiisParameters::dot_product(const DiisParametersPtr& params)
{
    iterator itl = residual_begin();
    iterator stop = residual_end();
    iterator itr = params->residual_begin();
    double E = 0;
    for ( ; itl != stop; ++itl, ++itr)
    {
        const YetiTensorPtr& l = *itl;
        const YetiTensorPtr& r = *itr;
        E += l * r;
    }
    return E;
}

DiisParameters::iterator
DiisParameters::tensor_begin() const
{
    return tensors_.begin();
}

DiisParameters::iterator
DiisParameters::tensor_end() const
{
    return tensors_.end();
}

DiisParameters::iterator
DiisParameters::residual_begin() const
{
    return residuals_.begin();
}

DiisParameters::iterator
DiisParameters::residual_end() const
{
    return residuals_.end();
}

DiisExtrapolation::DiisExtrapolation(usi max_nparams)
    : 
    max_nparams_(max_nparams),
    norm_(0),
    ferr(0),
    berr(0),
    work(0),
    iwork(0),
    A(0),
    AP(0),
    ipiv(0),
    B(0),
    X(0)
{
    int nparams = max_nparams;
    int n = nparams + 1;

    // convert A to packed upper triangular form
    int ntri = n*(n+1)/2;
    A = new double[ntri];
    AP = new double[ntri];
    ipiv = new int[n];
    B = new double[n];
    X = new double[n];
    ferr = new double[1];
    berr = new double[1];
    work = new double[3*n];
    iwork = new int[n];
}

DiisExtrapolation::~DiisExtrapolation()
{
    delete[] A;
    delete[] AP;
    delete[] ipiv;
    delete[] B;
    delete[] X;
    delete[] ferr;
    delete[] berr;
    delete[] work;
    delete[] iwork;
}

void
DiisExtrapolation::add(
    const YetiTensorPtr& t1,
    const YetiTensorPtr& r1,
    const YetiTensorPtr& t2,
    const YetiTensorPtr& r2,
    const YetiTensorPtr& t3,
    const YetiTensorPtr& r3,
    const YetiTensorPtr& t4,
    const YetiTensorPtr& r4
)
{
    if (params_.size() >= max_nparams_)
        params_.pop_front();

    DiisParametersPtr params = new DiisParameters;
    params->add(t1,r1);
    if (t2) params->add(t2,r2);
    if (t3) params->add(t3,r3);
    if (t4) params->add(t4,r4);

    params_.push_back(params);
}

void
DiisExtrapolation::_extrapolate()
{
    int nparams = params_.size();
    int n = nparams + 1;

    iterator itrow = params_.begin();

    // convert A to packed upper triangular form
    int ntri = n*(n+1)/2;
    double* aptr = A;
    double norm = 1;
    for (int row=0; row < nparams + 1; ++row, ++itrow)
    {
        const DiisParametersPtr& diisrow = *itrow;
        iterator itcol = params_.begin();
        for (int col=0; col <= row; ++col, ++itcol, ++aptr)
        {
            if (row == nparams && col == nparams)
            {
                *aptr = 0;
            }
            else if (row == nparams)
            {
                *aptr = -1;
            }
            else
            {
                const DiisParametersPtr& diiscol = *itcol;
                double E = diisrow->dot_product(diiscol);
                if (row == 0 && col == 0)
                    norm = E;
                *aptr = E / norm;
            }
        }
    }

    for (usi i=0; i < nparams; ++i)
        B[i] = 0;
    B[nparams] = -1;


    char fact = 'N';
    char uplo = 'U';
    int nrhs = 1;
    double rcond = 0.0;
    int info = 0;

    F_DSPSVX(&fact, &uplo, &n, &nrhs, A, AP, ipiv, B, &n, X, &n, &rcond, ferr, berr, work, iwork, &info);
}

bool
DiisExtrapolation::extrapolate(
    const YetiTensorPtr& t1,
    const YetiTensorPtr& t2,
    const YetiTensorPtr& t3,
    const YetiTensorPtr& t4
)
{
    int nparams = params_.size();
    int n = nparams + 1;

    int nerr = 1;
    while (nerr)
    {
        _extrapolate();

        nerr = 0;
        for (int i=0; i < nparams; ++i)
        {
            double err = ferr[i];
            if ( fabs(err) > 0.1 )
                ++nerr;
        }

        if (nerr) //throw out the oldest one
        {
            cout << "reduced DIIS dimensionality" << endl;
            params_.pop_front();
        }
    }

    nparams = params_.size();
    n = nparams + 1;

    if (nparams == 1) //failed!
        return false;


    iterator itparams = params_.begin();

    t1->zero();
    if (t2) t2->zero();
    if (t3) t3->zero();
    if (t4) t4->zero();

    for (usi i=0; i < nparams; ++i, ++itparams)
    {
        DiisParametersPtr params = *itparams;
        double coef = X[i];

        DiisParameters::iterator it = params->tensor_begin();

        if (t1)
        {
            const YetiTensorPtr& _t = *it;
            t1->accumulate(_t, coef);
        }
        if (t2)
        {
            ++it;
            const YetiTensorPtr& _t = *it;
            t2->accumulate(_t, coef);
        }
        if (t3)
        {
            ++it;
            const YetiTensorPtr& _t = *it;
            t3->accumulate(_t, coef);
        }
        if (t4)
        {
            ++it;
            const YetiTensorPtr& _t = *it;
            t4->accumulate(_t, coef);
        }

    }
}
