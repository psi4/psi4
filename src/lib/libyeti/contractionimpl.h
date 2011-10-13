#ifndef yeti_contraction_impl_h
#define yeti_contraction_impl_h

#include "data.h"
#include "contraction.h"
#include "tensor.h"
#include "runtime.h"

#include <libsmartptr/printstream.h>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define NO_DGEMM 0
#define DGEMM_CUTOFF 100

#include "blas.h"
extern "C" {

extern void F_DGEMM(const char*, const char*, const int*,
  const int*, const int*, const double*, const double*, const int*,
  const double*, const int*, const double*, double*, const int*);

}//EndExternC

namespace yeti {


#define class_type_switch(data_type, cls, fxn, ...) \
    switch(data_type) \
    { \
    case TemplateInfo::double_type: \
        cls<double>::fxn(__VA_ARGS__); \
        break; \
    case TemplateInfo::integer_type: \
        cls<int>::fxn(__VA_ARGS__); \
        break; \
    case TemplateInfo::float_type: \
        cls<float>::fxn(__VA_ARGS__); \
        break; \
    case TemplateInfo::quad_type: \
        cls<quad>::fxn(__VA_ARGS__); \
        break; \
    }

/**
    @class Contraction_tn Contraction. Both link indices are row indices.
            C = A^T * B
*/
template <typename data_t>
struct __Contraction_tn {

        static void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            data_t scale
        )
        {
            if (YetiRuntime::print_cxn)
                std::cout << "TN" << std::endl;
            //now accumulate the contraction
            for (uli link=0; link < nlink; ++link, ldata += nrows, rdata += ncols)
            {
                data_t* lptr = ldata;
                data_t* pptr = pdata;
                for (uli row=0; row < nrows; ++row, ++lptr)
                {
                    data_t* rptr = rdata;
                    for (uli col=0; col < ncols; ++col, ++rptr, ++pptr)
                    {
                        data_t l = *lptr;
                        data_t r = *rptr;
                        (*pptr) += scale * l * r;
                        if (YetiRuntime::print_cxn)
                        {
                            if ( fabs(*lptr) > 1e-4 && fabs(*rptr) > 1e-4)
                            {
                                std::cout << std::stream_printf(
                                 "%18.12f += %8.4f * %18.12f * %18.12f",
                                *pptr, *lptr, *rptr
                                ) << std::endl;

                            }
                        }
                    }
                }
            }
        }

};

/**
    @class Contraction_tt Contraction. Link indices are rows on the left,
            columns on the right.
            C = A^T * B^T
*/
template <typename data_t>
struct __Contraction_tt {

        static void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (uli row=0; row < nrows; ++row, ++ldata)
            {
                data_t* rtmp = rdata;
                for (uli col=0; col < ncols; ++col, rtmp += nlink, ++pptr)
                {
                    data_t* rptr = rtmp;
                    data_t* lptr = ldata;
                    for (uli link=0; link < nlink; ++link, ++rptr, lptr += nrows)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }


};

/**
    @class Contraction_nn Contraction. Link indices are cols on the left,
            rows on the right.
            C = A * B
*/
template <typename data_t>
struct __Contraction_nn {

        static void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (uli row=0; row < nrows; ++row, ldata += nlink)
            {
                data_t* rtmp = rdata;
                for (uli col=0; col < ncols; ++col, ++rtmp, ++pptr)
                {
                    data_t* rptr = rtmp;
                    data_t* lptr = ldata;
                    for (uli link=0; link < nlink; ++link, rptr += ncols, ++lptr)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }
};

/**
    @class Contraction_nt Contraction. Link indices are both cols.
*/
template <typename data_t>
struct __Contraction_nt {

        static void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (uli row=0; row < nrows; ++row, ldata += nlink)
            {
                data_t* rptr = rdata;
                for (uli col=0; col < ncols; ++col, ++pptr)
                {
                    data_t* lptr = ldata;
                    for (uli link=0; link < nlink; ++link, ++rptr, ++lptr)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }
};

template <typename data_t>
struct Contraction_tt {

    static void contract(
        data_t* ldata,
        data_t* rdata,
        data_t* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        data_t scale
    )
    {
        __Contraction_tt<data_t>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
    }
};

template <typename data_t>
struct Contraction_nt {

    static void contract(
        data_t* ldata,
        data_t* rdata,
        data_t* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        data_t scale
    )
    {
        __Contraction_nt<data_t>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
    }
};

template <typename data_t>
struct Contraction_tn {

    static void contract(
        data_t* ldata,
        data_t* rdata,
        data_t* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        data_t scale
    )
    {
        __Contraction_tn<data_t>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
    }
};

template <typename data_t>
struct Contraction_nn {

    static void contract(
        data_t* ldata,
        data_t* rdata,
        data_t* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        data_t scale
    )
    {
        __Contraction_nn<data_t>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
    }
};


template <>
struct Contraction_nn<double> {

    static void contract(
        double* ldata,
        double* rdata,
        double* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        double scale
    )
    {
#if NO_DGEMM
        if (0)
#else
        if (nlink * nrows * ncols > DGEMM_CUTOFF)
#endif
        {
            const char* opl = "N";
            const char* opr = "N";
            int nrow_ = (int) ncols;
            int ncol_ =  (int) nrows;
            int nlink_ = (int) nlink;
            double beta = 1.0;
            F_DGEMM(opl, opr, &nrow_, &ncol_, &nlink_, &scale, rdata, &nrow_, ldata,
                  &nlink_, &beta, pdata, &nrow_);

        }
        else
        {
           __Contraction_nn<double>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
        }
    }
};




template <>
struct Contraction_nt<double> {

    static void contract(
        double* ldata,
        double* rdata,
        double* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        double scale
    )
    {
#if NO_DGEMM
        if (0)
#else
        if (nlink * nrows * ncols > DGEMM_CUTOFF)
#endif
        {
            const char* opl = "T";
            const char* opr = "N";
            int nrow_ = (int) ncols;
            int ncol_ =  (int) nrows;
            int nlink_ = (int) nlink;
            double beta = 1.0;
            F_DGEMM(opl, opr, &nrow_, &ncol_, &nlink_, &scale, rdata, &nlink_, ldata,
                  &nlink_, &beta, pdata, &nrow_);

        }
        else
        {
           __Contraction_nt<double>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
        }
    }
};

template <>
struct Contraction_tn<double> {

    static void contract(
        double* ldata,
        double* rdata,
        double* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        double scale
    )
    {
#if NO_DGEMM
        if (0)
#else
        if (nlink * nrows * ncols > DGEMM_CUTOFF)
#endif
        {
            const char* opl = "N";
            const char* opr = "T";
            int nrow_ = (int) ncols;
            int ncol_ =  (int) nrows;
            int nlink_ = (int) nlink;
            double beta = 1.0;
            F_DGEMM(opl, opr, &nrow_, &ncol_, &nlink_, &scale, rdata, &nrow_, ldata,
                  &ncol_, &beta, pdata, &nrow_);

        }
        else
        {
           __Contraction_tn<double>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
        }
    }
};

template <>
struct Contraction_tt<double> {

    static void contract(
        double* ldata,
        double* rdata,
        double* pdata,
        uli nrows,
        uli ncols,
        uli nlink,
        double scale
    )
    {
#if NO_DGEMM
        if (0)
#else
        if (nlink * nrows * ncols > DGEMM_CUTOFF)
#endif
        {
            const char* opl = "T";
            const char* opr = "T";
            int nrow_ = (int) ncols;
            int ncol_ =  (int) nrows;
            int nlink_ = (int) nlink;
            double beta = 1.0;
            F_DGEMM(opl, opr, &nrow_, &ncol_, &nlink_, &scale, rdata, &nlink_, ldata,
                  &ncol_, &beta, pdata, &nrow_);

        }
        else
        {
           __Contraction_tt<double>::contract(ldata, rdata, pdata, nrows, ncols, nlink, scale);
        }
    }
};




template <
    template <typename data_t> class cxn_t
>
struct ContractionTemplate :
    public ContractionEngine
{

    public:


        template <class data_t>
        void
        contract(
            char* ldata,
            char* rdata,
            char* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            double scale
        )
        {
            data_t* larr = reinterpret_cast<data_t*>(ldata);
            data_t* rarr = reinterpret_cast<data_t*>(rdata);
            data_t* parr = reinterpret_cast<data_t*>(pdata);
            cxn_t<data_t>::contract(larr, rarr, parr, nrows, ncols, nlink, scale);
        }

        void
        contract(
            DataNode* ldata,
            DataNode* rdata,
            DataNode* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            double scale,
            TemplateInfo::type_t cxn_data_type
        )
        {
            data_type_switch(
                cxn_data_type,
                contract,
                ldata->data(),
                rdata->data(),
                pdata->data(),
                nrows,
                ncols,
                nlink,
                scale
            );
        }

};


}

#ifdef redefine_size_t
#undef size_t
#endif


#endif

