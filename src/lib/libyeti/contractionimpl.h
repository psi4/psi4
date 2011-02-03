#ifndef yeti_contraction_impl_h
#define yeti_contraction_impl_h

#include "data.h"
#include "contraction.h"

#include <libsmartptr/printstream.h>

namespace yeti {


template <class data_t>
void
transpose_print(
        data_t* data,
        size_t nrows,
        size_t ncols
)
{
    std::cout << "Matrix " << ncols << " x " << nrows << std::endl;
    for (size_t c=0; c < ncols; ++c, ++data)
    {
        data_t* dptr = data;
        for (size_t r=0; r < nrows; ++r, dptr += ncols)
        {
            std::cout << " " << stream_printf(TypeInfo<data_t>::printf_str, *dptr);
        }
        std::cout << std::endl;
    }
}

template <class data_t>
void
normal_print(
        data_t* data,
        size_t nrows,
        size_t ncols
)
{
    std::cout << "Matrix " << nrows << " x " << ncols << std::endl;
    for (size_t r=0; r < nrows; ++r)
    {
        for (size_t c=0; c < ncols; ++c, ++data)
        {
            std::cout << " " << stream_printf(TypeInfo<data_t>::printf_str, *data);
        }
        std::cout << std::endl;
    }
}

/**
    @class Contraction_tt Contraction. Both link indices are row indices.
*/
struct Contraction_tt {

        template <
            class data_t
        >
        void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            data_t scale
        )
        {
            //now accumulate the contraction
            for (size_t link=0; link < nlink; ++link, ldata += nrows, rdata += ncols)
            {
                data_t* lptr = ldata;
                data_t* pptr = pdata;
                for (size_t row=0; row < nrows; ++row, ++lptr)
                {
                    data_t* rptr = rdata;
                    for (size_t col=0; col < ncols; ++col, ++rptr, ++pptr)
                    {
                        data_t l = *lptr;
                        data_t r = *rptr;
                        //(*pptr) += nisotropy * (*lptr) * (*rptr);
                        (*pptr) += scale * l * r;
                    }
                }
            }
        }

        template <
            class data_t
        >
        void print_left_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Left Transpose" << std::endl;
            transpose_print<data_t>(data, nrows, ncols);
        }

        template <
            class data_t
        >
        void print_right_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Right Transpose" << std::endl;
            transpose_print<data_t>(data, nrows, ncols);
        }


};

struct Contraction_tn {

        template <
            class data_t
        >
        void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (size_t row=0; row < nrows; ++row, ++ldata)
            {
                data_t* rtmp = rdata;
                for (size_t col=0; col < ncols; ++col, rtmp += nlink, ++pptr)
                {
                    data_t* rptr = rtmp;
                    data_t* lptr = ldata;
                    for (size_t link=0; link < nlink; ++link, ++rptr, lptr += nrows)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }

        template <
            class data_t
        >
        void print_left_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Left Transpose" << std::endl;
            transpose_print<data_t>(data, nrows, ncols);
        }

        template <
            class data_t
        >
        void print_right_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Right Normal" << std::endl;
            normal_print<data_t>(data, nrows, ncols);
        }


};

struct Contraction_nt {

        template <
            class data_t
        >
        void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (size_t row=0; row < nrows; ++row, ldata += nlink)
            {
                data_t* rtmp = rdata;
                for (size_t col=0; col < ncols; ++col, ++rtmp, ++pptr)
                {
                    data_t* rptr = rtmp;
                    data_t* lptr = ldata;
                    for (size_t link=0; link < nlink; ++link, rptr += ncols, ++lptr)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }

        template <
            class data_t
        >
        void print_left_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Left Normal" << std::endl;
            normal_print<data_t>(data, nrows, ncols);
        }

        template <
            class data_t
        >
        void print_right_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Right Transpose" << std::endl;
            transpose_print<data_t>(data, nrows, ncols);
        }

};

struct Contraction_nn {

        template <
            class data_t
        >
        void contract(
            data_t* ldata,
            data_t* rdata,
            data_t* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            data_t scale
        )
        {
            data_t* pptr = pdata;
            for (size_t row=0; row < nrows; ++row, ldata += nlink)
            {
                data_t* rptr = rdata;
                for (size_t col=0; col < ncols; ++col, ++pptr)
                {
                    data_t* lptr = ldata;
                    for (size_t link=0; link < nlink; ++link, ++rptr, ++lptr)
                    {
                        (*pptr) += scale * (*lptr) * (*rptr);
                    }
                }
            }
        }

        template <
            class data_t
        >
        void print_left_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Left Normal" << std::endl;
            normal_print<data_t>(data, nrows, ncols);
        }

        template <
            class data_t
        >
        void print_right_matrix(
            data_t* data,
            size_t nrows,
            size_t ncols
        )
        {
            std::cout << "Right Normal" << std::endl;
            normal_print<data_t>(data, nrows, ncols);
        }

};

struct ContractionEngine {

    public:
        virtual void contract(
            Data* ldata,
            Data* rdata,
            Data* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            void* ltmp_void,
            void* rtmp_void,
            void* ptmp_void,
            double scale,
            TemplateInfo::type_t dtype
        ) = 0;

        virtual void print_left_matrix(
            Data* data,
            size_t nrows,
            size_t ncols,
            void* tmp
        ) = 0;

        virtual void print_right_matrix(
            Data* data,
            size_t nrows,
            size_t ncols,
            void* tmp
        ) = 0;

};

template <
    class cxn_t
>
struct ContractionTemplate : public ContractionEngine {

    private:
        cxn_t cxn_;

    public:

        template <class data_t>
        void
        contract(
            Data* ldata,
            Data* rdata,
            Data* pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            void* ltmp_void,
            void* rtmp_void,
            void* ptmp_void,
            double scale
        )
        {
            data_t* larr = ldata->cast(reinterpret_cast<data_t*>(ltmp_void));
            data_t* rarr = rdata->cast(reinterpret_cast<data_t*>(rtmp_void));
            data_t* parr; pdata->get(&parr);
            cxn_.contract<data_t>(larr, rarr, parr, nrows, ncols, nlink, scale);
        }

        void
        contract(
            Data *ldata,
            Data *rdata,
            Data *pdata,
            size_t nrows,
            size_t ncols,
            size_t nlink,
            void *ltmp_void,
            void *rtmp_void,
            void *ptmp_void,
            double scale,
            TemplateInfo::type_t cxn_data_type
        )
        {
            data_type_switch(
                cxn_data_type,
                contract,
                ldata,
                rdata,
                pdata,
                nrows,
                ncols,
                nlink,
                ltmp_void,
                rtmp_void,
                ptmp_void,
                scale
            );
        }

        template <typename data_t>
        void
        print_left_matrix_tmpl(
            Data *data,
            size_t nrows,
            size_t ncols,
            void *tmp
        )
        {
            data_t* arr = data->cast(reinterpret_cast<data_t*>(tmp));
            cxn_.template print_left_matrix<data_t>(arr, nrows, ncols);
        }

        template <typename data_t>
        void
        print_right_matrix_tmpl(
            Data *data,
            size_t nrows,
            size_t ncols,
            void *tmp
        )
        {
            data_t* arr = data->cast(reinterpret_cast<data_t*>(tmp));
            cxn_.template print_right_matrix<data_t>(arr, nrows, ncols);
        }

        void print_left_matrix(
            Data* data,
            size_t nrows,
            size_t ncols,
            void* tmp
        )
        {
            data_type_switch(data->type(), print_left_matrix_tmpl, data, nrows, ncols, tmp);
        }



        void print_right_matrix(
            Data* data,
            size_t nrows,
            size_t ncols,
            void* tmp
        )
        {
            data_type_switch(data->type(), print_right_matrix_tmpl, data, nrows, ncols, tmp);
        }

};


}




#endif

