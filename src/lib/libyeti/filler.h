#ifndef yeti_filler_h
#define yeti_filler_h

#include "class.h"

#include "permutation.hpp"
#include "sort.hpp"
#include "index.hpp"
#include "filler.hpp"
#include "data.hpp"
#include "gigmatrix.hpp"

#include "mallocimpl.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

enum ScreeningScheme {
    NoScreening,
    CauchySchwarzScreening
};

class TensorValueEstimater :
    public smartptr::Countable
{

    private:
        static UnitEstimaterPtr unit_;

    public:
        virtual float max_log(const uli* indices) const = 0;

        static UnitEstimater* get_unit_estimater();

};


class MultipleEstimaterMinimumEstimater :
    public TensorValueEstimater
{
    private:

        std::vector<TensorValueEstimaterPtr> estimators_;

    public:

        MultipleEstimaterMinimumEstimater();

        ~MultipleEstimaterMinimumEstimater();

        void add_estimator(TensorValueEstimaterPtr estimator);

        virtual float max_log(const uli *indices) const;
};


class UnitEstimater :
    public TensorValueEstimater
{

    public:
        float max_log(const uli* indices) const;

};

class TensorElementFilter :
    public smartptr::Countable
{
    protected:
        TensorIndexDescr* descr_;

    public:
        TensorElementFilter();

        void set_index_descr(TensorIndexDescr* descr);

        virtual bool is_zero(const uli* indices) const = 0;
};

class TensorSymmetryFilter :
    public TensorElementFilter
{

    public:
        bool is_zero(const uli *indices) const;

};

class TensorElementComputer :
    public smartptr::Countable
{

    private:
        void* buffer_;

    protected:
        TensorIndexDescrPtr descr_;

    public:
        TensorElementComputer();

        virtual ~TensorElementComputer();

        /**
            @return A copy of the element computer.  This must return a "thread-safe"
            copy for which calls to compute do not interfere between parent and copy.
        */
        virtual TensorElementComputer* copy() const = 0;

        virtual TemplateInfo::type_t element_type(const uli* indices, usi depth) = 0;

        virtual void compute(const uli* indices, double* data, uli n);

        virtual void compute(const uli* indices, quad* data, uli n);

        virtual void compute(const uli* indices, int* data, uli n);

        virtual void compute(const uli* indices, float* data, uli n);

        virtual TensorValueEstimater* get_estimater(usi depth) const;

        void allocate_buffer(uli maxblocksize);

        void copy_index_descr(TensorIndexDescr* descr);

        void set_index_descr(const TensorIndexDescrPtr& descr);

        virtual void sort(Permutation* p);

};

/**
  Abstract superclass for all computer objects that involve two electron integrals.
  Mostly allows for the standardization across backends of TensorElementEstimater
  assignment.
  */
class TwoElectronEstimableComputer
    : public TensorElementComputer
{
    protected:

        std::vector<TensorValueEstimater*> estimators_;

        TwoElectronEstimableComputer();

        virtual void init_cauchy_schwarz();

    public:

        virtual TensorValueEstimater* get_estimater(usi depth) const;

        virtual void set_screening_scheme(ScreeningScheme scheme);

};

class TEIShellComputeFunctor :
    public smartptr::Countable
{
    public:

        virtual const double* buffer() const = 0;

        virtual void operator()(uli, uli, uli, uli) const = 0;
};

class MemsetElementComputer :
    public TensorElementComputer
{

    private:
        template <typename data_t>
        void
        memset(uli n, data_t* data);

        TemplateInfo::type_t element_type_;

    public:
        MemsetElementComputer(TemplateInfo::type_t elem_type);

        void compute(const uli* indices, double* data, uli n);

        void compute(const uli* indices, quad* data, uli n);

        void compute(const uli* indices, int* data, uli n);

        void compute(const uli* indices, float* data, uli n);

        TensorElementComputer* copy() const;

        TemplateInfo::type_t element_type(const uli *indices, usi depth);

};

class ThreadedTensorElementComputer :
    public smartptr::Countable
{

    protected:
        std::vector<TensorElementComputerPtr> fillers_;

        usi mindepth_;

    public:
        ThreadedTensorElementComputer(const TensorElementComputerPtr& comp);

        virtual ~ThreadedTensorElementComputer();

        template <typename data_t>
        void compute(
            const uli* indices,
            data_t* data,
            uli n,
            uli threadnum
        )
        {
            fillers_[threadnum]->compute(indices, data, n);
        }

        void allocate_buffer(uli maxblocksize);

        void operator=(const TensorElementComputerPtr& comp);

        TensorElementComputer* get_computer(uli threadnum) const;

        usi mindepth() const;

        void set_mindepth(usi depth);

        void sort(Permutation* p);

};

class DoubleArrayElementComputer :
    public TensorElementComputer
{
    private:
        const double* data_;

        const double* dataptr_;

    public:
        using TensorElementComputer::compute;

        DoubleArrayElementComputer(const double* data);

        TensorElementComputer* copy() const;

        void compute(const uli* indices, double* data, uli n);

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth){
            return yeti::TemplateInfo::double_type;
        }
};


class YetiMatrixElementComputer :
    public TensorElementComputer
{
    private:
        MatrixPtr M_;

        int rowstart_;

        int colstart_;

        int nrow_;

        int ncol_;

    public:
        using TensorElementComputer::compute;

        YetiMatrixElementComputer(const MatrixPtr& m);

        TensorElementComputer* copy() const;

        void compute(const uli* indices, double* data, uli n);

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth){
            return yeti::TemplateInfo::double_type;
        }

};


class DiagonalMatrixEstimater :
    public TensorValueEstimater
{

    public:
        float max_log(const uli* indices) const;

};

class DiagonalMatrixElementComputer :
    public TensorElementComputer
{

    private:
        double* values_;

        uli nvalues_;

        TensorValueEstimaterPtr diagonal_matrix_estimater_;

    public:
        using TensorElementComputer::compute;

        DiagonalMatrixElementComputer(
            const double* vals,
            uli nvals
        );

        ~DiagonalMatrixElementComputer();

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth){
            return yeti::TemplateInfo::double_type;
        }

        virtual TensorValueEstimater* get_estimater(usi depth) const;

        TensorElementComputer* copy() const;

        void compute(const uli* indices, double* data, uli n);
};


class IdentityMatrixElementComputer :
    public DiagonalMatrixElementComputer
{

    public:
        IdentityMatrixElementComputer(uli nvalues);

};

class Diagonal_IJIJ_ValueEstimater :
    public TensorValueEstimater
{
    private:
        float iiii_bound_;

        float ijij_bound_;

        float ijji_bound_;

    public:
        Diagonal_IJIJ_ValueEstimater(
            double iiii,
            double ijij,
            double ijji
        );

    float max_log(const uli* indices) const;
};

class DenominatorElementComputer :
    public TensorElementComputer
{

    private:
        double* evals_;

    public:
       using TensorElementComputer::compute;

       DenominatorElementComputer(
           double* evals
       );

       ~DenominatorElementComputer();

       yeti::TemplateInfo::type_t
       element_type(const uli* indices, usi depth){
           return yeti::TemplateInfo::double_type;
       }

       void compute(const uli* indices, double* data, uli n);

       TensorElementComputer* copy() const;
};

class Diagonal_IJIJ_ElementComputer :
    public TensorElementComputer
{

    private:
        double ijij_;

        double ijji_;

        double iiii_;

        Diagonal_IJIJ_ValueEstimater* estimater_;

    public:
       using TensorElementComputer::compute;

       Diagonal_IJIJ_ElementComputer(
           double iiii,
           double ijij,
           double ijji
       );

       ~Diagonal_IJIJ_ElementComputer();

       yeti::TemplateInfo::type_t
       element_type(const uli* indices, usi depth){
           return yeti::TemplateInfo::double_type;
       }

       virtual TensorValueEstimater* get_estimater(usi depth) const;

       void compute(const uli* indices, double* data, uli n);

       TensorElementComputer* copy() const;
};


class UnitElementComputer :
    public TensorElementComputer
{

    public:
       using TensorElementComputer::compute;

       yeti::TemplateInfo::type_t
       element_type(const uli* indices, usi depth){
           return yeti::TemplateInfo::double_type;
       }

       void compute(const uli* indices, double* data, uli n);

       TensorElementComputer* copy() const;
};




} //end namespace yeti

#ifdef redefine_size_t
#undef size_t
#endif

#endif

