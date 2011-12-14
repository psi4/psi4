#ifndef yeti_ELEMENTOP_H
#define yeti_ELEMENTOP_H

#include "class.h"
#include "mapimpl.h"

#include "data.hpp"
#include "index.hpp"
#include "elementop.hpp"
#include "tensorblock.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class ElementOp :
    public smartptr::Countable
{

    public:
        virtual ~ElementOp();

        virtual void configure(Tensor* tensor);

        // The element_op methods are _basically_ the
        // equivalent of pure virtual methods that
        // throw an error at runtime instead of compile
        // time, allowing subclasses to implement some but
        // not all of the element_op varieties

        virtual void element_op(
            uli nblock,
            const uli* indices,
            const uli* sizes,
            double* data
        );

        virtual void element_op(
            uli nblock,
            const uli* indices,
            const uli* sizes,
            int* data
        );

        virtual void element_op(
            uli nblock,
            const uli* indices,
            const uli* sizes,
            float* data
        );

        virtual void element_op(
            uli nblock,
            const uli* indices,
            const uli* sizes,
            quad* data
        );

        virtual void retrieve(TensorBlock* block) const = 0;

        virtual void release(TensorBlock* block) const = 0;

        virtual bool do_update_after() const;


};

class DiaOp :
    public ElementOp
{

    private:
        const double* evals_i_;
        const double* evals_a_;
        uli istart_;
        uli istop_;
        uli astart_;
        uli astop_;

    public:
        using ElementOp::element_op;

        DiaOp(
            const double* ei,
            const double* ea
        );

        /**
          Carries out the operation characteristic of the DiaOp class.
          Specifically, it multiplies each element of the block from data by 1/(ei - ea)
          (Note that this should never be called at the user level).

          @param nblock The number of items in the being operated on. (Unused in current implementation).
          @param index_starts Pointer to a two element array containing the starting indices for i and a, in that order.
          @param sizes Pointer to a two element array containing the number of i and a elements in the block to be operated on.
          @param data Pointer to the data block being operated on.

          */
        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;


};

class DijabOp :
    public ElementOp
{

    private:
        const double* evals_i_;
        const double* evals_j_;
        const double* evals_a_;
        const double* evals_b_;
        double ci_;
        double cj_;
        double ca_;
        double cb_;

    public:
        using ElementOp::element_op;

        DijabOp(
            const double* ei,
            const double* ej,
            const double* ea,
            const double* eb,
            double ci = 1.0,
            double cj = 1.0,
            double ca = -1.0,
            double cb = -1.0
        );

        /**
          Carries out the operation characteristic of the DiaOp class.
          Specifically, it multiplies each element of the block from data by
          1/(ci*ei + cj*ej + ca*ea + cb*eb).
          (Note that this should never be called at the user level).

          @param nblock The number of items in the block being operated on. (Unused in current implementation).
          @param index_starts Pointer to a four element array containing the starting indices for i, j, a, and b in that order.
          @param sizes Pointer to a four element array containing the number of i, j, a, and b elements in the block to be operated on.
          @param data Pointer to the data block being operated on.

          */
        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;


};

class ZeroOp :
    public ElementOp
{
    private:
        template <typename data_t>
        void
        _zero(
            uli nblock,
            data_t* data
        );

    public:
        using ElementOp::element_op;

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            int* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            float* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            quad* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;

        bool do_update_after() const;

};

class ScaleOp :
    public ElementOp
{
    private:
        double scale_;

        template <typename data_t>
        void
        _scale(
            uli nblock,
            data_t* data
        );

    public:
        using ElementOp::element_op;

        ScaleOp(double scale);

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            int* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            float* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            quad* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;


};

class NormElementOp :
    public ElementOp
{
    private:
        double normsq_;


    public:
        using ElementOp::element_op;

        NormElementOp();

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            int* data
        );

        double norm() const;

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;
};

class Diagonalize_IJIJ_Op :
    public ElementOp
{

    public:
        using ElementOp::element_op;

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;
};

class UpperTriangleGetOp :
    public ElementOp
{
    
    private:
        double* utri_;

        uli nindex_;

        uli offset_p_;

        uli offset_q_;

    public:
        using ElementOp::element_op;

        UpperTriangleGetOp();

        ~UpperTriangleGetOp();

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;

        void configure(Tensor* tensor);

        const double* data() const;

        bool do_update_after() const {return false;}
};

class DoubleArrayGetOp :
    public ElementOp
{
    private:
        double* data_;

        double* dataptr_;

        TensorIndexDescr* descr_;

        uli offset_p_;

        uli offset_q_;

        uli offset_r_;

        uli np_;

        uli nq_;

        uli nr_;
    
    public:
        using ElementOp::element_op;
        
        DoubleArrayGetOp();

        ~DoubleArrayGetOp();

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void increment_ptr(uli nelements);

        void configure(Tensor* tensor);

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;

        const double* data() const;

        bool do_update_after() const {return false;}
};


} //namespace yeti

#ifdef redefine_size_t
#undef size_t
#endif

#endif // ELEMENTOP_H
