#ifndef yeti_MPQC_H
#define yeti_MPQC_H



#include "tile.h"
#include "index.h"
#include "class.h"

#if HAVE_MPQC
#include <util/misc/autovec.h>
#include <util/class/scexception.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/basis/gpetite.h>
#include <chemistry/qc/basis/intdescr.h>

#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H
#include <chemistry/qc/zaptr12/rohfwfn.h>
#include <chemistry/qc/zaptr12/zapt2f12.h>
#include <chemistry/qc/zaptr12/matrixkits.h>
#include <chemistry/qc/zaptr12/matrixoperations.h>
#include <chemistry/qc/zaptr12/tbintmerge.h>
#include <chemistry/qc/zaptr12/integralset.h>
#endif


namespace mpqc {

using sc::Ref;
using sc::TwoBodyInt;
using sc::Integral;
using sc::GaussianBasisSet;
using sc::RefDiagSCMatrix;
using sc::RefSCMatrix;
using sc::RefSymmSCMatrix;
using sc::IntegralTransform;
using sc::TwoBodyIntDescr;
using yeti::Tile;
using yeti::TileElementComputer;
using yeti::IndexRangeTuplePtr;
using yeti::uli;
using yeti::IndexRange;

IndexRange*
get_ao_range(const Ref<GaussianBasisSet>& obs);

class TwoElectronIntegralComputer :
    public TileElementComputer
{
    private:
        const double* buffer_;

        Ref<GaussianBasisSet> bs1_;

        Ref<GaussianBasisSet> bs2_;

        Ref<GaussianBasisSet> bs3_;

        Ref<GaussianBasisSet> bs4_;

        Ref<TwoBodyInt> tbint_;

        Ref<TwoBodyIntDescr> descr_;

        Ref<Integral> IF_;

        uli nblock_;

        void init();

    public:
        TwoElectronIntegralComputer(
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2 = 0,
            const Ref<GaussianBasisSet>& bs3 = 0,
            const Ref<GaussianBasisSet>& bs4 = 0
        );

        TwoElectronIntegralComputer(
            const Ref<TwoBodyIntDescr>& descr,
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet> &bs1,
            const Ref<GaussianBasisSet> &bs2 = 0,
            const Ref<GaussianBasisSet> &bs3 = 0,
            const Ref<GaussianBasisSet> &bs4 = 0
       );

        ~TwoElectronIntegralComputer();

        void compute(Tile* tile, double* data);

        TileElementComputer* copy() const;

};

class DijabElementComputer :
    public TileElementComputer
{

    private:
        RefDiagSCMatrix ei_;

        RefDiagSCMatrix ej_;

        RefDiagSCMatrix ea_;

        RefDiagSCMatrix eb_;

        size_t istart_;

        size_t ni_;

        size_t jstart_;

        size_t nj_;

        size_t astart_;

        size_t na_;

        size_t bstart_;

        size_t nb_;

    public:
        DijabElementComputer(
            const RefDiagSCMatrix& ei,
            const RefDiagSCMatrix& ej,
            const RefDiagSCMatrix& ea,
            const RefDiagSCMatrix& eb
        );

        void compute(Tile* tile, double* data);

        TileElementComputer* copy() const;

};

#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H
class IntegralTransformElementComputer :
    public TileElementComputer
{

    private:
        Ref<IntegralTransform> tform_;

        uli offset_i_;

        uli offset_j_;

        uli offset_a_;

        uli offset_b_;

        uli istart_;

        uli ni_;

        uli jstart_;

        uli nj_;

        uli astart_;

        uli na_;

        uli bstart_;

        uli nb_;

        size_t nindex4_;

    public:
        IntegralTransformElementComputer(
            const Ref<IntegralTransform>& tform,
            uli offset_i,
            uli offset_j,
            uli offset_a,
            uli offset_b,
            size_t nindex4
        );

        TileElementComputer* copy() const;

        void compute(Tile* tile, double* data);


};
#endif


template <class matrix_t>
class SCMatrixElementComputer :
    public TileElementComputer
{

    private:
        matrix_t M_;

        int rowstart_;

        int colstart_;

        int nrow_;

        int ncol_;

    public:
        SCMatrixElementComputer(
            const matrix_t& M
        );

        TileElementComputer* copy() const;

        void compute(Tile* tile, double* data);

};

template <class matrix_t>
SCMatrixElementComputer<matrix_t>::SCMatrixElementComputer(
   const matrix_t& M
) :
    M_(M),
    rowstart_(0),
    colstart_(0),
    nrow_(0),
    ncol_(0)
{
}

template <class matrix_t>
void
SCMatrixElementComputer<matrix_t>::compute(Tile* tile, double* data)
{
    IndexRangeTuplePtr tuple(tile->get_index_ranges());
    rowstart_ = tuple->get(0)->start();
    nrow_ = tuple->get(0)->n();
    colstart_ = tuple->get(1)->start();
    ncol_ = tuple->get(1)->n();

    double* dptr = data;
    for (int r=rowstart_; r < rowstart_ + nrow_; ++r)
    {
        for (int c=colstart_; c < colstart_ + ncol_; ++c, ++dptr)
        {
            (*dptr) = M_.get_element(r, c);
        }
    }
}

template <class matrix_t>
TileElementComputer*
SCMatrixElementComputer<matrix_t>::copy() const
{
    return new SCMatrixElementComputer<matrix_t>(M_);
}

typedef SCMatrixElementComputer<RefSCMatrix> RectMatrixElementComputer;
typedef SCMatrixElementComputer<RefSymmSCMatrix> SymmMatrixElementComputer;

} //end namespace

#endif

#endif // MPQC_H header guard
