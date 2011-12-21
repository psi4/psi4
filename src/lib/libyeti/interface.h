#ifndef yeti_INTERFACE_H
#define yeti_INTERFACE_H

/*
 *  Defines the interface between YETI and the external quantum chemistry code
 *  that YETI is being linked into.
 */ 
#include <libyeti/class.h>
#include <libyeti/sort.hpp>
#include <libyeti/index.hpp>

#if HAVE_MPQC
#include <libyeti/aobasis.h>
#include <libyeti/mobasis.h>
#include <libyeti/data.h>
#include <libyeti/env.h>
#include <libyeti/exception.h>
#include <libyeti/index.h>
#include <libyeti/matrix.h>
#include <libyeti/permutation.h>
#include <libyeti/permutationimpl.h>
#include <libyeti/sort.h>
#include <libyeti/tensor.h>
#include <libyeti/dataimpl.h>
#include <libyeti/elementop.h>
#include <libyeti/mobasis.h>
#include <libyeti/filler.h>
#include <libyeti/mpqc.hpp>
#include <libyeti/aobasis.hpp>

#include <util/misc/autovec.h>
#include <util/class/scexception.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/basis/gpetite.h>
#include <chemistry/qc/basis/intdescr.h>

#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/zaptr12/rohfwfn.h>
#include <chemistry/qc/zaptr12/zapt2f12.h>
#include <chemistry/qc/zaptr12/matrixkits.h>
#include <chemistry/qc/zaptr12/matrixoperations.h>
#include <chemistry/qc/zaptr12/tbintmerge.h>
#include <chemistry/qc/zaptr12/integralset.h>
#endif

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace sc {

class MPQCMole : 
    public Wavefunction {

    Ref<CLHF> ref_;
    Ref<GaussianBasisSet> obs_;
    int nrepeat_;
    int sleep_;
    int memory_;
    int nthread_;

  public: 
    MPQCMole(const Ref<KeyVal> &);
    MPQCMole(StateIn &);
    void save_data_state(StateOut &);
    void compute(void);
    void obsolete(void);
    int nelectron(void);
    RefSymmSCMatrix density(void);
    int spin_polarized(void);
    int value_implemented(void) const;
    Ref<CLHF> ref(void);

};

Ref<MOIndexSpace> get_ao_space(
    Ref<GaussianBasisSet> obs,
    Ref<Integral> IF,
    const std::string& id
);

}

namespace mpqc {

class MPQC {

    public:
        typedef enum { hcore, overlap } obint_type_t;

        static sc::Ref<sc::GaussianBasisSet> unit_basis;

};

using sc::Ref;
using sc::TwoBodyInt;
using sc::OneBodyInt;
using sc::Integral;
using sc::GaussianBasisSet;
using sc::GaussianShell;
using sc::RefDiagSCMatrix;
using sc::RefSCMatrix;
using sc::RefSymmSCMatrix;
using sc::IntegralTransform;
using sc::TwoBodyIntDescr;
using yeti::TensorElementComputer;
using yeti::TensorValueEstimater;
using yeti::TensorIndexDescr;
using yeti::TensorIndexDescrPtr;
using yeti::IndexDescr;
using yeti::uli;
using yeti::usi;
using yeti::IndexRange;
using yeti::MultiShellMap;
using yeti::MultiShellMapPtr;
using yeti::AOBasisPtr;
using yeti::ScreeningScheme;
using yeti::PartitioningPolicyPtr;
using yeti::TwoElectronEstimableComputer;
using yeti::TEIShellComputeFunctor;

double
min_exponent(const GaussianShell& shell);

void
build_ao_range(
    const Ref<GaussianBasisSet>& bs,
    const std::string& descr,
    uli nfxn_min_per_tile,
    const char* id1,
    const char* id2 = 0,
    const char* id3 = 0,
    const char* id4 = 0,
    PartitioningPolicyPtr policy = 0
);

class MPQCShellComputeFunctor :
    public yeti::TEIShellComputeFunctor
{
    public:
        typedef enum { chemist, physicist } notation_type_t;

    private:

        Ref<TwoBodyInt> tbint_;

        const double* buffer_;

    public:


        MPQCShellComputeFunctor(
            Ref<TwoBodyInt> tbint
        );

        virtual void operator()(uli i, uli j, uli a, uli b) const;

        virtual inline const double* buffer() const { return buffer_; }


};

class TwoElectronIntegralComputer :
    public TwoElectronEstimableComputer

{
    protected:
        bool multishell_;

        const double* buffer_;

        Ref<GaussianBasisSet> bs1_;

        Ref<GaussianBasisSet> bs2_;

        Ref<GaussianBasisSet> bs3_;

        Ref<GaussianBasisSet> bs4_;

        AOBasisPtr aobs1_;

        AOBasisPtr aobs2_;

        AOBasisPtr aobs3_;

        AOBasisPtr aobs4_;

        MultiShellMapPtr shmap1_;

        MultiShellMapPtr shmap2_;

        MultiShellMapPtr shmap3_;

        MultiShellMapPtr shmap4_;

        Ref<TwoBodyInt> tbint_;

        Ref<TwoBodyIntDescr> tbint_descr_;

        Ref<Integral> IF_;

        void init();

        void init_cauchy_schwarz();

        yeti::Sort* sort_;

        uli permuted_indices_[NINDEX];

        uli nelements_[NINDEX];

    public:
        using TensorElementComputer::compute;
        using TwoElectronEstimableComputer::get_estimater;



        TwoElectronIntegralComputer(
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2,
            const Ref<GaussianBasisSet>& bs3,
            const Ref<GaussianBasisSet>& bs4
        );



        TwoElectronIntegralComputer(
            const Ref<TwoBodyIntDescr>& descr,
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet> &bs1,
            const Ref<GaussianBasisSet> &bs2,
            const Ref<GaussianBasisSet> &bs3,
            const Ref<GaussianBasisSet> &bs4
        );

        ~TwoElectronIntegralComputer();

        void compute(const uli* indices, double* data, uli n);

        TensorElementComputer* copy() const;

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

        void sort(yeti::Permutation* p);
};

class TwoElectronIntegralComputer2Index :
    public TensorElementComputer
{
    private:
        bool multishell_;

        const double* buffer_;

        Ref<GaussianBasisSet> bs1_;

        Ref<GaussianBasisSet> bs2_;

        AOBasisPtr aobs1_;

        AOBasisPtr aobs2_;

        MultiShellMapPtr shmap1_;

        MultiShellMapPtr shmap2_;

        Ref<TwoBodyInt> tbint_;

        Ref<TwoBodyIntDescr> tbint_descr_;

        Ref<Integral> IF_;

        void init();

    public:
        using TensorElementComputer::compute;

        TwoElectronIntegralComputer2Index(
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2
         );

        TwoElectronIntegralComputer2Index(
            const Ref<TwoBodyIntDescr>& descr,
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2
         );

        ~TwoElectronIntegralComputer2Index();

        void compute(const uli* indices, double* data, uli n);

        TensorElementComputer* copy() const;

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

};

class TwoElectronIntegralComputer3Index :
    public TensorElementComputer
{
    private:
        bool multishell_;

        const double* buffer_;

        Ref<GaussianBasisSet> bs1_;

        Ref<GaussianBasisSet> bs2_;

        Ref<GaussianBasisSet> bs3_;

        AOBasisPtr aobs1_;

        AOBasisPtr aobs2_;

        AOBasisPtr aobs3_;

        MultiShellMapPtr shmap1_;

        MultiShellMapPtr shmap2_;

        MultiShellMapPtr shmap3_;

        Ref<TwoBodyInt> tbint_;

        Ref<TwoBodyIntDescr> tbint_descr_;

        Ref<Integral> IF_;

        void init();

    public:
        using TensorElementComputer::compute;

        TwoElectronIntegralComputer3Index(
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2,
            const Ref<GaussianBasisSet>& bs3
        );

        TwoElectronIntegralComputer3Index(
            const Ref<TwoBodyIntDescr>& descr,
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2,
            const Ref<GaussianBasisSet>& bs3
        );

        ~TwoElectronIntegralComputer3Index();

        void compute(const uli* indices, double* data, uli n);

        TensorElementComputer* copy() const;

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

};

class OneElectronIntegralComputer :
    public TensorElementComputer
{
    private:
        bool multishell_;

        const double* buffer_;

        Ref<GaussianBasisSet> bs1_;

        Ref<GaussianBasisSet> bs2_;

        AOBasisPtr aobs1_;

        AOBasisPtr aobs2_;

        MultiShellMapPtr shmap1_;

        MultiShellMapPtr shmap2_;

        Ref<OneBodyInt> obint_;

        Ref<Integral> IF_;

        MPQC::obint_type_t obint_type_;

    public:
        using TensorElementComputer::compute;

        OneElectronIntegralComputer(
            MPQC::obint_type_t obint_type,
            const Ref<Integral>& IF,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2
        );

        ~OneElectronIntegralComputer();

        void compute(const uli* indices, double* data, uli n);

        TensorElementComputer* copy() const;

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

};

class DijabElementComputer :
    public TensorElementComputer
{

    private:
        RefDiagSCMatrix ei_;

        RefDiagSCMatrix ej_;

        RefDiagSCMatrix ea_;

        RefDiagSCMatrix eb_;

        uli istart_;

        uli ni_;

        uli jstart_;

        uli nj_;

        uli astart_;

        uli na_;

        uli bstart_;

        uli nb_;

    public:
        using TensorElementComputer::compute;

        DijabElementComputer(
            const RefDiagSCMatrix& ei,
            const RefDiagSCMatrix& ej,
            const RefDiagSCMatrix& ea,
            const RefDiagSCMatrix& eb
        );

        void compute(const uli* indices, double* data, uli n);

        TensorElementComputer* copy() const;

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

};

#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H
class IntegralTransformElementComputer :
    public TensorElementComputer
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

        uli nindex4_;

        TwoBodyInt::tbint_type inttype_;

    public:
        using TensorElementComputer::compute;

        IntegralTransformElementComputer(
            const Ref<IntegralTransform>& tform,
            uli offset_i,
            uli offset_j,
            uli offset_a,
            uli offset_b,
            uli nindex4,
            TwoBodyInt::tbint_type inttype = TwoBodyInt::eri
        );

        TensorElementComputer* copy() const;

        void compute(const uli* indices, double* data, uli n);

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth);

};
#endif //end have zapt library


template <class matrix_t>
class SCMatrixElementComputer :
    public TensorElementComputer
{

    private:
        matrix_t M_;

        int rowstart_;

        int colstart_;

        int nrow_;

        int ncol_;

        uli rowoffset_;

        uli coloffset_;

    public:
        using TensorElementComputer::compute;

        SCMatrixElementComputer(
            const matrix_t& M,
            uli rowoffset = 0,
            uli coloffset = 0
        );

        TensorElementComputer* copy() const;

        void compute(const uli* indices, double* data, uli n);

        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth){
            return yeti::TemplateInfo::double_type;
        }

};

template <class matrix_t>
SCMatrixElementComputer<matrix_t>::SCMatrixElementComputer(
    const matrix_t& M,
    uli rowoffset,
    uli coloffset
) :
    M_(M),
    rowstart_(0),
    colstart_(0),
    nrow_(0),
    ncol_(0),
    rowoffset_(rowoffset),
    coloffset_(coloffset)
{
}

template <class matrix_t>
void
SCMatrixElementComputer<matrix_t>::compute(
    const uli* indices,
    double* data,
    uli n
)
{

    rowstart_ = descr_->get(0)->index_start(indices[0]);
    nrow_ = descr_->get(0)->nelements(indices[0]);
    colstart_ = descr_->get(1)->index_start(indices[1]);
    ncol_ = descr_->get(1)->nelements(indices[1]);

    double* dptr = data;
    for (int r=rowstart_; r < rowstart_ + nrow_; ++r)
    {
        for (int c=colstart_; c < colstart_ + ncol_; ++c, ++dptr)
        {
            double d = M_.get_element(r - rowoffset_, c - coloffset_);
            (*dptr) = d;
        }
    }
}

template <class matrix_t>
TensorElementComputer*
SCMatrixElementComputer<matrix_t>::copy() const
{
    TensorElementComputer* cpy = new SCMatrixElementComputer<matrix_t>(M_);
    cpy->set_index_descr(descr_);
    return cpy;
}

typedef SCMatrixElementComputer<RefSCMatrix> RectMatrixElementComputer;
typedef SCMatrixElementComputer<RefSymmSCMatrix> SymmMatrixElementComputer;

} //end namespace
#elif HAVE_PSI

#include "libmints/typedefs.h"

using namespace yeti;

namespace psi{

    /// The type of notation used in two electron integrals
    enum IntegralNotation {Physicist, Chemist};

    struct _dpdbuf4;
    typedef _dpdbuf4 dpdbuf4;

    double min_exponent(const boost::shared_ptr<GaussianShell>& shell);
    
    void build_ao_range(const boost::shared_ptr<BasisSet>& obs,
                   const char* id1, const char* id2, const char* id3, const char* id4);
    
    class TwoElectronIntegralComputer: public TensorElementComputer {
      protected:
        /// The integral generation object
        boost::shared_ptr<TwoBodyAOInt> eri_;
        /// The integral buffer
        const double *buffer_;
        /// Whether multiple shells are grouped together or not
        bool multishell_;
        /// The map of shells grouped to form index 1
        MultiShellMapPtr shmap1_;
        /// The map of shells grouped to form index 2
        MultiShellMapPtr shmap2_;
        /// The map of shells grouped to form index 3
        MultiShellMapPtr shmap3_;
        /// The map of shells grouped to form index 4
        MultiShellMapPtr shmap4_;
        /// The AO basis object, describing shell groupings for index 1
        AOBasisPtr aobs1_;
        /// The AO basis object, describing shell groupings for index 2
        AOBasisPtr aobs2_;
        /// The AO basis object, describing shell groupings for index 3
        AOBasisPtr aobs3_;
        /// The AO basis object, describing shell groupings for index 4
        AOBasisPtr aobs4_;
      public:
        /**
         * Makes a new two-electron integral computer.
         * @param eri a TwoBodyAOInt object to compute the desired integrals
         */
        TwoElectronIntegralComputer(boost::shared_ptr<TwoBodyAOInt> eri);
    
        /// The copy constructor
        TensorElementComputer* copy() const;
    
        yeti::TemplateInfo::type_t
        element_type(const uli* indices, usi depth) { return yeti::TemplateInfo::double_type; }
    
        /// Performs the integral computation
        void compute(const uli* indices, double* data, uli nelements);
    };
    
    class MatrixFiller : public TensorElementComputer {
        private:
            int callcount_;
            SharedMatrix matrix_;
            yeti::uli poff_;
            yeti::uli qoff_;
        public:
            MatrixFiller(SharedMatrix matrix, yeti::uli poff, yeti::uli qoff);

            TensorElementComputer* copy() const;

            yeti::TemplateInfo::type_t
            element_type(const uli* indices, usi depth) { return yeti::TemplateInfo::double_type; }
    
            void compute(const uli* indices, double* data, uli n);
    };
    
    class VectorFiller : public TensorElementComputer {
        private:
            int callcount_;
            SharedVector vector_;
            yeti::uli offset_;
        public:
            VectorFiller(SharedVector vector, yeti::uli offset);
            TensorElementComputer* copy() const;

            yeti::TemplateInfo::type_t
            element_type(const uli* indices, usi depth) { return yeti::TemplateInfo::double_type; }
    
            void compute(const uli* indices, double* data, uli n);
    };
    
    
    class DPDFiller : public TensorElementComputer {
        private:
            int callcount_;
            dpdbuf4 *buf_;
            yeti::uli poff_;
            yeti::uli qoff_;
            yeti::uli roff_;
            yeti::uli soff_;
        public:
            DPDFiller(dpdbuf4 *buf, yeti::uli poff, yeti::uli qoff, yeti::uli roff, yeti::uli soff);

            TensorElementComputer* copy() const;
    
            yeti::TemplateInfo::type_t
            element_type(const uli* indices, usi depth) { return yeti::TemplateInfo::double_type; }

            void compute(const uli* indices, double* data, yeti::uli n);
    };

} // Namespace psi

#else  // have package
#error You must link to an external package and defining HAVE_PSI, or HAVE_MPQC
#endif //end have package

#ifdef redefine_size_t
#undef size_t
#endif

#endif // INTERFACE_H header guard
