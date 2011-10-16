#include "tensorparser.h"
#include "elementop.h"
#include "tensoraction.h"
#include "index.h"
#include "runtime.h"
#include "permutation.h"
#include "tensor.h"
#include "node.h"
#include "tensorblock.h"
#include "malloc.h"
#include "matrix.h"
#include "exception.h"
#include "dataimpl.h"
#include "filler.h"
#include "contraction.h"

#include "blas.h"
extern "C" {

extern void dspsvx(const char* fact, const char* uplo, const int* n, const int* nrhs,
                       const double* AP, double* AFP, int* ipiv, const double* BB, const int* nb,
                       double* XX, const int* nx, double* rcond, double* ferr, double* berr,
                       double* work, int* iwork, int* info);

extern void DSPTRS(const char* uplo, const int* n, const int* nrhs,
                       const double* AP, int* ipiv, double* BB, const int* nb,
                       int* info);

extern void dsptrs_(const char* uplo, const int* n, const int* nrhs,
                       const double* AP, int* ipiv, double* BB, const int* nb,
                       int* info);

extern void dsptrs(const char* uplo, const int* n, const int* nrhs,
                       const double* AP, int* ipiv, double* BB, const int* nb,
                       int* info);
}//EndExternC

#include <vector>

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define ALLOW_PERMUTATIONS 1

#define DEBUG_CXN_SORT 0

static std::list<PermutationRuntimeParserPtr> cxn_symmetrization_list;

StringParser::StringParser(const std::string &str)
    : str_(str)
{
}

const std::string&
StringParser::str() const
{
    return str_;
}

uli
StringParser::nentries() const
{
    return entries_.size();
}

void
StringParser::tokenize(const char *token)
{
    entries_.clear();
    if (str_.size()) //only create entries if a string exists
        smartptr::split(entries_, str_, token);
}

StringParser::iterator
StringParser::begin() const
{
    return entries_.begin();
}

StringParser::iterator
StringParser::end() const
{
    return entries_.end();
}

void
StringParser::count(std::map<std::string, uli> &counts) const
{
    iterator it(begin());
    iterator stop(end());
    for ( ; it != stop; ++it)
        ++counts[*it];
}

bool
StringParser::contains(const std::string &str) const
{
    iterator it(begin());
    iterator stop(end());
    for ( ; it != stop; ++it)
    {
        if (str == *it)
            return true;
    }
    return false;
}

uli
StringParser::index(const std::string &str) const
{
    iterator it(begin());
    iterator stop(end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        if (*it == str)
            return idx;
    }

    return nentries();
}

const std::vector<std::string>&
StringParser::entries() const
{
    return entries_;
}

void
StringParser::str(const std::string &str)
{
    entries_.clear();
    str_ = str;
}

const std::string&
StringParser::substr(uli index) const
{
    return entries_.at(index);
}

Permutation*
StringParser::get_permutation(const PermutationRuntimeParserPtr& perm)
{
    usi indexmap[NINDEX];
    usi nindex = nentries();
    for (usi i=0; i < nindex; ++i)
        indexmap[i] = i;

    usi nentries = this->nentries();

    if (perm->nentries() == 2) //transposition
    {
        usi i = index(perm->substr(0));
        usi j = index(perm->substr(1));
        indexmap[i] = j;
        indexmap[j] = i;

        if (i >= nentries || j >= nentries)
        {
            return 0;
        }
    }
    else
    {
        //vector<std::string> indices;
        vector<std::string>::const_iterator it(perm->begin());
        vector<std::string>::const_iterator stop(perm->end());
        usi idx = 0;
        for ( ; it != stop; ++it, ++idx)
        {
            const std::string& str = *it;
            usi i = index(str);
            indexmap[idx] = i;

            if (i >= nentries)
            {
                return 0;
            }
        }

    }

    Permutation* p = Permutation::get_permutation(perm->sign(), nindex, indexmap);
    return p;
}

usi
StringParser::get_nindex_front(const std::string &str)
{
    StringParser parser(str);
    parser.tokenize(", ");
    usi nindex = parser.nentries();

    for (usi i=0; i < parser.nentries(); ++i)
    {
        usi idxcheck = this->index(parser.substr(i));
        if (idxcheck != i)
        {
            cerr << str << "does not specify front indices" << endl;
            abort();
        }
    }

    return nindex;
}

YetiTensor::YetiTensor(
    const std::string& name,
    const std::string& str,
    const PermutationRuntimeParserPtr& p1,
    const PermutationRuntimeParserPtr& p2,
    const PermutationRuntimeParserPtr& p3,
    const PermutationRuntimeParserPtr& p4,
    const PermutationRuntimeParserPtr& p5
)
    : StringParser(str),
    tensor_(0),
    scale_(1.0),
    perm_(0)
{
    init(name, str, p1, p2, p3, p4, p5);
}

YetiTensor::YetiTensor()
    : StringParser(""),
    tensor_(0),
    scale_(1.0),
    perm_(0)
{
}

YetiTensor::~YetiTensor()
{
    tensor_ = 0;
}

void
YetiTensor::init(
    const std::string& name,
    const std::string& a1,
    const std::string& a2,
    const std::string& a3,
    const std::string& a4
)
{
    std::string str = a1 + "," + a2 + "," + a3 + "," + a4;
    init(name, str);
}

void
YetiTensor::add(
    const PermutationRuntimeParserPtr& perm
)
{
    Permutation* p = StringParser::get_permutation(perm);
    if (!p)
    {
        cerr << "Invalid permutation parser " << perm->str() << " for tensor " << str() << endl;
        abort();
    }
    tensor_->get_tensor_grp()->add(p);
}

void
YetiTensor::init(
    const std::string& name,
    const std::string& a1,
    const std::string& a2,
    const std::string& a3
)
{
    std::string str = a1 + "," + a2 + "," + a3;
    init(name, str);
}

void
YetiTensor::init(
    const std::string& name,
    const std::string& str,
    const PermutationRuntimeParserPtr& p1,
    const PermutationRuntimeParserPtr& p2,
    const PermutationRuntimeParserPtr& p3,
    const PermutationRuntimeParserPtr& p4,
    const PermutationRuntimeParserPtr& p5
)
{
    str_ = str;

    tokenize(", ");
    usi nindex = entries_.size();
    TensorIndexDescrPtr indexdescr = new TensorIndexDescr(nindex);
    for (usi i=0; i < nindex; ++i)
    {
        IndexDescrPtr descr = YetiRuntime::get_descr(entries_[i]);
        indexdescr->set(i, descr);
    }

    PermutationGroupPtr grp = new PermutationGroup(nindex);

    tensor_ = new Tensor(
                    name,
                    indexdescr,
                    grp
                  );
    perm_ = grp->get_identity();

    scale_ = 1.0;

#if ALLOW_PERMUTATIONS
    incref();
    if (p1) add(p1);
    if (p2) add(p2);
    if (p3) add(p3);
    if (p4) add(p4);
    if (p5) add(p5);
    decref();
    grp->close();
    tensor_->recompute_permutation();
#endif

#if 0
    ulli totalsize = tensor_->get_totalsize();
    if (totalsize != 1)
    {
        size_t billion = 1e9;
        size_t million = 1e6;
        size_t thousand = 1e3;
        size_t nbillions = totalsize / billion;
        size_t remainder = totalsize - billion * nbillions;
        size_t nmillions = remainder / million;
        remainder -= million * nmillions;
        size_t nthousands = remainder / thousand;
        remainder -= thousand * nthousands;


        dout << stream_printf("Initialized tensor %30s of total size ", name.c_str());
        if (nbillions)
            dout << nbillions << ",";
        if (nmillions)
            dout << nmillions << ",";
        if (nthousands)
            dout << nthousands << ",";
        dout << remainder << endl;
    }
#endif
}

YetiTensor::YetiTensor(
    const TensorPtr& tensor,
    Permutation* perm,
    double scale,
    const std::string& str
) : tensor_(tensor),
    perm_(perm),
    scale_(scale),
    StringParser(str)
{
    tokenize(", ");
    bool subtensor = false;
    usi nindex = tensor->get_block_descr()->nindex();
    TensorIndexDescrPtr indexdescr = new TensorIndexDescr(nindex);

    if (nentries() != nindex)
    {
        cerr << "Incorrect number of indices for tensor "
             << tensor->get_name() << " with string " << str << endl;
        abort();
    }

    for (usi i=0; i < nindex; ++i)
    {
        const std::string& label = substr(i);
        IndexDescr* descr = YetiRuntime::get_descr(label);
        IndexDescr* mydescr = tensor->get_block_descr()->get(i);
        if (mydescr != descr)
        {
            subtensor = true;
            //validate the subtensor
            if (!YetiRuntime::is_valid_subrange(
                mydescr->get_parent_range(), 
                descr->get_parent_range()
            ))
            {
                cerr << "Invalid index " << label << " at index " << i
                    << " passed to tensor " << tensor->get_name()
                     << endl;
                cerr << "This is probably a type-o in the code.  Alternatively, you "
                        "may have sorted a subtensor and forgotten to sort it back."
                        << endl;
                abort();
            }
        }
        indexdescr->set(i, descr);
    }

    if (subtensor)
    {
        //create a subtensor object    
        tensor_ = new Tensor(
                    indexdescr,
                    tensor.get()
                  );
    }
}

YetiTensor::YetiTensor(
    const std::string &name,
    const std::string &idxA,
    const std::string &idxB,
    const std::string &idxC,
    const std::string &idxD
) :
    StringParser(""),
    scale_(1.0)
{
    std::string str = idxA + "," + idxB + "," + idxC + "," + idxD;
    this->init(name, str);
}

YetiTensor::YetiTensor(
    const std::string &name,
    const std::string &idxA,
    const std::string &idxB,
    const std::string &idxC
) :
    StringParser(""),
    scale_(1.0)
{
    std::string str = idxA + "," + idxB + "," + idxC;
    this->init(name, str);
}

YetiTensor::YetiTensor(const YetiTensorPtr& yt)
    : StringParser(yt->str()),
      tensor_(yt->get_tensor()),
      scale_(yt->get_scale()),
      perm_(yt->get_permutation())
{
}

YetiTensorPtr
YetiTensor::copy() const
{
    TensorPtr cpy = tensor_->copy();
    YetiTensorPtr yt = new YetiTensor(cpy, perm_, scale_, str_);
    return yt;
}

Tensor*
YetiTensor::operator ->() const
{
    return tensor_.get();
}

void
YetiTensor::print(std::ostream& os)
{
    tensor_->print(os);
}

void
YetiTensor::distribute(const std::string& str)
{
    StringParser parser(str);
    parser.tokenize(",");
    usi distr_indices[NINDEX];
    usi nindex_distr = parser.nentries();
    for (usi i=0; i < nindex_distr; ++i)
    {
        const std::string& substr = parser.substr(i);
        distr_indices[i] = index(substr);
    }

    tensor_->distribute(distr_indices, nindex_distr);
}

void
YetiTensor::accumulate(
    const YetiTensorPtr& parser,
    double scale
)
{
    tensor_->accumulate(
            parser->get_tensor(),
            parser->scale_ * scale,
            parser->perm_
          );

}

void
YetiTensor::accumulate(
    const YetiContractionPtr& yeticxn,
    double scale
)
{
    usi ncxn_entries = yeticxn->nentries();
    if (ncxn_entries != this->nentries())
    {
        if (ncxn_entries != 0) //make sure not a dot product
        {
            cerr << "Contraction " << yeticxn->str()
                 << " has " << yeticxn->nentries() << " target indices "
             << "but tensor " << str()
              << " has " << this->nentries() << endl;
            abort();
        }
    }
    

    //validate that the contraction has the correct index ranges
    Tensor* ltensor = yeticxn->ltensor->get_tensor();
    Tensor* rtensor = yeticxn->rtensor->get_tensor();
    Tensor* ptensor = tensor_.get();

    //determine if there's a presort needed on the product tensor
    PermutationRuntimeParserPtr product_perm_parser = SymmPermute(yeticxn->str());
    Permutation* product_presort = 0;
    if (this->str() != yeticxn->str())
    {
        product_presort = StringParser::get_permutation(product_perm_parser);
    }

    TensorIndexDescr* descr = ptensor->get_block_descr();
    if (ncxn_entries == 0)
    {
        //make sure dot product
        IndexDescr* mydescr = descr->get(0);
        IndexDescr* dotprod_descr = YetiRuntime::get_descr("dotproduct");
        if (mydescr != dotprod_descr)
        {
            cerr << "No cxn entries, but contraction is not a dot product" << endl;
            abort();
        }
    }

    Permutation *left_postsort = 0, *right_postsort = 0, *product_postsort = 0;

    if (yeticxn->left_presort)
    {
        yeticxn->ltensor->get_tensor()->sort(yeticxn->left_presort);
        left_postsort = yeticxn->left_presort->inverse();
    }

    if (yeticxn->right_presort)
    {
        if (yeticxn->rtensor->get_tensor()->intersects(yeticxn->ltensor->get_tensor()))
        {
            if (!yeticxn->left_presort ||
                !yeticxn->left_presort->equals(yeticxn->right_presort))
            {
                cerr << "self-contracting tensor " << yeticxn->rtensor->get_tensor()->get_name()
                        << " but sorts are not equal" << endl;
                yeticxn->left_presort->print(cerr); cerr << endl;
                yeticxn->right_presort->print(cerr); cerr << endl;
                abort();
            }
        }
        //go ahead and sort
        else
        {
            yeticxn->rtensor->get_tensor()->sort(yeticxn->right_presort);
            right_postsort = yeticxn->right_presort->inverse();
        }
    }

    if (product_presort && !product_presort->is_identity())
    {
        ptensor->sort(product_presort);
        product_postsort = product_presort->inverse();
    }

#if DEBUG_CXN_SORT
    if (yeticxn->left_presort)
    {
        dout << "left presort" << endl;
        yeticxn->left_presort->print(dout); dout << endl;
    }
    if (yeticxn->right_presort)
    {
        dout << "right presort" << endl;
        yeticxn->right_presort->print(dout); dout << endl;
    }
    if (product_presort)
    {
        dout << "product presort" << endl;
        product_presort->print(dout); dout << endl;
    }
#endif

    Permutation* lperm = yeticxn->ltensor->get_permutation();
    if (!lperm->is_identity())
    {
        //sort the left tensor
        if (yeticxn->left_presort) //can't build an automatic sort and have a custom sort
        {
            cerr << "Automatically generated sort and custom sort for tensor "
                  << yeticxn->ltensor->get_tensor()->get_name() << " with indices "
                  << yeticxn->ltensor->str()
                  << endl;
        }
        yeticxn->ltensor->get_tensor()->sort(lperm);
    }

    Permutation* rperm = yeticxn->rtensor->get_permutation();
    if (!rperm->is_identity())
    {
        //sort the left tensor
        if (yeticxn->right_presort) //can't build an automatic sort and have a custom sort
        {
            cerr << "Automatically generated sort and custom sort for tensor "
                  << yeticxn->rtensor->get_tensor()->get_name() << " with indices "
                  << yeticxn->rtensor->str()
                  << endl;
        }
        yeticxn->rtensor->get_tensor()->sort(rperm);
    }
    double contraction_scale
        = yeticxn->contraction_scale
        * yeticxn->ltensor->scale_
        * yeticxn->rtensor->scale_
        * scale;

    Tensor* acc_tensor = 0;
    PermutationSetPtr symmetrization_set = yeticxn->symmetrization_set;
    if (symmetrization_set->order() == 0)
    {
        acc_tensor = ptensor;
    }
    else
    {
        PermutationGroup* grp = new PermutationGroup(ptensor->get_block_descr()->nindex());
        acc_tensor = new Tensor("temporary", ptensor->get_block_descr(), grp);
    }

    ContractionPtr cxn = new Contraction(
        contraction_scale,
        ltensor,
        rtensor,
        acc_tensor,
        yeticxn->lindex,
        yeticxn->rindex,
        yeticxn->ltensor->get_permutation(),
        yeticxn->rtensor->get_permutation()
    );

    if (symmetrization_set->order() != 0)
    {
        acc_tensor->get_tensor_grp()->add(cxn->get_default_grp());
        acc_tensor->get_tensor_grp()->close();
    }

    for (usi i=0; i < ncxn_entries; ++i)
    {

        const std::string& label = yeticxn->substr(i);
        IndexDescr* cxndescr = YetiRuntime::get_descr(label);
        IndexDescr* mydescr = descr->get(i);
        const std::string& mylabel = substr(i);
        if (!mydescr->is_equivalent(cxndescr))
        {
            cerr << "Invalid index " << label << " passed to contraction to form " << yeticxn->str()
                << endl;
            abort();
        }
    }

    cxn->build();
    cxn->run();

    if (symmetrization_set->order() != 0)
    {
        symmetrization_set->close();
        PermutationSet::iterator it = symmetrization_set->begin();
        PermutationSet::iterator stop = symmetrization_set->end();
        for ( ; it != stop; ++it)
        {
            ptensor->accumulate(acc_tensor, yeticxn->post_scale, *it);
        }
        delete acc_tensor;
    }
    else if (yeticxn->post_scale != 1.0)
    {
        cerr << "Non-unit post-contraction scale specified with & operator"
              << ", but no symmetrization set has been given!" << endl;
        abort();
    }

    if (left_postsort)
    {
        yeticxn->ltensor->get_tensor()->sort(left_postsort);
    }
    if (right_postsort)
    {
        yeticxn->rtensor->get_tensor()->sort(right_postsort);
    }
    if (product_postsort)
    {
        ptensor->sort(product_postsort);
    }

    if (!lperm->is_identity())
    {
        Permutation* linv = lperm->inverse();
        yeticxn->ltensor->get_tensor()->sort(linv);
    }
    if (!rperm->is_identity())
    {
        Permutation* rinv = rperm->inverse();
        yeticxn->rtensor->get_tensor()->sort(rinv);
    }

}

void
YetiTensor::element_op(ElementOp* op)
{
    tensor_->element_op(op);
}

bool
YetiTensor::equals(const YetiTensorPtr& yt) const
{
    null_check();
    return tensor_->equals(yt->get_tensor());
}

void
YetiTensor::free()
{
    tensor_ = 0;
    perm_ = 0;
}

void
YetiTensor::clear()
{
    free();
}

void
YetiTensor::zero()
{
    null_check();

    tensor_->zero();
}

bool
YetiTensor::is_initialized() const
{
    return (bool) tensor_;
}

bool
YetiTensor::is_null() const
{
    bool nonnull = tensor_;
    return !nonnull;
}

bool
YetiTensor::is_nonnull() const
{
    return tensor_;
}

double
YetiTensor::norm()
{
    return tensor_->norm();
}

Permutation*
YetiTensor::get_permutation() const
{
    return perm_;
}

Tensor*
YetiTensor::get_tensor() const
{
    return tensor_.get();
}

double
YetiTensor::get_scale() const
{
    return scale_;
}

void
YetiTensor::null_check() const
{
    if (!tensor_)
    {
        cerr << "Tensor has not been initialized" << endl;
        abort();
    }
}

bool
YetiTensor::nonzero() const
{
    return tensor_ && tensor_->nonzero();
}

bool
YetiTensor::unique_nonzero() const
{
    return tensor_ && tensor_->nonzero();
}

YetiTensorPtr
YetiTensor::operator()(const std::string& str) const
{
    null_check();

    return new YetiTensor(
        tensor_,
        perm_,
        scale_,
        str
    );
}

YetiTensorPtr
YetiTensor::operator()(
    const std::string& idxA,
    const std::string& idxB
) const
{
    null_check();

    std::string str = idxA + "," + idxB;

    return new YetiTensor(
        tensor_,
        perm_,
        scale_,
        str
    );
}

YetiTensorPtr
YetiTensor::operator()(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC
) const
{
    null_check();
    std::string str = idxA + "," + idxB + "," + idxC;
    return new YetiTensor(
        tensor_,
        perm_,
        scale_,
        str
    );
}

YetiTensorPtr
YetiTensor::operator()(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC,
    const std::string& idxD
) const
{
    null_check();
    std::string str = idxA + "," + idxB + "," + idxC + "," + idxD;
    return new YetiTensor(
        tensor_,
        perm_,
        scale_,
        str
    );
}

void
YetiTensor::set_permutation(Permutation* perm)
{
    perm_ = perm;
}

void
YetiTensor::set_scale(double scale)
{
    scale_ = scale;
}

void
YetiTensor::sort(const std::string& str)
{
    StringParser parser(str);
    parser.tokenize(", ");
    if (parser.nentries() != nentries())
    {
        cerr << "Invalid string specified for sort.  Sort string has wrong number of indices." << endl;
        abort();
    }

    usi indexmap[NINDEX];
    usi idx=0;

    vector<string>::const_iterator it(parser.begin());
    vector<string>::const_iterator stop(parser.end());
    for ( ; it != stop; ++it, ++idx)
    {
        indexmap[idx] = StringParser::index(*it);
    }
    short plus = 1;
    Permutation* p = Permutation::get_permutation(plus, parser.nentries(), indexmap);
    tensor_->sort(p);

    str_ = str;
    StringParser::tokenize(", ");
}

void
YetiTensor::sort(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC
)
{
    std::string str = idxA + "," + idxB + "," + idxC;
    sort(str);
}

void
YetiTensor::sort(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC,
    const std::string& idxD
)
{
    std::string str = idxA + "," + idxB + "," + idxC + "," + idxD;
    sort(str);
}

void 
YetiTensor::fill(TensorElementComputer* filler)
{
    tensor_->fill(filler);
}

void
YetiTensor::fill(const MatrixPtr& matrix)
{
    tensor_->fill(matrix);
}

bool
YetiTensor::equals(const void* vals)
{
    return tensor_->equals(vals);
}

TensorPtr&
YetiTensor::operator *()
{
    return tensor_;
}

void
YetiTensor::operator=(const YetiTensorPtr& yt)
{
    tensor_ = yt->get_tensor();
    scale_ = yt->get_scale();
    perm_ = yt->get_permutation();
    str_ = yt->str();
    tokenize(",");
}

void
YetiTensor::operator=(const YetiTensor& yt)
{
    tensor_ = yt.get_tensor();
    scale_ = yt.get_scale();
    perm_ = yt.get_permutation();
    str_ = yt.str();
    tokenize(",");
}

YetiTensorPtr::YetiTensorPtr()
    : boost::intrusive_ptr<YetiTensor>(0)
{
}

YetiTensorPtr::YetiTensorPtr(YetiTensor* t)
    : boost::intrusive_ptr<YetiTensor>(t)
{
}

void
YetiTensorPtr::assign(const YetiTensorPtr& tensor)
{
    boost::intrusive_ptr<YetiTensor> me = get();
    me = tensor.get();
}

void
YetiTensorPtr::operator*=(double factor)
{
    get()->get_tensor()->scale(factor);
}

void
YetiTensorPtr::operator-=(const YetiTensorPtr& parser)
{
    get()->accumulate(parser, -1.0);
}

void
YetiTensorPtr::operator+=(const YetiTensorPtr& parser)
{
    get()->accumulate(parser, 1.0);
}

void
YetiTensorPtr::operator+=(const PermutationRuntimeParserPtr& perm)
{
#if ALLOW_PERMUTATIONS
    StringParser* parser = get();
    Permutation* p = parser->get_permutation(perm);
    if (!p)
    {
        cerr << "Invalid permutation parser " << perm->str() << " for tensor " << parser->str() << endl;
        abort();
    }

    PermutationGroup* grp = get()->get_tensor()->get_tensor_grp();

    grp->add(p);
    grp->close();
    get()->get_tensor()->recompute_permutation();
#endif
}

void
YetiTensorPtr::operator|=(const YetiTensorPtr& tensor)
{
    Tensor* me = get()->get_tensor();
    BlockRetrieveAction* action = new AccumulateBlockRetrieveAction(tensor);
    me->configure(action);
}

void
YetiTensorPtr::operator|=(const YetiContractionPtr& cxn)
{
    Tensor* me = get()->get_tensor();
    BlockRetrieveAction* action = new ContractionBlockRetrieveAction(cxn, me);
    me->configure(action);
}

void
YetiTensorPtr::operator|=(ElementOp* op)
{
    Tensor* me = get()->get_tensor();
    op->configure(me);
    BlockRetrieveAction* action = new ElementOpRetrieveAction(op);
    me->configure(action);
}

void
YetiTensorPtr::operator+=(const YetiContractionPtr& yeticxn)
{
    get()->accumulate(yeticxn, 1.0);
}

void
YetiTensorPtr::operator-=(const YetiContractionPtr& yeticxn)
{
    get()->accumulate(yeticxn, -1.0);
}

YetiTensor*
YetiTensorPtr::operator->() const
{
    return get();
}

void
YetiTensorPtr::operator<<=(const YetiTensorPtr& yt)
{
    //validate index ranges
    YetiTensor& dst = *(get());
    YetiTensor& src = *(yt.get());
    usi ntarget_idx = dst.nentries();
    usi nidx = src.nentries();
    usi ncxn_idx = nidx - ntarget_idx;

    if (nidx <= ntarget_idx)
    {
        cerr << "Invalid internal contraction" << endl;
        abort();
    }

    //ensure the indices are the same between tensors
    for (usi i=0; i < ntarget_idx; ++i)
    {
        IndexDescr* dst_descr = dst->get_block_descr()->get(i);
        IndexDescr* src_descr = src->get_block_descr()->get(i);
        if (dst_descr != src_descr)
        {
            cerr << "Invalid internal contraction" << endl;
            abort();
        }
    }

    IndexDescr* cxn_descr = src->get_block_descr()->get(ntarget_idx);
    const std::string& cxn_idx_str = src.substr(ntarget_idx);
    for (usi i=ntarget_idx; i < nidx; ++i)
    {
        IndexDescr* descr = src->get_block_descr()->get(i);
        if (descr != cxn_descr)
        {
            cerr << "Invalid internal contraction" << endl;
            abort();
        }

        if (cxn_idx_str != src.substr(i))
        {
            cerr << "Invalid internal contraction" << endl;
            abort();
        }
    }

    MatrixIndex* mindex = new MatrixIndex(ntarget_idx, ncxn_idx, false);
    MatrixConfiguration* config = new MatrixConfiguration(
        mindex,
        src->get_tensor_grp()
    );
    src->internal_contraction(dst.get_tensor(), config);
}


YetiContraction::YetiContraction(
    const std::string &cxnstring,
    const YetiTensorPtr& _ltensor,
    const YetiTensorPtr& _rtensor
)
    :
    StringParser(cxnstring),
    contraction_scale(1.0),
    ltensor(_ltensor),
    rtensor(_rtensor),
    left_presort(0),
    right_presort(0),
    post_scale(1.0)
{
    tokenize(", ");
}

YetiContraction::~YetiContraction()
{
    ltensor = 0;
    rtensor = 0;

}

YetiContractionPtr::YetiContractionPtr(YetiContraction *cxn)
  :   boost::intrusive_ptr<YetiContraction>(cxn)
{
}

void
YetiContractionPtr::operator=(const YetiTensorPtr& tensor)
{
    Tensor* Atensor = get()->ltensor->get_tensor();
    Tensor* Xtensor = get()->rtensor->get_tensor();
    Tensor* Btensor = tensor->get_tensor();
    //solve a linear system of equations
    DoubleArrayGetOp* Bop = new DoubleArrayGetOp;
    Btensor->element_op(Bop);

    UpperTriangleGetOp* Aop = new UpperTriangleGetOp;
    Atensor->element_op(Aop);

    int n = Atensor->get_block_descr()->get(0)->nelements_data();
    int nelements_B = Btensor->get_totalsize();
    int nrhs = nelements_B / n;
    int ntri = n*(n+1)/2;
    const double* A = Aop->data();
    const double* B = Bop->data();
    double* AP = new double[ntri];
    int* ipiv = new int[n];
    double* X = new double[n*nrhs];
    ::memset(X,0,n*nrhs*sizeof(double));
    double* ferr = new double[nrhs];
    double* berr = new double[nrhs];
    double* work = new double[3*n];
    int* iwork = new int[n];

    char fact = 'N';
    char uplo = 'U';
    double rcond = 0.0;
    int info = 0;


    
    /** Both of these routines don't fucking work for some cases. Fuck LAPACK. How the fuck
        do these routines not fucking work for aug-TZ but do for TZ? */
    //dspsvx(&fact, &uplo, &n, &nrhs, A, AP, ipiv, B, &n, X, &n, &rcond, ferr, berr, work, iwork, &info);
    //dsptrs_(&uplo, &n, &nrhs, A, ipiv, X, &n, &info);

    if (info != 0)
    {
        raise(SanityCheckError, "DSPSVX argument error");
    }

#if 0
    uli idx = 0;
    for (int p=0; p < n; ++p)
    {
        for (int q=0; q <= p; ++q, ++idx)
        {
            dout << stream_printf("J(%d,%d)->J(%d) = %14.10f", p, q, idx, A[idx]) << endl;
        }
    }
#endif

    for (int i=0; i < nrhs; ++i)
    {
        double f = ferr[i];
        double b = berr[i];
        if (fabs(f) > 1e-8 )//|| fabs(b) > 1e-8)
        {
            cerr << stream_printf("FERR[%d] = %12.4e  BERR[%d] = %12.4e", i, f, b) << endl;
            raise(SanityCheckError, "DSPSVX numerical stability error");
        }
    }

    int ni = 1;
    int na = 4;
#if 0
    const double* Bptr = B;
    for (int i=0; i < ni; ++i)
    {
        for (int a=0; a < na; ++a)
        {
            for (int p=0; p < n; ++p, ++Bptr)
            {
                dout << stream_printf("B(%d|%d,%d) = %14.10f",p,i,a,*Bptr) << endl;
            }
        }
    }
#endif
#if 0
    const double* Xptr = X;
    for (int i=0; i < ni; ++i)
    {
        for (int a=0; a < na; ++a)
        {
            for (int p=0; p < n; ++p, ++Xptr)
            {
                dout << stream_printf("X(%d|%d,%d) = %14.10f",p,i,a,*Xptr) << endl;
            }
        }
    }
#endif

    Xtensor->fill(X);

    delete[] AP;
    delete[] ipiv;
    //delete[] X;
    delete[] ferr;
    delete[] berr;
    delete[] work;
    delete[] iwork;
}

template <typename data_t>
data_t
YetiContractionPtr::dot_product() const
{
    //build a tensor with the dot product index descr
    //this builds only a single element!
    YetiTensor dotprod("dot product tensor", "dotproduct");
    dotprod->configure(Tensor::in_core);




    //create an empty tensor
    get()->incref();
    dotprod.accumulate(get(), 1.0);
    get()->decref();

    YetiContraction* cxn = get();
    double scale = cxn->contraction_scale;
    double lscale = cxn->ltensor->get_scale();
    double rscale = cxn->rtensor->get_scale();

    uli idx = 0;
    TensorBlock* block = dotprod->get_block(idx);
    if (!block) //hmm... nothing
        return 0;
        
    DataNode* data = block->get_first_data_node();
    if (data)
    {
        data_t* dataptr = reinterpret_cast<data_t*>(data->data());
        data_t value = *dataptr;
        return value;// * scale * lscale * rscale;
    }
    else
    {
        //no data generated... sooo.... zero
        return 0;
    }
}

YetiContractionPtr::operator int() const
{
    return dot_product<int>();
}

YetiContractionPtr::operator double() const
{
    return dot_product<double>();
}

YetiContractionPtr::operator quad() const
{
    return dot_product<quad>();
}

YetiContractionPtr
yeti::operator*(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor)
{
    std::string cxnstring;
    bool build_left_matrix_as_transpose = false;
    MatrixIndexPtr lindex, rindex;
    Permutation *left_presort, *right_presort;

    get_matrix_index(
        ltensor, rtensor,
        cxnstring,
        lindex, rindex,
        left_presort, right_presort
    );

    usi nprodidx_left = lindex->is_transpose() ?
                           lindex->nrowindex() :
                           lindex->ncolindex();

    usi nprodidx_right = rindex->is_transpose() ?
                           rindex->ncolindex() :
                           rindex->nrowindex();

    YetiContraction* cxn = new YetiContraction(cxnstring, ltensor, rtensor);

    usi nprodidx = cxn->nentries();
    PermutationSet* symmetrization_set = new PermutationSet(nprodidx);

    cxn->lindex = lindex;
    cxn->rindex = rindex;
    cxn->symmetrization_set = symmetrization_set;
    cxn->left_presort = left_presort;
    cxn->right_presort = right_presort;

    return cxn;
}

YetiTensorPtr
yeti::operator->*(
    const PermutationRuntimeParserPtr& perm,
    const YetiTensorPtr& tensor
)
{
    StringParser* strparser = tensor.get();
    Permutation* newperm = strparser->get_permutation(perm);
    tensor->set_permutation(newperm);
    return tensor;
}

YetiContractionPtr
yeti::operator&(
    const PermutationRuntimeParserPtr& perm,
    const YetiContractionPtr& cxn
)
{
    StringParser* strparser = cxn.get();
    Permutation* newperm = strparser->get_permutation(perm);
    if (!newperm)
    {
        cerr << "Invalid permutation parser " << perm->str() << " for string " << cxn->str() << endl;
        abort();
    }

    cxn->symmetrization_set->add(newperm);
    return cxn;
}

std::list<PermutationRuntimeParserPtr>&
yeti::operator&(
    const PermutationRuntimeParserPtr& p,
    const PermutationRuntimeParserPtr& q
)
{
    cxn_symmetrization_list.clear();
    cxn_symmetrization_list.push_back(p);
    cxn_symmetrization_list.push_back(q);
    return cxn_symmetrization_list;
}

YetiContractionPtr
yeti::operator &(
    const std::list<PermutationRuntimeParserPtr>& symm_set,
    const YetiContractionPtr& cxn
)
{
    std::list<PermutationRuntimeParserPtr>::const_iterator it = symm_set.begin();
    std::list<PermutationRuntimeParserPtr>::const_iterator stop = symm_set.end();
    StringParser* strparser = cxn.get();
    for ( ; it != stop; ++it)
    {
        Permutation* newperm = strparser->get_permutation(*it);
        cxn->symmetrization_set->add(newperm);
    }
    return cxn;
}

YetiTensorPtr
yeti::operator*(double scale, const YetiTensorPtr& tensor)
{
    double new_scale = tensor->get_scale() * scale;
    tensor->set_scale(new_scale);
    return tensor;
}

YetiContractionPtr
yeti::operator*(double scale, const YetiContractionPtr& cxn)
{
    cxn->contraction_scale *= scale;
    return cxn;
}

YetiContractionPtr
yeti::operator&&(double scale, const YetiContractionPtr& cxn)
{
    cxn->post_scale *= scale;
    return cxn;
}

PermutationRuntimeParser::PermutationRuntimeParser(const std::string& str)
    : StringParser(str.substr(1)), sign_(0)
{
    string signstr = str.substr(0,1);
    if (signstr == "+")
        sign_ = +1;
    else if (signstr == "-")
        sign_ = -1;
    else
        raise(SanityCheckError, "Invalid string in permutation parser. Need a plus or minus sign.");

    tokenize(",");
}

short
PermutationRuntimeParser::sign() const
{
    return sign_;
}

PermutationRuntimeParserPtr
yeti::Permute(const std::string &str)
{
    return new PermutationRuntimeParser(str);
}

PermutationRuntimeParserPtr
yeti::Permute(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC,
    const std::string& idxD
)
{
    std::string sign = "+";
    std::string str = sign + idxA + "," + idxB + "," + idxC + "," + idxD;
    return new PermutationRuntimeParser(str);
}


PermutationRuntimeParserPtr
yeti::SymmPermute(const std::string &str)
{
    std::string plus_str = "+";
    plus_str += str;
    return new PermutationRuntimeParser(plus_str);
}

void
yeti::get_matrix_index(
    const YetiTensorPtr& ltensor,
    const YetiTensorPtr& rtensor,
    std::string& cxnstring,
    MatrixIndexPtr& lindex,
    MatrixIndexPtr& rindex,
    Permutation*& left_presort,
    Permutation*& right_presort
)
{
    usi nidx_target_left = 0, nidx_target_right = 0;
    usi nidx_cxn = 0;

    usi cxn_indices_left[NINDEX];
    usi cxn_indices_right[NINDEX];
    usi target_indices_left[NINDEX];
    usi target_indices_right[NINDEX];

    vector<string>::const_iterator it(ltensor->begin());
    vector<string>::const_iterator stop(ltensor->end());
    usi lidx = 0, ridx = 0;
    for ( ; it != stop; ++it, ++lidx)
    {
        ridx = rtensor->index(*it);
        if (ridx >= rtensor->nentries())
        {
            target_indices_left[nidx_target_left] = lidx;
            ++nidx_target_left;
            if (cxnstring.size() == 0)
                cxnstring += *it;
            else
                cxnstring += "," + *it;
        }
        else
        {

            cxn_indices_left[nidx_cxn] = lidx;
            cxn_indices_right[nidx_cxn] = ridx;
            ++nidx_cxn;
        }
    }

    it = rtensor->begin();
    stop = rtensor->end();
    ridx = 0;
    for ( ; it != stop; ++it, ++ridx)
    {
        uli lidx = ltensor->index(*it);
        if (lidx >= ltensor->nentries())
        {
            target_indices_right[nidx_target_right] = ridx;
            ++nidx_target_right;

            if (cxnstring.size() == 0)
                cxnstring += *it;
            else
                cxnstring += "," + *it;
        }
    }

#if 0
    if (nidx_cxn == 0)
    {
        cerr << "No contraction indices in contraction of "
            << "(" << ltensor->str() << ")"
            << " ->* "
            << "(" << rtensor->str() << ")" << endl;

        abort();
    }
#endif

    bool sort_left = false;
    for (usi i=1; i < nidx_cxn; ++i)
    {
        if (cxn_indices_left[i] != cxn_indices_left[i-1] + 1) //we must sort left tensor
        {
            sort_left = true;
            break;
        }
    }
    for (usi i=1; i < nidx_target_left; ++i)
    {
        if (target_indices_left[i] != target_indices_left[i-1] + 1) //we must sort left tensor
        {
            sort_left = true;
            break;
        }
    }

    bool sort_right = false;
    for (usi i=1; i < nidx_cxn; ++i)
    {
        if (cxn_indices_right[i] != cxn_indices_right[i-1] + 1) //we must sort left tensor
        {
            sort_right = true;
            break;
        }
    }
    for (usi i=1; i < nidx_target_right; ++i)
    {
        if (target_indices_right[i] != target_indices_right[i-1] + 1) //we must sort left tensor
        {
            sort_right = true;
            break;
        }
    }

    usi pmap[NINDEX];
    short plus = 1;
    if (sort_left)
    {
        usi idx = 0;
        for (usi i=0; i < nidx_target_left; ++i, ++idx)
        {
            pmap[idx] = target_indices_left[i];
        }
        for (usi i=0; i < nidx_cxn; ++i, ++idx)
        {
            pmap[idx] = cxn_indices_left[i];
        }
        usi nidx_left = nidx_target_left + nidx_cxn;
        left_presort = Permutation::get_permutation(plus, nidx_left, pmap);
        lindex = new MatrixIndex(nidx_target_left, nidx_cxn, false);
    }
    else if (nidx_target_left == 0)
    {
        lindex = new MatrixIndex(0, nidx_cxn, false);
        left_presort = 0;
    }
    else if (nidx_cxn == 0)
    {
        lindex = new MatrixIndex(nidx_target_left, 0);
        left_presort = 0;
    }
    else
    {
        usi first_cxn_idx = cxn_indices_left[0];
        if (first_cxn_idx == 0)
            lindex = new MatrixIndex(nidx_cxn, nidx_target_left, true);
        else
            lindex = new MatrixIndex(nidx_target_left, nidx_cxn, false);
        left_presort = 0;
    }


    if (sort_right)
    {
        usi idx = 0;
        for (usi i=0; i < nidx_target_right; ++i, ++idx)
        {
            pmap[idx] = target_indices_right[i];
        }
        for (usi i=0; i < nidx_cxn; ++i, ++idx)
        {
            pmap[idx] = cxn_indices_right[i];
        }
        usi nidx_right = nidx_target_right + nidx_cxn;
        right_presort = Permutation::get_permutation(plus, nidx_right, pmap);
        rindex = new MatrixIndex(nidx_target_right, nidx_cxn, true);
    }
    else if (nidx_target_right == 0)
    {
        rindex = new MatrixIndex(0, nidx_cxn, true);
        right_presort = 0;
    }
    else if (nidx_cxn == 0)
    {
        rindex = new MatrixIndex(nidx_target_right, 0, true);
        right_presort = 0;
    }
    else
    {
        usi first_cxn_idx = cxn_indices_right[0];
        if (first_cxn_idx == 0)
            rindex = new MatrixIndex(nidx_cxn, nidx_target_right, false);
        else
            rindex = new MatrixIndex(nidx_target_right, nidx_cxn, true);
        right_presort = 0;
    }
}

YetiMatrix::YetiMatrix(
    const YetiTensor &yt,
    const std::string &rows,
    const std::string &cols
) :
    tensor_(0),
    config_(0),
    rowparser_(rows),
    colparser_(cols)
{
    std::string total = rows + "," + cols;
    tensor_ = yt(total);

    rowparser_.tokenize(",");
    colparser_.tokenize(",");

    //not transposed
    uli nrows = rowparser_.nentries();
    uli ncols = colparser_.nentries();

    MatrixIndex* index = new MatrixIndex(nrows, ncols, false);
    config_ = new MatrixConfiguration(
                    index,
                    tensor_->get_tensor()->get_tensor_grp()
              );

    TensorIndexDescr* tensor_descr = tensor_->get_tensor()->get_block_descr();
    for (usi i=0; i < nrows; ++i)
    {
        IndexDescr* mydescr = YetiRuntime::get_descr(rowparser_.substr(i));
        IndexDescr* test_descr = tensor_descr->get(i);
        if (mydescr != test_descr)
        {
            cerr << "Invalid matrix configuration in building yeti matrix" << endl;
            abort();
        }
    }

    for (usi i=0; i < ncols; ++i)
    {
        IndexDescr* mydescr = YetiRuntime::get_descr(colparser_.substr(i));
        IndexDescr* test_descr = tensor_descr->get(i+nrows);
        if (mydescr != test_descr)
        {
            cerr << "Invalid matrix configuration in building yeti matrix" << endl;
            abort();
        }
    }
}

void
YetiMatrix::get_matrix(RectMatrixPtr& matrix)
{
    tensor_->get_tensor()->get_matrix(matrix, config_.get());
}

void
YetiMatrix::get_matrix(SymmMatrixPtr& matrix)
{
    tensor_->get_tensor()->get_matrix(matrix, config_.get());
}

void
YetiMatrix::accumulate(RectMatrixPtr& matrix)
{
    tensor_->get_tensor()->accumulate(matrix, config_.get());
}

void
YetiMatrix::accumulate(SymmMatrixPtr& matrix)
{
    tensor_->get_tensor()->accumulate(matrix, config_.get());
}


ostream&
yeti::operator<<(std::ostream& os, const YetiTensorPtr& yt)
{
    yt->get_tensor()->print(os);
    return os;
}

ostream&
yeti::operator<<(std::ostream& os, const YetiContractionPtr& ytcxn)
{
    YetiTensor t("temporary intermediate", ytcxn->str());
    t(ytcxn->str()) += ytcxn;
    t.get_tensor()->print(os);
    return os;
}

double&
yeti::operator+=(double& e, const YetiContractionPtr& cxn)
{
    double tmp = cxn;
    e += tmp;
    return e;
}

double&
yeti::operator-=(double& e, const YetiContractionPtr& cxn)
{
    double tmp = cxn;
    e -= tmp;
    return e;
}

YetiContractionPtr
yeti::operator*(const YetiContractionPtr& cxn, const YetiTensorPtr& rtensor)
{
    YetiTensorPtr ltensor = new YetiTensor("temporary cxn tensor", cxn->str());
    ltensor += cxn;
    YetiContractionPtr newcxn = ltensor * rtensor;
    return newcxn;
}

YetiContractionPtr
yeti::operator*(const YetiTensorPtr& ltensor, const YetiContractionPtr& cxn)
{
    YetiTensorPtr rtensor = new YetiTensor("temporary cxn tensor", cxn->str());
    rtensor += cxn;
    YetiContractionPtr newcxn = ltensor * rtensor;
    return newcxn;
}
