#include "tensorparser.h"
#include "index.h"
#include "runtime.h"
#include "permutation.h"
#include "tensor.h"
#include "malloc.h"
#include "matrix.h"
#include "exception.h"
#include "contraction.h"
#include "dataimpl.h"

#include <vector>

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

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

    cerr << "index string " << str << " is not contained in " << str_ << endl;
    abort();
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
    perm_(0),
    tensor_name_(name),
    iter_(0)
{
    tokenize(", ");
    IndexRangeTuple* tuple = YetiRuntime::get_index_tuple(str_);

    usi nindex = tuple->nindex();
    PermutationGroup* grp = new PermutationGroup(nindex);
    perm_ = grp->get_identity();

    incref();
    if (p1) grp->add(p1->get_permutation(this));
    if (p2) grp->add(p2->get_permutation(this));
    if (p3) grp->add(p3->get_permutation(this));
    if (p4) grp->add(p4->get_permutation(this));
    if (p5) grp->add(p5->get_permutation(this));
    decref();

    grp->close();

    tensor_ = new Tensor(tensor_name_, tuple, grp);
}

YetiTensor::YetiTensor(
    const TensorPtr& tensor,
    const PermutationPtr& perm,
    double scale,
    const std::string& str
) : tensor_(tensor),
    perm_(perm),
    scale_(scale),
    StringParser(str),
    tensor_name_(tensor->config()->name),
    iter_(0)
{
}

YetiTensor::YetiTensor(
    const YetiTensorPtr& tensor,
    const PermutationPtr& p
) : tensor_(tensor->tensor_),
    perm_(p),
    scale_(1.0),
    StringParser(tensor->str_),
    tensor_name_(tensor->tensor_name_),
    iter_(0)
{
    if (!tensor->perm_->is_identity())
    {
        cerr << "Permutation specifed for tensor which is"
                "already permuted.  Permutations should only"
                "be specified for base tensor" << endl;
        abort();
    }

    tokenize(", ");
}

YetiTensor::YetiTensor(const TensorPtr &tensor)
    : StringParser(""),
    tensor_(tensor),
    perm_(tensor->get_permutation_grp()->get_identity()),
    scale_(1.0),
    tensor_name_(tensor->config()->name)
{
}

YetiTensor::YetiTensor(const YetiTensorPtr& tensor)
    : StringParser(tensor->str()),
    tensor_(0),
    perm_(0),
    scale_(1.0),
    tensor_name_(tensor->tensor_name_)
{
    tokenize(", ");
    IndexRangeTuple* tuple = YetiRuntime::get_index_tuple(str_);

    usi nindex = tuple->nindex();
    PermutationGroupPtr grp = tensor->tensor_->get_permutation_grp();
    perm_ = grp->get_identity();
    TensorConfigurationPtr config
        = tensor->tensor_->config()->copy(tensor_name_);
    tensor_ = new Tensor(tensor_name_, tuple, grp, config);
    tensor_->accumulate(
                tensor->tensor_,
                tensor->perm_,
                tensor->scale_
              );
}

YetiTensor::YetiTensor()
    :
    tensor_(0),
    StringParser(""),
    perm_(0),
    scale_(1.0),
    tensor_name_(""),
    iter_(0)
{
}

void
YetiTensor::operator=(const YetiTensorPtr &tensor)
{
    if (tensor_name_ == tensor->tensor_name_)
        raise(SanityCheckError, "cannot copy tensors to same name");

    tensor_->accumulate(tensor->tensor_, tensor->perm_, tensor->scale_);
}

void
YetiTensor::operator=(const YetiTensor& tensor)
{
    raise(SanityCheckError, "assignment operator should never be called on YetiTensor type");
}


YetiTensor::YetiTensor(const YetiTensor &tensor)
    : StringParser("")
{
    raise(SanityCheckError, "Copy constructor for YetiTensor should never be called");
}

YetiTensor::YetiTensor(
    const std::string &str,
    YetiTensor* t
) : StringParser(str),
    tensor_(t->tensor_),
    tensor_name_(t->tensor_name_),
    perm_(t->perm_),
    scale_(1.0),
    iter_(0)
{
    init_runtime_tensor();
}

YetiTensor::YetiTensor(
    const std::string& name,
    const std::string &str,
    YetiTensor* t
) : StringParser(str),
    tensor_(t->tensor_),
    tensor_name_(name),
    perm_(t->perm_),
    scale_(1.0),
    iter_(0)
{
    init_runtime_tensor();
}

YetiTensor::YetiTensor(
    YetiTensor* t,
    const std::string& idx1,
    const std::string& idx2,
    const std::string& idx3,
    const std::string& idx4
)
    : StringParser(""),
    tensor_(t->tensor_),
    tensor_name_(t->tensor_name_),
    perm_(t->perm_),
    scale_(1.0),
    iter_(0)
{
    str_ = idx1 + "," + idx2 + "," + idx3 + "," + idx4;
    entries_.push_back(idx1);
    entries_.push_back(idx2);
    entries_.push_back(idx3);
    entries_.push_back(idx4);
}

YetiTensor::YetiTensor(
    YetiTensor* t,
    const std::string& idx1,
    const std::string& idx2
)
    : StringParser(""),
    tensor_(t->tensor_),
    tensor_name_(t->tensor_name_),
    perm_(t->perm_),
    scale_(1.0),
    iter_(0)
{
    str_ = idx1 + "," + idx2;
    entries_.push_back(idx1);
    entries_.push_back(idx2);
}

void
YetiTensor::init_runtime_tensor()
{
    tokenize(", ");
     if (entries_.size() != tensor_->nindex())
        raise(SanityCheckError, "invalid string passed to yeti tensor");

    //verify the index ranges
    IndexRangeTuple* tuple = tensor_->get_index_ranges();

    vector<string>::const_iterator itstring(entries_.begin());
    IndexRangeTuple::iterator it(tuple->begin());
    IndexRangeTuple::iterator stop(tuple->end());
    bool subtensor = false; //assume we are asking for the exact tensor
    IndexRangeTuplePtr newtuple = new IndexRangeTuple(tuple->nindex());
    usi idx = 0;
    for ( ; it != stop; ++it, ++itstring, ++idx)
    {
        IndexRange* subrange = YetiRuntime::get_range(*itstring);
        IndexRange* parent_range = *it;
        newtuple->set(idx, subrange);
        if (subrange != parent_range)
        {
            subtensor = true;
        }
    }

    if (subtensor)
    {
        if (tensor_->name() != tensor_name_)
        {
            cerr << "subtensor name " << tensor_name_
                << " must be the same as parent tensor name " << tensor_->name()
                << endl;
           abort();
        }

        TensorPtr parent = tensor_;
        tensor_ = new Tensor(
                        tensor_name_,
                        newtuple,
                        parent->get_permutation_grp(),
                        tensor_->config()
                       );
        tensor_->config()->tile_map_builder = new SubtensorTileMapBuilder(parent.get(), newtuple);
    }
}


ContractionPtr
YetiTensor::accumulate_tasks(const YetiContractionPtr& yeticxn)
{
    if (yeticxn->nentries() != this->nentries())
    {
        cerr << "Contraction has " << yeticxn->nentries() << " target indices "
             << "but tensor has " << this->nentries() << endl;
        abort();
    }

    //validate that the contraction has the correct index ranges
    IndexRangeTuple* tuple = tensor_->get_index_ranges();
    for (usi i=0; i < yeticxn->nentries(); ++i)
    {
        const std::string& label = yeticxn->substr(i);
        IndexRange* cxnrange = YetiRuntime::get_range(label);
        IndexRange* myrange = tuple->get(i);
        if (myrange != cxnrange)
        {
            cerr << "Invalid index " << label << " passed to contraction to form " << yeticxn->str()
                << endl;
            abort();
        }
    }

    //create a default tensor configuration for now
    double scale_factory
        = yeticxn->scale * yeticxn->ltensor->scale_ * yeticxn->rtensor->scale_;
    ContractionPtr cxn = new Contraction(
        scale_factory,
        yeticxn->ltensor->get_tensor(),
        yeticxn->rtensor->get_tensor(),
        yeticxn->lindex,
        yeticxn->rindex,
        yeticxn->required_grp,
        yeticxn->ltensor->get_permutation(),
        yeticxn->rtensor->get_permutation(),
        tensor_
    );

    MatrixPtr lmatrix = cxn->get_left_matrix();
    MatrixPtr rmatrix = cxn->get_right_matrix();
    MatrixPtr pmatrix = cxn->get_product_matrix();

    str_ = yeticxn->str();

    perm_ = tensor_->get_permutation_grp()->get_identity();

    return cxn;
}

void
YetiTensor::allocate_iterator()
{
    iter_ = tensor_->get_iterator();
}

Tile**
YetiTensor::begin() const
{
    if (!iter_)
    {
        cerr << "iterator not yet allocated for yeti tensor" << endl;
        abort();
    }
    return iter_->begin();
}

Tile**
YetiTensor::end() const
{
    if (!iter_)
    {
        cerr << "iterator not yet allocated for yeti tensor" << endl;
        abort();
    }
    return iter_->end();
}

usi
YetiTensor::get_nindex_front(const std::string &str)
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

Tensor*
YetiTensor::get_tensor() const
{
    return tensor_.get();
}

PermutationPtr
YetiTensor::get_permutation() const
{
    return perm_;
}

TensorPtr
YetiTensor::parameterize(const std::string& str)
{
    usi nindex = get_nindex_front(str);
    TensorPtr param_tensor = tensor_->get_parameterized_tensor(nindex);
    std::string newstr = register_new_tuple(param_tensor->get_index_ranges());
    return param_tensor;
}

std::string
YetiTensor::register_new_tuple(const IndexRangeTuplePtr &tuple)
{
    std::string new_tensor_str;
    for (usi i=0; i < tuple->nindex(); ++i)
    {
        IndexRange* range = tuple->get(i);
        std::string new_idx_str = stream_printf("index%p", range);
        YetiRuntime::register_index_range(new_idx_str, range);

        if (i==0)
        {
            new_tensor_str += new_idx_str;
        }
        else
        {
            new_tensor_str += ",";
            new_tensor_str += new_idx_str;
        }
    }

    return new_tensor_str;
}

void
YetiTensor::scale(double scale)
{
    scale_ *= scale;
}

void
YetiTensor::sort(const std::string &str)
{
    if (str.size() != str_.size())
    {
        cerr << "Invalid string specified for sort.  Sort string has wrong number of indices." << endl;
        abort();
    }

    vector<string> index_list;
    smartptr::split(index_list, str, ",");
    sort(index_list);
}

void
YetiTensor::sort(
    const std::string& idx1,
    const std::string& idx2,
    const std::string& idx3,
    const std::string& idx4
)
{
    vector<string> index_list;
    index_list.push_back(idx1);
    index_list.push_back(idx2);
    index_list.push_back(idx3);
    index_list.push_back(idx4);
    sort(index_list);
}

void
YetiTensor::sort(const std::vector<std::string>& index_list)
{
    usi* indexmap = yeti_malloc_perm();
    usi idx=0;
    vector<string>::const_iterator it(index_list.begin());
    vector<string>::const_iterator stop(index_list.end());
    for ( ; it != stop; ++it, ++idx)
    {
        indexmap[idx] = StringParser::index(*it);
    }

    PermutationPtr p = new Permutation(indexmap, tensor_->nindex(), 1);
    tensor_->sort(p);
}

TensorPtr
YetiTensor::split(const std::string &str)
{
    usi nindex = get_nindex_front(str);

    TensorPtr split_tensor = tensor_->split(nindex);
    IndexRangeTuplePtr newtuple = split_tensor->get_index_ranges();
    std::string new_tensor_str = register_new_tuple(newtuple);
    return split_tensor;
}

Tensor*
YetiTensor::operator ->() const
{
    return tensor_.get();
}

void
YetiTensor::operator*=(double factor)
{
    tensor_->scale(factor);
}

void
YetiTensor::operator=(const YetiContractionPtr& yeticxn)
{
    //do not allocate yet... let the contraction allocate
    ContractionPtr cxn = accumulate_tasks(yeticxn);

    cxn->run();
    tensor_->update();
}

void
YetiTensor::operator+=(const YetiContractionPtr& yeticxn)
{
    ContractionPtr cxn = accumulate_tasks(yeticxn);

    cxn->run();
    tensor_->update();
}

void
YetiTensor::operator-=(const YetiContractionPtr& yeticxn)
{
    yeticxn->scale *= -1;
    ContractionPtr cxn = accumulate_tasks(yeticxn);

    cxn->run();
    tensor_->update();
}

void
YetiTensor::free()
{
    scale_ = 0;
    tensor_ = 0;
    perm_ = 0;
    tensor_name_ = "";
}

void
YetiTensor::operator -=(const YetiTensorPtr& yt)
{
    tensor_->accumulate(yt->tensor_, yt->perm_, -1.0 * yt->scale_);
}

void
YetiTensor::operator +=(const YetiTensorPtr& yt)
{
    tensor_->accumulate(yt->tensor_, yt->perm_, yt->scale_);
}

YetiTensorPtr
YetiTensor::operator()(const std::string& str) const
{
    return new YetiTensor(
        str,
        const_cast<YetiTensor*>(this) //this is safe - it keeps from passing a bazillion arguments
    );
}

MatrixPtr
YetiTensor::matrix(const std::string& rowstr, const std::string& colstr) const
{
    StringParser rowparser(rowstr); rowparser.tokenize(", ");
    StringParser colparser(colstr); colparser.tokenize(", ");

    if (rowparser.nentries() == 0)
    {
        raise(SanityCheckError, "cannot build matrix with no rows");
    }
    if (colparser.nentries() == 0)
    {
        raise(SanityCheckError, "cannot build matrix with no cols");
    }

    usi rowstart = this->index(rowparser.substr(0));
    usi nrows = rowparser.nentries();
    usi prev = rowstart;
    for (usi i=1; i < nrows; ++i)
    {
        usi next = this->index(rowparser.substr(i));
        if (next != prev + 1)
        {
            raise(SanityCheckError, "non-consecutive matrix indices in row declaration");
        }
        prev = next;
    }


    usi colstart = this->index(colparser.substr(0));
    usi ncols = colparser.nentries();
    prev = colstart;
    for (usi i=1; i < ncols; ++i)
    {
        usi next = this->index(colparser.substr(i));
        if (next != prev + 1)
        {
            raise(SanityCheckError, "non-consecutive matrix indices in column declaration");
        }
        prev = next;
    }

    if (nrows + ncols != this->nentries())
        raise(SanityCheckError, "number of matrix indices is incorrect");


    MatrixIndex::matrix_index_t rowtype = rowstart == 0 ? MatrixIndex::front : MatrixIndex::back;
    MatrixIndex::matrix_index_t coltype = colstart == 0 ? MatrixIndex::front : MatrixIndex::back;

    MatrixIndex* index = new MatrixIndex(rowtype, nrows, coltype, ncols);
    MatrixConfiguration* config = new MatrixConfiguration(
                                            index,
                                            tensor_->get_permutation_grp(),
                                            tensor_->get_permutation_grp()->get_identity()
                                       );

    PermutationGroupPtr cxngrp = config->get_full_row_permutation_grp();
    config->configure_quotient_set(cxngrp);
    config->configure_isotropy(cxngrp);
    Matrix* matrix = new Matrix(tensor_, config);
    matrix->set_as_multiplicand();
    return matrix;
}

YetiTensorPtr::YetiTensorPtr(YetiTensor *t)
    : boost::intrusive_ptr<YetiTensor>(t)
{
}

YetiTensorPtr::YetiTensorPtr()
    : boost::intrusive_ptr<YetiTensor>(0)
{
}

void
YetiTensorPtr::operator=(const YetiContractionPtr& yeticxn)
{
    *(get()) = yeticxn;
}

void
YetiTensorPtr::operator+=(const YetiContractionPtr& yeticxn)
{
    *(get()) += yeticxn;
}

void
YetiTensorPtr::operator-=(const YetiContractionPtr& yeticxn)
{
    *(get()) -= yeticxn;
}

YetiTensorPtr&
YetiTensorPtr::operator[](const std::string& str)
{
    get()->tensor_name_ = str;
    return *this;
}

YetiSubsetTensor::YetiSubsetTensor(
    const YetiTensorPtr &tensor
) : YetiTensor(tensor->str(), tensor.get())
{
}

PermutationRuntimeParser::PermutationRuntimeParser(const std::string &str)
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

PermutationPtr
PermutationRuntimeParser::get_permutation(const StringParserPtr& parser)
{
    usi* indexmap = yeti_malloc_perm();
    usi nindex = parser->nentries();
    for (usi i=0; i < nindex; ++i)
        indexmap[i] = i;

    if (nentries() == 1) //transposition
    {
        vector<std::string> indices;
        smartptr::split(indices, *begin(), "->");
        usi i = parser->index(indices[0]);
        usi j = parser->index(indices[1]);
        indexmap[i] = j;
        indexmap[j] = i;
    }
    else
    {
        vector<std::string> indices;
        vector<std::string>::const_iterator it(begin());
        vector<std::string>::const_iterator stop(end());
        for ( ; it != stop; ++it)
        {
            indices.clear();
            smartptr::split(indices, *it, "->");
            usi i = parser->index(indices[0]);
            usi j = parser->index(indices[1]);
            indexmap[j] = i;
        }
    }

    Permutation* p = new Permutation(indexmap, nindex, sign_);
    return p;
}

YetiContraction::YetiContraction(const std::string &cxnstring)
    :
    StringParser(cxnstring),
    scale(1.0)
{
    tokenize(", ");
}

PermutationRuntimeParserPtr
yeti::P(const std::string &str)
{
    return new PermutationRuntimeParser(str);
}

PermutationRuntimeParserPtr
yeti::S(const std::string &str)
{
    std::string plus_str = "+";
    plus_str += str;
    return new PermutationRuntimeParser(plus_str);
}

MatrixIndexPtr
yeti::get_matrix_index(
    const YetiTensorPtr& ltensor,
    const YetiTensorPtr& rtensor,
    std::string& cxnstring
)
{
    usi nidx_target = 0;
    usi nidx_cxn = 0;
    usi idx = 0;

    const vector<string>& entries = ltensor->entries();

    usi nindex = ltensor->nentries();
    vector<string> permuted_indices(ltensor->nentries());

    PermutationPtr p = ltensor->get_permutation();
    const usi* indexmap = p->indexmap();
    for (usi i=0; i < nindex; ++i)
        permuted_indices[i] = entries[indexmap[i]];

    vector<string>::iterator it(permuted_indices.begin());
    vector<string>::iterator stop(permuted_indices.end());

    usi* cxn_indices = yeti_malloc_perm();

    for ( ; it != stop; ++it, ++idx)
    {
        if (rtensor->contains(*it))
        {
            cxn_indices[nidx_cxn] = idx;
            ++nidx_cxn;
        }
        else
        {
            if (cxnstring.size() == 0)
            {
                cxnstring += *it;
            }
            else
            {
                cxnstring += ",";
                cxnstring += *it;
            }
            ++nidx_target;
        }
    }

    for (usi i=1; i < nidx_cxn; ++i)
    {
        if (cxn_indices[i] != cxn_indices[i-1] + 1)
        {
            cerr << "Contraction indices not consecutive in contraction of "
            << "(" << ltensor->str() << ")" << " ->* " << "(" << rtensor->str() << ")" << endl;


            if (!ltensor->get_permutation()->is_identity())
                cerr << "Left Permutation: " << ltensor->get_permutation() << endl;
            if (!rtensor->get_permutation()->is_identity())
                cerr << "Right Permutation: " << rtensor->get_permutation() << endl;

            abort();

        }



    }

    if (nidx_cxn == 0)
    {
        cerr << "No contraction indices in contraction of "
            << "(" << ltensor->str() << ")" << " ->* " << "(" << rtensor->str() << ")" << endl;

        abort();
    }

    MatrixIndex::matrix_index_t rowtype, coltype;
    if (cxn_indices[0] == 0)
    {
        rowtype = MatrixIndex::front;
        coltype = MatrixIndex::back;
    }
    else
    {
        rowtype = MatrixIndex::back;
        coltype = MatrixIndex::front;
    }

    usi nrows = nidx_cxn;
    usi ncols = nidx_target;

    //if the tensor has parameter indices, these must be included in the col indices
    MatrixIndex* mindex = new MatrixIndex(rowtype, nrows, coltype, ncols);

    yeti_free_perm(cxn_indices);

    return mindex;
}

YetiContractionPtr
yeti::operator*(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor)
{
    std::string cxnstring;
    MatrixIndexPtr lindex = get_matrix_index(ltensor, rtensor, cxnstring);
    MatrixIndexPtr rindex = get_matrix_index(rtensor, ltensor, cxnstring);
    usi nprodidx = lindex->ncolindex() + rindex->ncolindex();
    PermutationGroup* required_grp = new PermutationGroup(nprodidx);

    YetiContraction* cxn = new YetiContraction(cxnstring);
    cxn->ltensor = ltensor;
    cxn->rtensor = rtensor;
    cxn->lindex = lindex;
    cxn->rindex = rindex;
    cxn->required_grp = required_grp;
    cxn->nidx_cxn = lindex->nrowindex();
    cxn->nidx_target = lindex->ncolindex() + rindex->ncolindex();

    return cxn;
}

YetiContractionPtr
yeti::operator&(const PermutationRuntimeParserPtr& parser, const YetiContractionPtr& cxn)
{
    PermutationPtr p = parser->get_permutation(cxn);
    cxn->required_grp->add(p);
    cxn->required_grp->close();
    return cxn;
}

YetiTensorPtr
yeti::operator->*(
    const PermutationRuntimeParserPtr& perm,
    const YetiTensorPtr& tensor
)
{
    PermutationPtr p = perm->get_permutation(tensor);
    return new YetiTensor(tensor, p);
}

YetiTensorPtr
yeti::operator*(double scale, const YetiTensorPtr& tensor)
{
    tensor->scale(scale);
    return tensor;
}

YetiContractionPtr
yeti::operator*(double scale, const YetiContractionPtr& cxn)
{
    cxn->scale *= scale;
    return cxn;
}

std::ostream&
yeti::operator<<(std::ostream& os, const YetiTensor& tensor)
{
    tensor->print(os);
    return os;
}

YetiContractionPtr::YetiContractionPtr(YetiContraction *cxn)
  :   boost::intrusive_ptr<YetiContraction>(cxn)
{
}

void
YetiContractionPtr::operator=(const YetiTensorPtr& tensor)
{
    //do nothing
}

template <typename data_t>
data_t
YetiContractionPtr::dot_product()
{
    YetiContractionPtr yeticxn = get();
    YetiTensor dotproduct("dot product", ""); //no index ranges
    dotproduct->configure(new MemoryBlockFactory<data_t>);

    dotproduct = yeticxn;

    Tile* data_tile = 0;
    if (dotproduct->depth() == 0)
    {
        data_tile = dotproduct.get_tensor();
    }
    else
    {
        TileIteratorPtr iter = dotproduct->get_iterator();
        if (iter->ntiles() != 1)
        {
            cerr << "Dot product produces multiple data tiles" << endl;
            abort();
        }
        data_tile = *(iter->begin());
    }

    DataBlock* dblock = data_tile->get_data();
    dblock->retrieve(NOT_THREADED);
    data_t* dataptr = dblock->data()->template get<data_t>();
    dblock->release(NOT_THREADED);
    data_t e = *dataptr;
    return e;
}

YetiContractionPtr::operator double()
{
    return dot_product<double>();
}

YetiContractionPtr::operator quad()
{
    return dot_product<quad>();
}
