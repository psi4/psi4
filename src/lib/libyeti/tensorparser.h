#ifndef _yeti_tensor_parser_h
#define _yeti_tensor_parser_h

#include <iostream>

#include "class.h"
#include "tensor.hpp"
#include "tile.hpp"
#include "matrix.hpp"
#include "tensorparser.hpp"
#include "contraction.hpp"
#include "permutation.hpp"
#include "index.hpp"

namespace yeti {

class YetiContractionPtr : public boost::intrusive_ptr<YetiContraction> {

    public:
        YetiContractionPtr(YetiContraction* cxn);

        void operator=(const YetiTensorPtr& tensor);

};

/**
    @class StringParser
    Class for parsing/tokenizing strings
*/
class StringParser : public smartptr::Countable {

    protected:
        std::string str_;

        std::vector<std::string> entries_;

    public:
        typedef std::vector<std::string>::const_iterator iterator;

        StringParser(const std::string& str);

        /**
            @return The number of entries after splitting the string on a given token
        */
        uli nentries() const;

        /**
            @param The substring to locate after splitting
            @return The index of the substring in the list of substrings
        */
        uli index(const std::string& str) const;

        /**
            @param str The substring to locate after tokenizing
            @return Whether the substring is contained in the list of substrings
        */
        bool contains(const std::string& str) const;

        /**
            @return The underlying string
        */
        const std::string& str() const;

        void str(const std::string& str);

        const std::string& substr(uli index) const;

        /**
            Split the underlying string into a list of substrings
            based on a given token
            @param token
        */
        void tokenize(const char* token);

        iterator begin() const;

        iterator end() const;

        /**
            Count the number of ocurrences of individual substrings
            @param counts Reference return.  Keys are substrings, values
                          are the number of times string appears
        */
        void count(
            std::map<std::string,uli>& counts
        ) const;

        /**
            @return The vector of substrings after tokenizing
        */
        const std::vector<std::string>& entries() const;
};


/**
    @class YetiTensor
    String parser designed for generating tensor operations
*/
class YetiTensor : public StringParser {

    private:
        TensorPtr tensor_;

        std::string tensor_name_;

        /**
            Permutation of the tensor values that will eventually be used
            in a contraction
        */
        PermutationPtr perm_;

        /**
            Scale factor the tensor that will eventually be used in a contraction
        */
        double scale_;

        ContractionPtr accumulate_tasks(const YetiContractionPtr& cxn);

    protected:
        YetiTensor(
            const std::string &str,
            YetiTensor* t
        );

        YetiTensor(
            const std::string& name,
            const std::string& str,
            YetiTensor* t
        );

        YetiTensor(
            YetiTensor* t,
            const std::string& idx1,
            const std::string& idx2
        );

        YetiTensor(
            YetiTensor* t,
            const std::string& idx1,
            const std::string& idx2,
            const std::string& idx3,
            const std::string& idx4
        );

        void init_runtime_tensor();

        void sort(const std::vector<std::string>& index_list);

    public:
        /**
            @name
            @str The list of indices, e.g. T("i,j,a,b") comma separated
                 that define the tensor
        */
        YetiTensor(
            const std::string& name,
            const std::string& str
        );

        YetiTensor(
            const std::string& name,
            const std::string& idx1,
            const std::string& idx2
        );

        YetiTensor(
            const std::string& name,
            const std::string& idx1,
            const std::string& idx2,
            const std::string& idx3,
            const std::string& idx4
        );

        YetiTensor();

        YetiTensor(const TensorSubsetPtr& subset);

        /**
            Link to a tensor, but with a different scale factor.  This merely
            configures the scale factor, and does not actually scale the tensor
            @param yt
            @param scale
        */
        YetiTensor(
            const YetiTensorPtr& yt,
            double scale
        );

        /**
            Link to a tensor, but with a permutation of the underlying values.
            This just configures the permutation, and does not actually sort the tensor
            @param tensor
            @param perm
        */
        YetiTensor(
            const YetiTensorPtr& tensor,
            const PermutationPtr& perm
        );

        void operator=(const YetiTensorPtr& tensor);

        void operator=(const YetiTensor& tensor);

        /**
            Override of copy constructor.  This should never be used.
            Allow YetiTensorPtr objects should ever be passed around.
            @param tensor
            @throw If called, ever.
        */
        YetiTensor(
            const YetiTensor& tensor
        );

        /**
            Pointer dereference operator for accessing the underlying tensor
        */
        Tensor* operator->() const;

        Tensor* get_tensor() const;

        PermutationPtr get_permutation() const;


        void operator*=(double scale);

        /**
            Add a permutation to the tensor permutation group.
            This should be called before any major tensor operations
            have been performed
            @param p
        */
        void operator+=(const PermutationRuntimeParserPtr& p);

        /**
            Accumulate a contraction of two tensors into
            the given tensor
            @param cxn
        */
        void operator+=(const YetiContractionPtr& cxn);

        /**
            Accumulate a contraction of two tensors into
            the given tensor
            @param cxn
        */
        void operator-=(const YetiContractionPtr& cxn);

        /**
            Subtract the values from the parameter tensor from this tensor
            The parameter tensor may be permuted and scale during
            the accumulation.
            @param tensor
        */
        void operator-=(const YetiTensorPtr& tensor);

        /**
            Add the values from the parameter tensor to this tensor.
            The parameter tensor may be permuted and scale during
            the accumulation.
            @param tensor
        */
        void operator+=(const YetiTensorPtr& tensor);

        /**
            Configure the string indices to be used for parsing
            tensor operations
            @param str
            @return A yeti tensor configured with the parameter string
        */
        YetiTensorPtr operator()(const std::string& str) const;

        MatrixPtr matrix(const std::string& rowstr, const std::string& colstr) const;

        YetiTensorPtr operator()(
            const std::string& idx1,
            const std::string& idx2
        ) const;

        YetiTensorPtr operator()(
            const std::string& idx1,
            const std::string& idx2,
            const std::string& idx3,
            const std::string& idx4
        ) const;

        /**
            Accumulate a contraction of two tensors into
            the given tensor
            @param cxn
        */
        void operator=(const YetiContractionPtr& cxn);

        /**
            Parse the parameter string, determine the permutation specified,
            and sort the underlying tensor
            @param str
        */
        void sort(const std::string& str);

        void sort(
            const std::string &idx1,
            const std::string &idx2,
            const std::string &idx3,
            const std::string &idx4
        );

        void scale(double scale);

        void operator<<(const IndexSubsetPtr& subset);

        void free();
};

class YetiTensorPtr :
    public boost::intrusive_ptr<YetiTensor>
{

    public:
        YetiTensorPtr(YetiTensor* t);

        YetiTensorPtr();

        void operator=(const YetiContractionPtr& cxn);

        void operator+=(const YetiContractionPtr& cxn);

        void operator-=(const YetiContractionPtr& cxn);
};

class YetiSubsetTensor :
    public YetiTensor
{
    private:
        friend class YetiTensor;

        YetiTensorPtr tensor_;

    public:
        YetiSubsetTensor(const YetiTensorPtr& tensor);

};

class PermutationRuntimeParser :
    public StringParser
{

    private:
        short sign_;

    public:
        PermutationRuntimeParser(const std::string& str);

        short sign() const;

        PermutationPtr get_permutation(
            const StringParserPtr& parser
        );

};


class YetiContraction :
    public StringParser
{

    public:
        YetiTensorPtr ltensor;

        YetiTensorPtr rtensor;

        MatrixIndexPtr lindex;

        MatrixIndexPtr rindex;

        PermutationGroupPtr required_grp;

        usi nidx_target;

        usi nidx_cxn;

        double scale;

        YetiContraction(const std::string& cxnstring);

};

class IndexSubset :
    public StringParser
{
    public:
        IndexSubset(const std::string& str);

        IndexRangeTuplePtr get_tuple() const;

};

IndexSubsetPtr
subset(const std::string& str);

PermutationRuntimeParserPtr
P(const std::string& str);

PermutationRuntimeParserPtr
S(const std::string& str);

YetiContractionPtr
operator*(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor);

YetiContractionPtr
operator&(const PermutationRuntimeParserPtr& p, const YetiContractionPtr& cxn);

YetiTensorPtr
operator->*(const PermutationRuntimeParserPtr& p, const YetiTensorPtr& tensor);

YetiContractionPtr
operator*(double scale, const YetiContractionPtr& cxn);

YetiTensorPtr
operator*(double scale, const YetiTensorPtr& yt);

MatrixIndexPtr
get_matrix_index(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor, std::string& cxnstring);


template <class T, class U, bool flag> class InvalidOperatorTimes {
    public: InvalidOperatorTimes() { FailCompileTime<flag> t; }
};

class InvalidContractionOperator_UseDereferenceArrow_NotStar {};

template<bool flag> class InvalidOperatorTimes<YetiTensorPtr,YetiTensorPtr,flag> {
    public: InvalidOperatorTimes() {
    FailCompileTime<flag,InvalidContractionOperator_UseDereferenceArrow_NotStar> t; }
};

template <class T, class U>
void
operator*(const T& t, const U& u)
{
    InvalidOperatorTimes<T,U,true> fail;
}

std::ostream&
operator<<(std::ostream& os, const YetiTensor& tensor);

}

#endif // Header Guard
