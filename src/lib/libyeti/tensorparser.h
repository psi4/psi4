#ifndef _yeti_tensor_parser_h
#define _yeti_tensor_parser_h

#include <iostream>

#include "class.h"
#include "tensor.hpp"
#include "matrix.hpp"
#include "tensorparser.hpp"
#include "permutation.hpp"
#include "index.hpp"
#include "filler.hpp"
#include "matrix.hpp"
#include "elementop.hpp"
#include "gigmatrix.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

/**
    @class StringParser
    Class for parsing/tokenizing strings
*/
class StringParser :
    public smartptr::Countable
{

    protected:
        std::string str_;

        std::vector<std::string> entries_;

        usi get_nindex_front(const std::string& str);

    public:
        typedef std::vector<std::string>::const_iterator iterator;

        StringParser(const std::string& str);

        virtual ~StringParser(){};

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

        Permutation* get_permutation(const PermutationRuntimeParserPtr& perm);
};

class YetiTensor :
    public StringParser
{
    private:
        TensorPtr tensor_;

        /**
            Permutation of the tensor values that will eventually be used
            in a contraction
        */
        Permutation* perm_;

        /**
            Scale factor the tensor that will eventually be used in a contraction
        */
        double scale_;

        void null_check() const;

        void add(const PermutationRuntimeParserPtr& perm);

    public:
        /**
            @name
            @str The list of indices, e.g. T("i,j,a,b") comma separated
                 that define the tensor
        */
        YetiTensor(
            const std::string& name,
            const std::string& str,
            const PermutationRuntimeParserPtr& p1 = 0,
            const PermutationRuntimeParserPtr& p2 = 0,
            const PermutationRuntimeParserPtr& p3 = 0,
            const PermutationRuntimeParserPtr& p4 = 0,
            const PermutationRuntimeParserPtr& p5 = 0
        );

        YetiTensor(
            const std::string& name,
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC
        );

        YetiTensor(
            const std::string& name,
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC,
            const std::string& idxD
        );

        YetiTensor(
            const TensorPtr& tensor,
            Permutation* p,
            double scale,
            const std::string& str
        );

        YetiTensor();

        ~YetiTensor();

        YetiTensor(const YetiTensorPtr& yt);

        void init(
            const std::string& name,
            const std::string& a1,
            const std::string& a2,
            const std::string& a3
        );

        void init(
            const std::string& name,
            const std::string& a1,
            const std::string& a2,
            const std::string& a3,
            const std::string& a4
        );

        void init(
            const std::string& name,
            const std::string& str,
            const PermutationRuntimeParserPtr& p1 = 0,
            const PermutationRuntimeParserPtr& p2 = 0,
            const PermutationRuntimeParserPtr& p3 = 0,
            const PermutationRuntimeParserPtr& p4 = 0,
            const PermutationRuntimeParserPtr& p5 = 0
        );

        void accumulate(const YetiTensorPtr& tensor, double scale);

        void accumulate(const YetiContractionPtr& cxn, double scale);

        bool is_initialized() const;

        bool equals(const YetiTensorPtr& yt) const;

        void free();

        void clear();

        /**
            Pointer dereference operator for accessing the underlying tensor
        */
        Tensor* operator->() const;

        YetiTensorPtr operator()(const std::string& str) const;

        YetiTensorPtr operator()(
            const std::string& idxA,
            const std::string& idxB
        ) const;

        YetiTensorPtr operator()(
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC
        ) const;

        YetiTensorPtr operator()(
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC,
            const std::string& idxD
        ) const;

        Permutation* get_permutation() const;

        TensorPtr& operator*();

        const TensorPtr& operator*() const;

        double get_scale() const;

        Tensor* get_tensor() const;

        double norm();

        bool nonzero() const;

        bool unique_nonzero() const;

        void set_scale(double scale);

        void set_permutation(Permutation* perm);

        void fill(TensorElementComputer* filler);

        void fill(const MatrixPtr& matrix);

        bool equals(const void* vals);

        void distribute(const std::string& str);

        void sort(const std::string& str);

        void sort(
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC
        );

        void sort(
            const std::string& idxA,
            const std::string& idxB,
            const std::string& idxC,
            const std::string& idxD
        );

        void zero();

        bool is_null() const;

        bool is_nonnull() const;

        void element_op(ElementOp* op);
        
        void print(std::ostream& os = std::cout);

        void operator=(const YetiTensorPtr& yt);

        void operator=(const YetiTensor& yt);

        YetiTensorPtr copy() const;

};

class YetiTensorPtr :
    public boost::intrusive_ptr<YetiTensor>
{

    public:
        YetiTensorPtr(YetiTensor* t);

        YetiTensorPtr();

        void operator+=(const YetiContractionPtr& cxn);

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

        void operator*=(double scale);

        void operator<<=(const YetiTensorPtr& tensor);

        void operator|=(const YetiTensorPtr& tensor);

        void operator|=(const YetiContractionPtr& tensor);

        void operator|=(ElementOp* op);

        YetiTensor* operator->() const;

        YetiTensorPtr& operator[](const std::string& str);

        void assign(const YetiTensorPtr& parser);

        void operator+=(const PermutationRuntimeParserPtr& parser);


};

class YetiMatrix :
    public smartptr::Countable
{

    private:
        YetiTensorPtr tensor_;

        StringParser rowparser_;

        StringParser colparser_;

        MatrixConfigurationPtr config_;

    public:
        YetiMatrix(
            const YetiTensor& yt,
            const std::string& rows,
            const std::string& cols
        );

        MatrixConfiguration* get_config() const;

        void get_matrix(RectMatrixPtr& matrix);

        void get_matrix(SymmMatrixPtr& matrix);

        void accumulate(RectMatrixPtr& matrix);

        void accumulate(SymmMatrixPtr& matrix);
};

class YetiContraction :
    public StringParser
{

    public:
        YetiTensorPtr ltensor;

        YetiTensorPtr rtensor;

        MatrixIndexPtr lindex;

        MatrixIndexPtr rindex;

        PermutationSetPtr symmetrization_set;

        Permutation* left_presort;

        Permutation* right_presort;

        double contraction_scale;

        double post_scale;

        YetiContraction(
            const std::string& cxnstring,
            const YetiTensorPtr& ltensor,
            const YetiTensorPtr& rtensor
        );

        ~YetiContraction();

};

class YetiContractionPtr :
    public boost::intrusive_ptr<YetiContraction>
{
    private:
        template <typename data_t>
        data_t
        dot_product() const;

    public:
        YetiContractionPtr(YetiContraction* cxn);

        void operator=(const YetiTensorPtr& tensor);

        operator double() const;

        operator quad() const;

        operator int() const;
};


class PermutationRuntimeParser :
    public StringParser
{

    private:
        short sign_;

    public:
        PermutationRuntimeParser(const std::string& str);

        short sign() const;

};


PermutationRuntimeParserPtr
Permute(const std::string& str);

PermutationRuntimeParserPtr
Permute(
    const std::string& idxA,
    const std::string& idxB,
    const std::string& idxC,
    const std::string& idxD
);

PermutationRuntimeParserPtr
SymmPermute(const std::string& str);

void
get_matrix_index(
    const YetiTensorPtr& ltensor,
    const YetiTensorPtr& rtensor,
    std::string& cxnstring,
    MatrixIndexPtr& lindex,
    MatrixIndexPtr& rindex,
    Permutation*& left_presort,
    Permutation*& right_presort
);

double&
operator+=(double& e, const YetiContractionPtr& cxn);

double&
operator-=(double& e, const YetiContractionPtr& cxn);


YetiTensorPtr
operator->*(const PermutationRuntimeParserPtr& p, const YetiTensorPtr& tensor);

YetiContractionPtr
operator*(double scale, const YetiContractionPtr& cxn);

YetiContractionPtr
operator&&(double scale, const YetiContractionPtr& cxn);

YetiTensorPtr
operator*(double scale, const YetiTensorPtr& yt);

YetiContractionPtr
operator&(const PermutationRuntimeParserPtr& p, const YetiContractionPtr& cxn);

YetiContractionPtr
operator&(const std::list<PermutationRuntimeParserPtr>& p, const YetiContractionPtr& cxn);

std::list<PermutationRuntimeParserPtr>&
operator&(const PermutationRuntimeParserPtr& p,const PermutationRuntimeParserPtr& q);

YetiContractionPtr
operator*(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor);

YetiContractionPtr
operator*(const YetiTensorPtr& ltensor, const YetiContractionPtr& rtensor);

YetiContractionPtr
operator*(const YetiContractionPtr& ltensor, const YetiTensorPtr& rtensor);

YetiTensorPtr
operator-(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor);

YetiTensorPtr
operator+(const YetiTensorPtr& ltensor, const YetiTensorPtr& rtensor);

std::ostream&
operator<<(std::ostream& os, const YetiTensorPtr& yt);

std::ostream&
operator<<(std::ostream& os, const YetiContractionPtr& yt);



}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // Header Guard
