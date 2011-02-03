#ifndef yeti_tensor_parser_hpp
#define yeti_tensor_parser_hpp

namespace yeti {

class PermutationRuntimeParser;
class YetiContraction;
class YetiTensor;
class YetiTensorPtr;
class StringParser;
class TensorSubset;
class IndexSubset;

typedef boost::intrusive_ptr<PermutationRuntimeParser> PermutationRuntimeParserPtr;
//typedef boost::intrusive_ptr<YetiContraction> YetiContractionPtr;
//typedef boost::intrusive_ptr<YetiTensor> YetiTensorPtr;
typedef boost::intrusive_ptr<StringParser> StringParserPtr;

typedef boost::intrusive_ptr<TensorSubset> TensorSubsetPtr;
typedef boost::intrusive_ptr<IndexSubset> IndexSubsetPtr;

}

#endif

