#ifndef yeti_tensor_parser_hpp
#define yeti_tensor_parser_hpp

namespace yeti {

class PermutationRuntimeParser;
class YetiContraction;
class YetiContractionPtr;
class YetiTensor;
class YetiTensorPtr;
class StringParser;
class TensorSubset;

typedef boost::intrusive_ptr<PermutationRuntimeParser> PermutationRuntimeParserPtr;
typedef boost::intrusive_ptr<StringParser> StringParserPtr;
typedef boost::intrusive_ptr<TensorSubset> TensorSubsetPtr;

}

#endif

