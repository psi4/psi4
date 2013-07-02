#ifndef yeti_tensor_hpp
#define yeti_tensor_hpp

namespace yeti {


class Tensor;
class TensorConfiguration;

typedef boost::intrusive_ptr<Tensor> TensorPtr;
typedef boost::intrusive_ptr<TensorConfiguration> TensorConfigurationPtr;

/**
    This specifies a location for a TensorDataController to start
    reading in its data block from a parent tensor
*/
typedef size_t tensor_vir_addr_unsigned_offset_t;

bool
tensor_less(Tensor* l, Tensor* r);


}

#endif

