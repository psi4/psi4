#ifndef yeti_tensor_hpp
#define yeti_tensor_hpp

namespace yeti {


class Tensor;
class Scalar;
class TensorConfiguration;

typedef boost::intrusive_ptr<Tensor> TensorPtr;
typedef boost::intrusive_ptr<const Tensor> constTensorPtr;

typedef boost::intrusive_ptr<Scalar> ScalarPtr;

typedef boost::intrusive_ptr<TensorConfiguration> TensorConfigurationPtr;

bool
tensor_less(const TensorPtr& l, const TensorPtr& r);

}

#endif

