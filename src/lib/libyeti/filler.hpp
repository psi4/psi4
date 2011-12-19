#ifndef yeti_filler_hpp
#define yeti_filler_hpp

namespace yeti {

class TensorElementComputer;
class ThreadedTensorElementComputer;
class TensorValueEstimater;
class UnitEstimater;
class TensorElementFilter;

typedef boost::intrusive_ptr<TensorValueEstimater> TensorValueEstimaterPtr;
typedef boost::intrusive_ptr<UnitEstimater> UnitEstimaterPtr;
typedef boost::intrusive_ptr<ThreadedTensorElementComputer> ThreadedTensorElementComputerPtr;
typedef boost::intrusive_ptr<TensorElementComputer> TensorElementComputerPtr;
typedef boost::intrusive_ptr<TensorElementFilter> TensorElementFilterPtr;

class TwoElectronEstimableComputer;
class TEIShellComputeFunctor;
typedef boost::intrusive_ptr<TEIShellComputeFunctor> TEIShellComputeFunctorPtr;

}

#endif

