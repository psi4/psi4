#ifndef yeti_contraction_hpp
#define yeti_contraction_hpp

namespace yeti {

class Contraction;
class ContractionTask;
class ContractionQueue;
class ContractionEngine;

typedef boost::intrusive_ptr<Contraction> ContractionPtr;
typedef boost::intrusive_ptr<const Contraction> constContractionPtr;
typedef boost::intrusive_ptr<ContractionTask> ContractionTaskPtr;
typedef boost::intrusive_ptr<const ContractionTask> constContractionTaskPtr;

typedef boost::intrusive_ptr<ContractionQueue> ContractionQueuePtr;


}

#endif
