#ifndef yeti_index_hpp
#define yeti_index_hpp

#include "tuple.h"

namespace yeti {

class IndexRange;
class IndexDescr;
class TensorIndexDescr;
class Indexer;

typedef boost::intrusive_ptr<IndexRange> IndexRangePtr;
typedef boost::intrusive_ptr<IndexDescr> IndexDescrPtr;
typedef boost::intrusive_ptr<TensorIndexDescr> TensorIndexDescrPtr;
typedef boost::intrusive_ptr<Indexer> IndexerPtr;

class SubindexTuple;
typedef boost::intrusive_ptr<SubindexTuple> SubindexTuplePtr;

typedef Tuple<IndexDescrPtr> IndexDescrTuple;
typedef boost::intrusive_ptr<IndexDescrTuple> IndexDescrTuplePtr;

}

#endif

