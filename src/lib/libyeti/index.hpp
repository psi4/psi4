#ifndef yeti_index_hpp
#define yeti_index_hpp

#include "tuple.h"

namespace yeti {

class IndexRange;
class IndexDescr;
class Indexer;
class IndexList;
class IndexRangeTuple;
class IndexRangeLocation;
class IndexRangeLocationCompare;

typedef boost::intrusive_ptr<IndexRange> IndexRangePtr;
typedef boost::intrusive_ptr<IndexDescr> IndexDescrPtr;
typedef boost::intrusive_ptr<const IndexDescr> constIndexDescrPtr;
typedef boost::intrusive_ptr<Indexer> IndexerPtr;


typedef boost::intrusive_ptr<IndexRangeTuple> IndexRangeTuplePtr;

typedef boost::intrusive_ptr<IndexRangeLocation> IndexRangeLocationPtr;

class SubindexTuple;
typedef boost::intrusive_ptr<SubindexTuple> SubindexTuplePtr;

typedef Tuple<IndexDescrPtr> IndexDescrTuple;
typedef boost::intrusive_ptr<IndexDescrTuple> IndexDescrTuplePtr;

}

#endif

