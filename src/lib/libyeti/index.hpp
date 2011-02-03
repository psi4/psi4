#ifndef yeti_index_hpp
#define yeti_index_hpp

#include "tuple.h"

namespace yeti {

class IndexRange;
class IndexDescr;
class Indexer;
class IndexRestriction;
class IndexList;
class IndexMap;
class IndexSet;
class IndexRangeTuple;
class IndexRangeLocation;
class IndexRangeLocationCompare;

typedef boost::intrusive_ptr<IndexRange> IndexRangePtr;
typedef boost::intrusive_ptr<IndexDescr> IndexDescrPtr;
typedef boost::intrusive_ptr<const IndexDescr> constIndexDescrPtr;
typedef boost::intrusive_ptr<IndexRestriction> IndexRestrictionPtr;
typedef boost::intrusive_ptr<const IndexRestriction> constIndexRestrictionPtr;
typedef boost::intrusive_ptr<IndexList> IndexListPtr;
typedef boost::intrusive_ptr<const IndexList> constIndexListPtr;
typedef boost::intrusive_ptr<Indexer> IndexerPtr;

typedef boost::intrusive_ptr<IndexMap> IndexMapPtr;
typedef boost::intrusive_ptr<const IndexMap> constIndexMapPtr;
typedef boost::intrusive_ptr<IndexSet> IndexSetPtr;
typedef boost::intrusive_ptr<const IndexSet> constIndexSetPtr;


typedef boost::intrusive_ptr<IndexRangeTuple> IndexRangeTuplePtr;

typedef boost::intrusive_ptr<IndexRangeLocation> IndexRangeLocationPtr;

//typedef Tuple<IndexRange*> SubindexTuple;
class SubindexTuple;
typedef boost::intrusive_ptr<SubindexTuple> SubindexTuplePtr;

typedef Tuple<IndexDescrPtr> IndexDescrTuple;
typedef boost::intrusive_ptr<IndexDescrTuple> IndexDescrTuplePtr;

}

#endif

