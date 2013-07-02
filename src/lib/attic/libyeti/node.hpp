#ifndef yeti_node_hpp
#define yeti_node_hpp

namespace yeti {


class TileNode;
class DataNode;
class MetaDataNode;

typedef boost::intrusive_ptr<TileNode> TileNodePtr;
typedef boost::intrusive_ptr<DataNode> DataNodePtr;
typedef boost::intrusive_ptr<MetaDataNode> MetaDataNodePtr;

typedef NodeMap<TileNode> TileMap;

/**
    This specifies a location relative to the beginning of a data block
    for a data node to start reading.
*/
typedef size_t data_block_vir_addr_unsigned_offset_t;

}

#endif

