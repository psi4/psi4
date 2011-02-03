#ifndef yeti_tile_hpp
#define yeti_tile_hpp

namespace yeti {

class TileMap;
class Tile;
class TileLocation;
class TileDistributer;
class TileIterator;
class TileMapBuilder;
class TileRegistry;
class TileFilter;

typedef boost::intrusive_ptr<Tile> TilePtr;
typedef boost::intrusive_ptr<TileMap> TileMapPtr;
typedef boost::intrusive_ptr<TileLocation> TileLocationPtr;
typedef boost::intrusive_ptr<TileMapBuilder> TileMapBuilderPtr;
typedef boost::intrusive_ptr<TileRegistry> TileRegistryPtr;
typedef boost::intrusive_ptr<TileFilter> TileFilterPtr;

typedef boost::intrusive_ptr<TileDistributer> TileDistributerPtr;
typedef boost::intrusive_ptr<TileIterator> TileIteratorPtr;


}

#endif

