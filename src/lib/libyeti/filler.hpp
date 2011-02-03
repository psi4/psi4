#ifndef yeti_filler_hpp
#define yeti_filler_hpp

namespace yeti {

class TileElementComputer;
class ThreadedTileElementComputer;
class TileEstimater;
class UnitEstimater;

typedef boost::intrusive_ptr<TileEstimater> TileEstimaterPtr;
typedef boost::intrusive_ptr<UnitEstimater> UnitEstimaterPtr;
typedef boost::intrusive_ptr<ThreadedTileElementComputer> ThreadedTileElementComputerPtr;
typedef boost::intrusive_ptr<TileElementComputer> TileElementComputerPtr;

}

#endif

