#include "permutation.h"
#include "sort.h"
#include "env.h"
#include "index.h"
#include "exception.h"
#include "malloc.h"
#include "filler.h"
#include "tile.h"
#include "data.h"
#include "runtime.h"

using namespace yeti;
using namespace std;

UnitEstimaterPtr TileEstimater::unit_ = new UnitEstimater;

ThreadedTileElementComputer::ThreadedTileElementComputer(const TileElementComputerPtr& comp)
    :
    mindepth_(0),
    fillers_(YetiRuntime::nthread(), 0) //create empty vector
{
    fillers_[0] = comp;
    for (uli i=1; i < YetiRuntime::nthread(); ++i)
        fillers_[i] = comp->copy();
}

ThreadedTileElementComputer::~ThreadedTileElementComputer()
{
}

TileEstimater*
ThreadedTileElementComputer::get_estimater(usi depth) const
{
    TileElementComputer* comp = fillers_[0].get();
    if (comp)
        return comp->get_estimater(depth);
    else
        return 0;
}

void
ThreadedTileElementComputer::compute(Tile *tile, Data *data, uli threadnum)
{
    data->compute(tile, fillers_[threadnum].get());
}

bool
ThreadedTileElementComputer::equals(Tile* tile, Data* data)
{
    return data->equals(tile, fillers_[0].get());
}

TileElementComputer*
ThreadedTileElementComputer::get_thread_computer(uli threadnum) const
{
    return fillers_[threadnum].get();
}

void
ThreadedTileElementComputer::set_mindepth(usi depth)
{
    mindepth_ = depth;
}

usi
ThreadedTileElementComputer::mindepth() const
{
    return mindepth_;
}

TileElementComputer::TileElementComputer()
    : buffer_(0)
{
}

TileElementComputer::~TileElementComputer()
{
    if (buffer_)
        ::free(buffer_);
}

TileEstimater*
TileElementComputer::get_estimater(usi depth) const
{
    return TileEstimater::get_unit_estimater();
}

template <class data_t>
bool
TileElementComputer::equals(Tile* tile, const data_t* data)
{
    if (buffer_ == 0)
    {
        raise(SanityCheckError, "filler equals called without allocating buffer");
    }

    //compute values into the buffer
    data_t* vals = reinterpret_cast<data_t*>(buffer_);
    compute(tile, vals);

    DataBlock* block = tile->get_data();
    uli n = block->n();
    for (size_t i=0; i < n; ++i)
    {
        if (!TestEquals<data_t>::equals(vals[i], data[i]))
        {
            cerr << vals[i] << " != " << data[i] << endl;
            cout << tile << endl;
            tile->get_data()->data()->reference(buffer_);
            cout << tile << endl;
            abort();
            return false;
        }
    }

    return true;
}

bool
TileElementComputer::equals(Tile* tile, const double* data)
{
    return this->equals<double>(tile, data);
}

bool
TileElementComputer::equals(Tile* tile, const float* data)
{
    return this->equals<float>(tile, data);
}

bool
TileElementComputer::equals(Tile* tile, const int* data)
{
    return this->equals<int>(tile, data);
}

bool
TileElementComputer::equals(Tile* tile, const quad* data)
{
    return this->equals<quad>(tile, data);
}

void
TileElementComputer::compute(Tile* tile, double *data)
{
    raise(SanityCheckError, "tile filler not configured to compute type double");
}

void
TileElementComputer::compute(Tile* tile, quad* data)
{
    raise(SanityCheckError, "tile filler not configured to compute type quad");
}

void
TileElementComputer::compute(Tile* tile, int *data)
{
    raise(SanityCheckError, "tile filler not configured to compute type int");
}

void
TileElementComputer::compute(Tile* tile, float *data)
{
    raise(SanityCheckError, "tile filler not configured to compute type float");
}

void
TileElementComputer::allocate_buffer(uli maxblocksize)
{
    if (buffer_) //clear the old one
        ::free(buffer_);
    buffer_ = ::malloc(maxblocksize);
}

UnitEstimater*
TileEstimater::get_unit_estimater()
{
    return unit_.get();
}

float
UnitEstimater::max_log(const uli *indices) const
{
    return 0;
}

void
MemsetElementComputer::compute(Tile *tile, double *data)
{
    ::memset(data, 0, tile->get_data()->data()->size());
}

void
MemsetElementComputer::compute(Tile *tile, int *data)
{
    ::memset(data, 0, tile->get_data()->data()->size());
}

void
MemsetElementComputer::compute(Tile *tile, float *data)
{
    ::memset(data, 0, tile->get_data()->data()->size());
}

void
MemsetElementComputer::compute(Tile *tile, quad *data)
{
    ::memset(data, 0, tile->get_data()->data()->size());
}

TileElementComputer*
MemsetElementComputer::copy() const
{
    return new MemsetElementComputer;
}

void
NullElementComputer::compute(Tile *tile, double *data)
{
}

void
NullElementComputer::compute(Tile *tile, int *data)
{
}

void
NullElementComputer::compute(Tile *tile, float *data)
{
}

void
NullElementComputer::compute(Tile *tile, quad *data)
{
}

TileElementComputer*
NullElementComputer::copy() const
{
    return new NullElementComputer;
}

