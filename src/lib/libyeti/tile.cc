#include "tile.h"
#include "data.h"
#include "index.h"
#include "tensor.h"
#include "class.h"
#include "exception.h"
#include "env.h"
#include "permutationimpl.h"
#include "sortimpl.h"
#include "dataimpl.h"
#include "malloc.h"
#include "contraction.h"
#include "matrix.h"
#include "threadimpl.h"
#include "opimpl.h"
#include "threadimpl.h"
#include "elementop.h"
#include "filler.h"

using namespace yeti;
using namespace std;

DECLARE_STATIC_THREAD_WORKSPACE(TileIteratorWorkspace);

DECLARE_MALLOC(Tile);
DECLARE_MALLOC(TileMap);


#define MAX_META_DEPTH 10

TileDistributer::TileDistributer(
    const usi* dist_indices,
    usi ndist_indices
) :
    sorted_task_list_(0),
    ntasks_(0),
    nworkers_(YetiRuntime::nproc()),
    ndistr_index_(ndist_indices),
    distr_indices_(yeti_malloc_perm())
{
    ::memcpy(distr_indices_, dist_indices, ndist_indices * sizeof(usi));
}

TileDistributer::~TileDistributer()
{
    if (sorted_task_list_)
        delete[] sorted_task_list_;
    yeti_free_perm(distr_indices_);
}

void
TileDistributer::tally(const TilePtr &tile)
{
}

void
TileDistributer::sort_tasks()
{
    if (!sorted_task_list_) //allocate
        sorted_task_list_ = new int[ntasks_];

    //create a workspace array for doing the sort
    uli* task_sizes = new uli[ntasks_];


    //now sort the array
    quicksort(task_sizes, sorted_task_list_, ntasks_);

    delete[] task_sizes;
}

PermutationalSymmetryFilter::PermutationalSymmetryFilter(
    const PermutationGroupPtr& pgrp
) :
    pgrp_(pgrp)
{
}

bool
PermutationalSymmetryFilter::is_zero(
    const uli *indices,
    usi depth
) const
{
    //return not unique
    if (pgrp_->order() > 1)
        return !Tile::is_unique(indices, pgrp_, depth);
}

TileFilter*
PermutationalSymmetryFilter::copy() const
{
    return new PermutationalSymmetryFilter(pgrp_);
}

size_t
TileDistributer::process_number(const TilePtr &tile) const
{
    uli index;
    uli tasknum = sorted_task_list_[index];

    //determine the remainder and pass number
    uli remainder = tasknum % nworkers_;
    uli pass = tasknum / nworkers_;

    uli num;
    if (pass % 2) //odd numbered pass, assign tasks starting from the end
    {
        num = nworkers_ - remainder - 1;
    }
    else //even numbered pass, assign tasks starting from the beginning
    {
        num = remainder;
    }

    return num;
}


TileMap::TileMap(
    const IndexRangeTuplePtr& tuple,
    const PermutationGroupPtr& pgrp,
    const TilePtr& parent
)
    : indexer_(0),
    tmpindices_(yeti_malloc_indexset()),
    image_(yeti_malloc_indexset()),
    offsets_(yeti_malloc_indexset()),
    nindices_(yeti_malloc_indexset()),
    //parent_(parent.get()),
    tuple_(tuple),
    pgrp_(pgrp),
    computed_(false),
    nindex_(tuple->size()),
    tiles_(0),
    tmptuple_(reinterpret_cast<IndexRange**>(yeti_malloc_indexptr())),
    filler_(0),
    allocated_(false),
    maxdepth_(tuple->maxdepth()),
    mindepth_(tuple->mindepth()),
    parents_(reinterpret_cast<Tile**>(yeti_malloc_indexptr())),
    index_strides_(yeti_malloc_indexset()),
    is_tiled_index_(yeti_malloc_perm()),
    rigorously_zero_tiles_(0),
    nparents_(1)
{
    parents_[0] = parent.get();

    uli size = 1;

    usi idx = 0;
    const uli* tile_indices = parent->indices();
    usi nactive_indices = 0;
    usi* indexmap = yeti_malloc_perm();

    foreach(range, tuple, IndexRange,
        if (range->depth() < maxdepth_)
        {
            offsets_[idx] = 0;
            nindices_[idx] = 1;
            is_tiled_index_[idx] = 0;
        }
        else
        {
            nindices_[idx] = range->n();
            offsets_[idx] = range->start();
            is_tiled_index_[idx] = 1;
            size *= range->n();
            //this should be indexed by the indexer
            indexmap[nactive_indices] = idx;
            ++nactive_indices;
        }
        ++idx;
    )

    Permutation* nullp = 0;
    indexer_ = new Indexer(nindices_, offsets_, indexmap, nactive_indices, nullp);
    yeti_free_perm(indexmap);

    tiles_ = new CountableArray<Tile>(size);
    rigorously_zero_tiles_ = new char[size];
    ::memset(rigorously_zero_tiles_, 0, size);

    uli total_stride = 1;
    for (int i=nindex_ - 1; i >=0; --i)
    {
        index_strides_[i] = total_stride;
        total_stride *= nindices_[i];
    }
}

TileMap::TileMap(
    Tile* tile,
    Tile* parent
)
    : indexer_(new Indexer),
    tmpindices_(0),
    image_(0),
    offsets_(0),
    nindices_(0),
    is_tiled_index_(0),
    index_strides_(0),
    parents_(reinterpret_cast<Tile**>(yeti_malloc_indexptr())),
    nparents_(1),
    tuple_(0),
    pgrp_(0),
    computed_(true), //nothing to compute
    nindex_(0),
    tiles_(0),
    tmptuple_(0),
    allocated_(false),
    maxdepth_(tile->depth() + 1),
    mindepth_(tile->depth() + 1),
    filler_(0),
    rigorously_zero_tiles_(new char[1])
{
    parents_[0] = parent;
    uli size = 1;
    tiles_ = new CountableArray<Tile>(size);
    tiles_->insert(0, tile);
    rigorously_zero_tiles_[0] = 0;
}


TileMap::TileMap()
    : indexer_(0),
    tmpindices_(0),
    image_(0),
    offsets_(0),
    nindices_(0),
    index_strides_(0),
    is_tiled_index_(0),
    parents_(reinterpret_cast<Tile**>(yeti_malloc_indexptr())),
    nparents_(0),
    tuple_(0),
    pgrp_(0),
    maxdepth_(0),
    mindepth_(0),
    computed_(false),
    nindex_(0),
    tiles_(0),
    tmptuple_(0),
    rigorously_zero_tiles_(0)
{
}

TileMap::~TileMap()
{
    if (tmpindices_)
        yeti_free_indexset(tmpindices_);
    if (image_)
        yeti_free_indexset(image_);
    if (nindices_)
        yeti_free_indexset(nindices_);
    if (offsets_)
        yeti_free_indexset(offsets_);
    if (tmptuple_)
        yeti_free_indexptr(tmptuple_);
    if (parents_)
        yeti_free_indexptr(parents_);
    if (is_tiled_index_)
        yeti_free_perm(is_tiled_index_);
    if (index_strides_)
        yeti_free_indexset(index_strides_);

    delete tiles_;
    delete[] rigorously_zero_tiles_;
}

TileMap::iterator
TileMap::begin() const
{
    return tiles_->begin();
}

TileMap*
TileMap::copy(Tile* parent)
{
    TileMap* copy = new TileMap;
    copy->incref();

    copy->indexer_ = indexer_->copy();
    copy->maxdepth_ = maxdepth_;
    copy->computed_ = true;
    copy->nindex_ = nindex_;
    copy->tuple_ = tuple_;
    copy->parents_[0] = parent;
    copy->pgrp_ = pgrp_;
    copy->allocated_ = allocated_;
    copy->nparents_ = 1;

    copy->tmptuple_ = reinterpret_cast<IndexRange**>(yeti_malloc_indexptr());
    copy->image_ = yeti_malloc_indexset();


    copy->tmpindices_ = yeti_malloc_indexset();
    ::memcpy(copy->tmpindices_, tmpindices_, nindex_ * sizeof(uli));

    copy->offsets_ = yeti_malloc_indexset();
    ::memcpy(copy->offsets_, offsets_, nindex_ * sizeof(uli));

    copy->nindices_ = yeti_malloc_indexset();
    ::memcpy(copy->nindices_, nindices_, nindex_ * sizeof(uli));

    copy->index_strides_ = yeti_malloc_indexset();
    ::memcpy(copy->index_strides_, index_strides_, nindex_ * sizeof(uli));

    copy->is_tiled_index_ = yeti_malloc_perm();
    ::memcpy(copy->is_tiled_index_, is_tiled_index_, nindex_ * sizeof(usi));

    copy->tiles_ = new CountableArray<Tile>(tiles_->size());

    Tile** copyptr = copy->tiles_->begin();
    uli n = 0;
    for (Tile** next = begin(); next != end(); ++next, ++n)
    {
        Tile* tile = *next;
        if (tile == 0)
            continue;

        tile->retrieve(NOT_THREADED);
        Tile* newtile = tile->copy(copy);
        tile->release(NOT_THREADED);
        copy->tiles_->insert(n, newtile);
    }

    copy->rigorously_zero_tiles_ = new char[tiles_->size()];
    ::memcpy(copy->rigorously_zero_tiles_, this->rigorously_zero_tiles_, tiles_->size());

    copy->decref();

    return copy;
}

TileMap::iterator
TileMap::end() const
{
    return tiles_->end();
}

void
TileMap::compute_indices(uli comp_index, uli *index_arr)
{
    uli remainder = comp_index;
    for (usi i=0; i < nindex_; ++i)
    {
        uli stride = index_strides_[i];
        uli index = remainder / stride;
        index_arr[i] = index + offsets_[i];
        remainder -= index * stride;
    }
}

IndexRange*
TileMap::get_subrange(usi index_number, uli index) const
{
    if (is_tiled_index_[index_number])
        return tuple_->get(index_number)->get_subindex(index);
    else //just return the parent index as it is not tiled
        return tuple_->get(index_number);
}

uli
TileMap::index(const uli *indices) const
{
    return indexer_->index(indices);
}

const uli*
TileMap::index_strides() const
{
    return index_strides_;
}

Tile*
TileMap::get(size_t index) const
{
    return tiles_->get(index);
}

Tile*
TileMap::get(const size_t *indices) const
{
    uli index = indexer_->index(indices);
#if YETI_SANITY_CHECK
    if (index >= tiles_->size())
    {
        cerr << "invalid index " << index << " in TileMap get" << endl;
        cerr << indexer_ << endl;
        abort();
    }
#endif
    return tiles_->get(index);
}

bool
TileMap::exists(const size_t *indices) const
{
    return tiles_->get(indexer_->index(indices));
}

TilePtr
TileMap::get(const constIndexSetPtr& indexset) const
{
    return get(indexset->data());
}

IndexRangeTuplePtr
TileMap::get_index_ranges() const
{
    return tuple_;
}

void
TileMap::insert(Tile* tile)
{
    //once a tile gets insert, this must be flagged as computed
    computed_ = true;

    if (!tile) //cannot insert null tile
    {
        raise(SanityCheckError, "cannot insert null tile into tile map");
    }

    uli idx = indexer_->index(tile->indices());
    tiles_->insert(idx, tile);
}

void
TileMap::insert(uli index, Tile* tile)
{
    tiles_->insert(index, tile);
}

const usi*
TileMap::is_tiled_index() const
{
    return is_tiled_index_;
}

void
TileMap::iterate(
    usi index
)
{
    if (index == nindex_)
    {
        insert_new_tile();
        return;
    }

    //fetch the object for estimating tile values
    TileEstimater* estimater = 0;
    if (filler_)
        estimater = filler_->get_estimater(maxdepth_, index + 1);

    uli idx = offsets_[index];
    uli stop = nindices_[index] + idx;
    IndexRange* range = tuple_->get(index);
    //maybe or maybe not recurse on this
    bool use_subindex = range->depth() >= maxdepth_;
    for ( ; idx != stop; ++idx) //loop al index values within the given range
    {
        //set the index number at the given index position
        tmpindices_[index] = idx;

        //set the IndexRange at the given index position
        tmptuple_[index] = use_subindex ? range->get_subindex(idx) : range;


        //if there is no object for estimating the given number of indices
        //and recursion depth, always build the tile
        if (estimater == 0)
        {
            iterate(index + 1);
        }
        //if we have an object for estimating and the estimate is below the cutoff,
        //do not create the next tile
        else if (estimater->max_log(tmpindices_) > YetiRuntime::matrix_multiply_cutoff)
        {
            iterate(index + 1);
        }
    }
}

void
TileMap::insert_new_tile()
{
    std::list<TileFilterPtr>::const_iterator it(config()->filters.begin());
    std::list<TileFilterPtr>::const_iterator stop(config()->filters.end());
    for ( ; it != stop; ++it)
    {
        bool is_zero = (*it)->is_zero(tmpindices_, maxdepth_);
        if (is_zero)
        {
            uli index = indexer_->index(tmpindices_);
            rigorously_zero_tiles_[index] = 1;
            return;
        }
    }

    //this corresponds to a permutationally unique block
    IndexRangeTuplePtr subtuple
        = new IndexRangeTuple(nindex_, tmptuple_);

    //create the tile
    uli* indexarr = yeti_malloc_indexset();
    ::memcpy(indexarr, tmpindices_, nindex_ * sizeof(size_t));

    //make(idxset, IndexSet, indexarr, nindex_);
    Tile* tile(new Tile(subtuple, indexarr, this));
    insert(tile);
}

bool
TileMap::is_rigorously_zero(Tile** tile) const
{
    uli ptrdiff = reinterpret_cast<size_t>(tile) - reinterpret_cast<size_t>(tiles_->begin());
    uli offset = ptrdiff / sizeof(Tile*);
    return rigorously_zero_tiles_[offset];
}

bool
TileMap::is_rigorously_zero(uli offset) const
{
    return rigorously_zero_tiles_[offset];
}

bool
TileMap::depths_aligned() const
{
    return mindepth_ == maxdepth_;
}

usi
TileMap::maxdepth() const
{
    return maxdepth_;
}

const uli*
TileMap::nindices() const
{
    return nindices_;
}

uli
TileMap::nrange(usi index) const
{
    return nindices_[index];
}

const uli*
TileMap::offsets() const
{
    return offsets_;
}

void
TileMap::print(ostream& os) const
{
    if (tiles_->n() == 0)
        return; //nothing to print

    os << Env::indent << "Tile Map at " << (void*) this;
    for (usi i=0; i < nparents_; ++i)
    {
        os << endl << Env::indent << "Parent " << i << ": " << (void*) parents_[i];
    }
    os << endl << Env::indent << "Sizes: " << ClassOutput<const uli*>::str(nindex_, nindices());
    os << endl << Env::indent << "Offsets: " << ClassOutput<const uli*>::str(nindex_, offsets());
    os << endl << indexer_;

    usi idx = 0;
    foreach_nonnull(tile, this, Tile,
        os << endl << Env::indent;
        ++Env::indent;
        tile->retrieve(NOT_THREADED);
        os << idx << " -> "; tile->print(os);
        tile->release(NOT_THREADED);
        --Env::indent;
        ++idx;
    )
}

void
TileMap::retrieve()
{
    if (computed_)
        return;

    //set up the filler object
    filler_ = get_parent_tile()->config()->filler;

    iterate(0);

    computed_ = true;
}

void
TileMap::release()
{
    //indexmap_->clear();
    //tiles_.clear();
    //computed_ = false;
}

size_t
TileMap::ntiles() const
{
    //return the number of nonzero entries in the tile
    return tiles_->n();
}

void
TileMap::set_max_log(float log, usi depth)
{
    for (usi i=0; i < nparents_; ++i)
    {
        parents_[i]->set_max_log(log, depth);
    }
}

Tile*
TileMap::get_parent_tile() const
{
    return *parents_;
}

void
TileMap::remove_parent(Tile* tile)
{
    usi idx = 0;
    for ( ; idx < nparents_; ++idx)
    {
        if (parents_[idx] == tile)
            break;
    }

    for ( ; idx < nparents_ - 1; ++idx)
        parents_[idx] = parents_[idx + 1];

    --nparents_;
}

void
TileMap::add_parent(Tile* tile)
{
    parents_[nparents_] = tile;
    ++nparents_;
}

TensorConfiguration*
TileMap::config() const
{
    return get_parent_tile()->config();
}

uli
TileMap::size() const
{
    return tiles_->size();
}

void
TileMap::allocate()
{
    if (allocated_)
        return;

    foreach_nonnull(tile, this, Tile,
        tile->retrieve(NOT_THREADED);
        tile->allocate();
        tile->release(NOT_THREADED);
    )

    allocated_ = true;
}

bool
TileMap::is_allocated() const
{
    return allocated_;
}

void
TileMap::sort(
    const SortPtr& sort,
    void* buffer
)
{
    sort->configure(nindices_);

    Tile** tiles = begin();
    Tile** tilebuf = reinterpret_cast<Tile**>(buffer);
    sort->sort_noscale<Tile*>(tiles, tilebuf);
    ::memcpy(tiles, tilebuf, size() * sizeof(Tile*));

    indexer_->permute(sort->get_permutation());

    uli* indexbuffer = yeti_malloc_indexset();

    ::memcpy(indexbuffer, nindices_, nindex_ * sizeof(uli));
    sort->get_permutation()->permute(indexbuffer, nindices_);

    ::memcpy(indexbuffer, offsets_, nindex_ * sizeof(uli));
    sort->get_permutation()->permute(indexbuffer, offsets_);

    yeti_free_indexset(indexbuffer);

    foreach_nonnull(tile, this, Tile,
        tile->retrieve(NOT_THREADED);
        tile->_sort(sort, buffer);
        tile->release(NOT_THREADED);
    )
}

void
Tile::init_estimate()
{
    if (is_parent())
    {
        //initialize the estimate
        estimate_.maxlog = 0;
        estimate_.depth = MAX_META_DEPTH;
    }
    else
    {
        //initialize the estimate
        estimate_.maxlog = 0;
        estimate_.depth = 0;
    }
}

void
Tile::init()
{
    if (is_parent())
    {
        //construct the tile map
        if (!tilemap_)
        {
            incref();
            tilemap_ = config_->tile_map_builder->build_map(
                                                    ranges_,
                                                    this
                                                 );
            decref();

        }
    }

    computed_ = true;
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices,
    //const PermutationGroupPtr& pgrp,
    const TensorConfigurationPtr& config
) : ranges_(ranges),
    computed_(false),
    //pgrp_(pgrp),
    indices_(indices),
    parent_(0),
    tilemap_(0),
    owner_process_(0),
    has_owner_process_(false), //top level tile, so no owner,
    data_(0),
    config_(config)
{
    init_estimate();
    yeti_register_new(lock_.get());
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices,
    //const PermutationGroupPtr& pgrp,
    TileMap* parent
) : ranges_(ranges),
    //pgrp_(pgrp),
    computed_(false),
    indices_(indices),
    parent_(parent),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(parent->config())
{
    yeti_register_new(lock_.get());
    init_estimate();
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices
) : ranges_(ranges),
    computed_(false),
    indices_(indices),
    parent_(0),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(0)
{
    init_estimate();
    yeti_register_new(lock_.get());
}

Tile::Tile(
    TensorConfiguration* config,
    const size_t *indices,
    const size_t *sizes,
    usi nidx
) :
    ranges_(0),
    //pgrp_(0),
    parent_(0),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(config),
    indices_(0),
    computed_(false)
{
    size_t* copy = yeti_malloc_indexset();
    ::memcpy(copy, indices, nidx * sizeof(size_t));
    indices_ = copy;

    ranges_ = new IndexRangeTuple(nidx);

    for (usi i=0; i < nidx; ++i)
    {
        size_t start = 0;
        size_t n = sizes[i];
        IndexRange* range(new IndexRange(start, n));
        ranges_->set(i, range);
    }

    estimate_.depth = 0;
    estimate_.maxlog = LOG_ZERO;

    yeti_register_new(lock_.get());
}

Tile::Tile(const TensorConfigurationPtr& config)
    :
    //pgrp_(0),
    ranges_(0),
    parent_(0),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(config),
    indices_(0),
    computed_(false)
{
}

Tile::~Tile()
{
    if (tilemap_)
        tilemap_->remove_parent(this);
    yeti_free_indexset(indices_);
}

void
Tile::allocate()
{
    if (tilemap_)
    {
        tilemap_->retrieve();
        tilemap_->allocate();
        tilemap_->release();
    }
    else
    {
        config_->data_factory->allocate(data_.get());
    }
}

void
Tile::accumulate(
    const TilePtr &tile,
    const SortPtr &sort,
    double scale
)
{
    Permutation* p = sort->get_permutation();
    if (tilemap_)
    {
        Tile** it(tilemap_->begin());
        Tile** stop(tilemap_->end());

        Tile** unsorted_tiles = tile->get_map()->begin();

        Tile** acc;
        if (sort->get_permutation()->is_identity())
        {
            acc = unsorted_tiles;
        }
        else
        {
            TileIteratorWorkspace* workspace = get_workspace<TileIteratorWorkspace>(NOT_THREADED);
            usi depth = tilemap_->maxdepth();
            Tile** sorted_tiles = reinterpret_cast<Tile**>(workspace->buffers[depth]);
            sort->configure(tile->get_map()->nindices());
            sort->sort_noscale<Tile*>(unsorted_tiles, sorted_tiles);
            acc = sorted_tiles;
        }

        for ( ; it != stop; ++it, ++acc)
        {
            Tile* tile_acc = *acc;
            if (!tile_acc) //nothing to accumulate
                continue;

            tile_acc->retrieve(NOT_THREADED);

            Tile* tile_me = *it;
            if (tile_me)
            {
                tile_me->retrieve(NOT_THREADED);
                tile_me->accumulate(tile_acc, sort, scale);
                tile_me->release(NOT_THREADED);
            }
            else if (tilemap_->is_rigorously_zero(it))
            {
                //this tile should never accumulate anything
                continue;
            }
            else
            {
                //copy the given tile and insert it
                tile_me = tile_acc->copy(tilemap_.get());
                tile_me->sort(p);
            }


            tile_acc->release(NOT_THREADED);
        }
    }
    else //data tile
    {
        Data* me = data_->data();
        Data* acc = tile->get_data()->data();
        if (p->is_identity())
        {
            me->accumulate(acc, scale, 0);
        }
        else
        {
            //The acc block is the thing being sorted in the
            //accumulation.  It is therefore very important
            //that here we configure the sort for the acc tile
            //and not for the this tile.
            sort->configure(tile->get_index_ranges());
            me->accumulate(acc, scale, sort);
        }
    }
}

TensorConfiguration*
Tile::config() const
{
    return config_.get();
}

Tile*
Tile::copy(TileMap* parent)
{
    Tile* copy = new Tile(parent->config());
    copy->parent_ = parent;
    copy_to(copy);

    return copy;
}

void
Tile::copy_to(Tile* copy)
{
    copy->ranges_ = ranges_;
    copy->owner_process_ = owner_process_;
    copy->has_owner_process_ = has_owner_process_;

    copy->indices_ = yeti_malloc_indexset();
    ::memcpy(copy->indices_, indices_, nindex() * sizeof(uli));

    if (tilemap_)
    {
        tilemap_->retrieve();
        copy->tilemap_ = tilemap_->copy(copy);
        tilemap_->release();
    }

    if (data_)
    {
        copy->data_ = copy->config_->data_factory->get_block(copy);
        copy->data_->allocate(data_->data()->type());
        copy->config_->data_factory->allocate(copy->data_.get());
        copy->data_->retrieve(NOT_THREADED);
        copy->data_->data()->assign(data_->data());
        copy->data_->release(NOT_THREADED);
    }

    copy->estimate_.depth = estimate_.depth;
    copy->estimate_.maxlog = estimate_.maxlog;

}

uli
Tile::total_data_block_size()
{
    if (has_owner_process_ && owner_process_ != YetiRuntime::me())
        return 0;

    if (tilemap_)
    {
        uli size = 0;
        tilemap_->retrieve();
        foreach_nonnull(tile, tilemap_, Tile,
            size += tile->total_data_block_size();
        )
        tilemap_->release();
        return size;
    }
    else
    {
        return data_->data()->size();
    }
}

template <class T>
void
Tile::data_action(const T &obj)
{
    if (tilemap_)
    {
        uli size = 0;
        tilemap_->retrieve();
        foreach_nonnull(tile, tilemap_, Tile,
            tile->retrieve(NOT_THREADED);
            tile->data_action(obj);
            tile->release(NOT_THREADED);
        )
        tilemap_->release();
    }
    else
    {
        obj.data_action(this, data_);
    }
}

usi
Tile::depth() const
{
    if (tilemap_)
    {
        return tilemap_->maxdepth();
    }
    else
    {
        return 0;
    }
}

void
Tile::distribute(const TileDistributerPtr &distr)
{
    if (!tilemap_)
        raise(SanityCheckError, "cannot distribute on tile without tile map");

    if (config_->map_storage_type == TileMap::temporary)
        return; //no distribution of storage on this

    tilemap_->retrieve();
    if (DISTRIBUTION_DEPTH == tilemap_->maxdepth()) //distribute here
    {
        foreach_nonnull(tile, tilemap_, Tile,
            size_t procnum = distr->process_number(tile);
            tile->set_owner_process(procnum);
        )
    }
    else //pass the distribution on
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->distribute(distr);
        )
    }
    tilemap_->release();
}

bool
Tile::_equals(const char** data)
{
    if (tilemap_)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->retrieve(NOT_THREADED);
            bool eq = tile->_equals(data);
            tile->release(NOT_THREADED);
            if (!eq) //equals failed
                return false;
        )
        return true; //success
    }
    else
    {
        if (data_)
        {
            data_->retrieve(NOT_THREADED);
            bool eq = data_->data()->equals(*data);
            (*data) += data_->data()->size();
            data_->release(NOT_THREADED);
            return eq;
        }
        else //must have data!
            return false;
    }
}

bool
Tile::equals(const ThreadedTileElementComputerPtr& filler)
{
    if (tilemap_)
    {
        foreach_nonnull(
            tile, tilemap_, Tile,

            tile->retrieve(NOT_THREADED);
            bool eq = tile->equals(filler);
            tile->release(NOT_THREADED);
            if (!eq) //equals failed
                return false;

        )
        return true; //success
    }
    else
    {
        if (data_)
        {
            data_->retrieve(NOT_THREADED);
            bool eq = filler->equals(this, data_->data());
            data_->release(NOT_THREADED);
            return eq;
        }
        else //must have data!
            return false;
    }
}

bool
Tile::equals(const void* data)
{
    const char* cd = reinterpret_cast<const char*>(data);
    return _equals(&cd);
}

void
Tile::fill()
{
    if (tilemap_)
    {
        tilemap_->retrieve();
        foreach_nonnull(tile, tilemap_, Tile,

            //I don't need to do this tile
            if ( tile->has_owner_process() && tile->owner_process() != YetiRuntime::me() )
                continue;

            tile->retrieve(NOT_THREADED);
            tile->fill();
            tile->release(NOT_THREADED);
        )
        tilemap_->release();
    }
    else //data tile
    {
        if (config_->data_factory->storage_type() != Data::recomputed
            && parent_->is_allocated()
            && config_->filler)
        {
            //compute the data for the tile
            data_->retrieve(NOT_THREADED);
            config_->filler->compute(this, data_->data(), 0);
            float maxlog = data_->data()->max_log();
            data_->release(NOT_THREADED);

            //setting log from min depth
            estimate_.maxlog = maxlog;
            estimate_.depth = 0;
            parent_->set_max_log(maxlog, 0);
        }
    }
}

DataBlock*
Tile::get_data() const
{
    return data_.get();
}

Tile*
Tile::get(const uli* indices) const
{
    return tilemap_->get(indices);
}

bool
Tile::exists(const uli *indices) const
{
    return tilemap_->exists(indices);
}

IndexRangeTuple*
Tile::get_index_ranges() const
{
    return ranges_.get();
}

const size_t*
Tile::indices() const
{
    return indices_;
}

TileMap*
Tile::get_map() const
{
    return tilemap_.get();
}

const tile_estimate_t&
Tile::get_max_log() const
{
    return this->estimate_;
}

bool
Tile::has_owner_process() const
{
    return has_owner_process_;
}

const uli*
Tile::index_sizes() const
{
    return tilemap_->nindices();
}

const uli*
Tile::index_offsets() const
{
    return tilemap_->offsets();
}

bool
Tile::is_aligned() const
{
    return tilemap_->depths_aligned();
}

bool
Tile::is_allocated() const
{
    if (tilemap_)
        return tilemap_->is_allocated();
    else
        return false; //no tile map, so no allocation
}

bool
Tile::is_equivalent(const TilePtr& tile) const
{
    if (type() != tile->type())
        return false;

    if (get_index_ranges()->size() != tile->get_index_ranges()->size())
        return false;

    {
        //validate that the index ranges are equals
        IndexRangeTuple::iterator itme(get_index_ranges()->begin());
        IndexRangeTuple::iterator itother(tile->get_index_ranges()->begin());
        for ( ; itme != get_index_ranges()->end(); ++itme, ++itother)
        {
            bool check = (*itme)->equals(*itother);
            if (!check)
                return false;
       }
    }
    {
        TileMapPtr my_map(get_map());
        TileMapPtr tile_map(tile->get_map());

        bool bool_me = my_map;
        bool bool_tile = tile_map;

        if (bool_me != bool_tile)
        {
            raise(SanityCheckError, "tiles of the same type must both have/have not tile maps");
        }

        if (!bool_me) //no tile map, don't bother continuing the check
            return true;

        //check that the tensors have the same number of tiles
        if (get_map()->ntiles() != tile->get_map()->ntiles())
            return false;

        //now iterate the tile map to determine if tiles are equivlant
        TileMap::iterator itme(tilemap_->begin());
        TileMap::iterator stop(tilemap_->end());
        TileMap::iterator itother(tile->get_map()->begin());
        for ( ; itme != stop; ++itme, ++itother)
        {
            Tile* tile_me = *itme;
            Tile* tile_other = *itother;

            bool bool_me = tile_me;
            bool bool_other = tile_other;
            if (bool_me != bool_other) //misaligned tiles
                return false;

            if (!bool_me) //null tiles
                continue;

            bool check =
                tile_me->is_equivalent(tile_other);

            if (!check) //tiles are not equivalent
                return false;
        }
    }
    return true;
}

bool
Tile::is_grand_parent() const
{
    IndexRange* child(ranges_->get(0)->get_first_child());
    return (child && child->is_parent());
}

bool
Tile::is_parent() const
{
    if (ranges_->size() == 0)
        return false;

    return ranges_->get(0)->is_parent();
}

bool
Tile::is_top_tile() const
{
    //if there is no parent, this is the top tile
    return (!parent_);
}

bool
Tile::is_unique(
    const size_t *indices,
    const PermutationGroupPtr& grp,
    usi depth
)
{
    if (depth == 1) //data block
    {
        return !grp->improves_sort(indices);
    }

    PermutationPtr p(grp->get_lowest_permutation(indices));
    if (p->is_identity())
        return true;

    usi nindex = grp->nindex();
    for (usi i=0; i < nindex; ++i)
    {
        usi idximg = p->image(i);
        if (idximg == i) //should not be factored in
            continue;

        uli obj = indices[i];
        uli img = indices[idximg];

        if (img < obj)
        {
            //there is a permutation that maps the indices
            //onto a "lower" or "canonical" index
            return false; //throw away
        }
        else
        {
           return true;
        }
    }

    return true; //keep, I guess
}

uli
Tile::max_blocksize() const
{
    if (tilemap_)
    {
        uli max = tilemap_->size() * sizeof(void*);
        foreach_nonnull(tile, tilemap_, Tile,
            uli size = tile->max_blocksize();
            if (size > max)
                max = size;
        )
        return max;
    }
    else
    {
        return data_->data()->size();
    }
}

float
Tile::max_log() const
{
    return estimate_.maxlog;
}

uli
Tile::ndata() const
{
    if (tilemap_)
        raise(SanityCheckError, "called ndata on non-data tile");

    size_t ntot = 1;
    foreach(range, ranges_, IndexRange,
        ntot *= range->n();
    )

    return ntot;
}

usi
Tile::nindex() const
{
    ranges_->size();
}

uli
Tile::nrange(usi index) const
{
    if (tilemap_)
        return tilemap_->nrange(index);
    else
        return ranges_->get(index)->n();
}

uli
Tile::ntiles(usi depth)
{
    if (!tilemap_)
        raise(SanityCheckError, "You cannot call ntiles without a tile map");

    tilemap_->retrieve();
    if (depth == tilemap_->maxdepth()) //this is the desired level
    {
        tilemap_->release();
        return tilemap_->ntiles();
    }

    //not yet the desired level
    //accumulate all of the tile sizes from the lower level
    size_t ntot = 0;
    foreach_nonnull(tile, tilemap_, Tile,
        ntot += tile->ntiles(depth);
    )

    tilemap_->release();
    return ntot;
}

uli
Tile::ntiles_max() const
{
    size_t nmax = 1;
    IndexRangeTuple::iterator it(ranges_->begin());
    IndexRangeTuple::iterator stop(ranges_->end());
    for ( ; it != stop; ++it)
    {
        nmax *= (*it)->n();
    }
    return nmax;
}

uli
Tile::ntiles_owned()
{
    return ntiles_owned(DISTRIBUTION_DEPTH);
}

uli
Tile::ntiles_owned(usi depth)
{
    if (!tilemap_)
        raise(SanityCheckError, "You cannot call ntiles_owned without a tile map");

    tilemap_->retrieve();
    if (depth == tilemap_->maxdepth()) //this is the desired level
    {
        size_t ntot = 0;
        foreach_nonnull(tile, tilemap_, Tile,
            if (tile->has_owner_process())
            {
                //check to make sure owner is me
                if ( tile->owner_process() == YetiRuntime::me() )
                    ++ntot;
            }
            else //no owner - belongs to everyone!
            {
                ++ntot;
            }
        )

        tilemap_->release();
        return ntot;
    }

    //loop the tile map and determine how many I own
    uli ntot = 0;
    foreach_nonnull(tile, tilemap_, Tile,

        if (tile->has_owner_process())
        {
            //check to make sure owner is me
            if ( tile->owner_process() == YetiRuntime::me() )
                ntot += tile->ntiles_owned(depth);
        }
        else //no owner - belongs to everyone!
        {
            ntot += tile->ntiles_owned(depth);
        }
    )

    tilemap_->release();
    return ntot;
}

uli
Tile::owner_process() const
{
    return owner_process_;
}

void
Tile::print(ostream& os) const
{
    os << Env::indent << ClassOutput<const size_t*>::str(nindex(), indices_) << "  " << ranges_;
    if (has_owner_process_)
        os << " Node " << owner_process_;

    os << Env::indent << "  Max: 10^(" << stream_printf("%8.4f", estimate_.maxlog) << ") from depth " << estimate_.depth    ;

    if (tilemap_ && tilemap_->ntiles() != 0)
    {
        ++Env::indent;
        os << endl << tilemap_;
        --Env::indent;
    }
    else if (config_->print_data && data_)
    {
        os << endl;
        data_->retrieve(NOT_THREADED);
        data_->print(os);
        data_->release(NOT_THREADED);
    }
}

void
Tile::_retrieve(uli threadnum)
{
    if (!computed_)
    {
        init();
    }

    if (tilemap_) //metadata tile
    {
        tilemap_->retrieve();
    }
    else //fetch the data block
    {
        if (!data_)
        {
            data_ = config_->data_factory->get_block(this);
        }
    }
}

void
Tile::_release(uli threadnum)
{
    if (tilemap_)
        tilemap_->release();
}

void
Tile::element_op(const ElementOpPtr& op)
{
    if (tilemap_)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->retrieve(NOT_THREADED);
            tile->element_op(op);
            tile->release(NOT_THREADED);
        )
    }
    else
    {
        data_->retrieve(NOT_THREADED);
        op->element_op(ranges_.get(), data_->data());
        data_->release(NOT_THREADED);
    }

    op->update(this);
}

void
Tile::register_tiles(
    const TileRegistryPtr& reg,
    usi depth
)
{
    if (!tilemap_)
    {
        cerr << "register tile called on data tile or tile without map" << endl;
        abort();
    }

    Tile** it = tilemap_->begin();
    Tile** stop = tilemap_->end();

    for ( ; it != stop; ++it)
    {
        Tile* tile = *it;
        if (tile == 0)
            continue;

        tile->retrieve(NOT_THREADED);
        if (tile->is_aligned() && tile->depth() == depth)
        {
            //register a pointer to the tile
            //so we can manipulate the tile map later
            reg->register_tile(it);
        }
        else if (tile->depth() > depth)
        {
            tile->register_tiles(reg, depth);
        }
        tile->release(NOT_THREADED);
    }
}

void
Tile::scale(double factor)
{
    element_op(new ScaleOp<double>(factor));
    scale_max_log(log10(factor));
}

void
Tile::scale(int factor)
{
    element_op(new ScaleOp<int>(factor));
    scale_max_log(log10(factor));
}

void
Tile::scale_max_log(float log_factor)
{
    estimate_.maxlog += log_factor;
}

void
Tile::set_max_log(float max, usi depth)
{
    //this is a better estimate than the one we currently have
    if (depth < estimate_.depth)
    {
        estimate_.maxlog = max;
        estimate_.depth = depth;
        if (parent_)
            parent_->set_max_log(max, depth);
    }
    else if (estimate_.maxlog < max)
    {
        estimate_.maxlog = max;
        if (parent_)
            parent_->set_max_log(max, depth);
    }
}

void
Tile::set_max_log(const tile_estimate_t &est)
{
    set_max_log(est.maxlog, est.depth);
}

void
Tile::set_owner_process(size_t owner)
{
    owner_process_ = owner;
    has_owner_process_ = true;
}

void
Tile::_sort(
    const SortPtr &sort,
    void *buffer
)
{
    if (tilemap_)
    {
        tilemap_->retrieve();
        tilemap_->sort(sort, buffer);
        if (tilemap_->get_index_ranges() != ranges_.get())
            tilemap_->get_index_ranges()->permute(sort->get_permutation());
        tilemap_->release();
    }
    else
    {
        sort->configure(ranges_);
        //data block
        data_->sort(sort, buffer);
    }

    //now sort the index ranges
    IndexRange** rangeptr = ranges_->begin();
    IndexRange** rangebuf = reinterpret_cast<IndexRange**>(buffer);
    sort->get_permutation()->permute(rangeptr, rangebuf);
    ::memcpy(rangeptr, rangebuf, ranges_->size() * sizeof(IndexRange*));

    //and the indices
    size_t* indexbuf = reinterpret_cast<size_t*>(buffer);
    sort->get_permutation()->permute(indices_, indexbuf);
    ::memcpy(indices_, indexbuf, nindex() * sizeof(size_t));
}

void
Tile::sort(const PermutationPtr& p)
{
    SortPtr s(new Sort(p, ranges_));
    void* buffer = malloc(this->max_blocksize());

    _sort(s, buffer);
    free(buffer);
}

void
Tile::tally(const TileDistributerPtr &distr)
{
    if (!tilemap_)
        raise(SanityCheckError, "cannot distribute on tile without data map");

    tilemap_->retrieve();
    if (DISTRIBUTION_DEPTH == tilemap_->maxdepth()) //distribute here
    {
        foreach_nonnull(tile, tilemap_, Tile,
            distr->tally(tile);
        )
    }
    else //pass the tallying on
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->tally(distr);
        )
    }
    tilemap_->release();
}

Tile::tile_t
Tile::type() const
{
    if (tilemap_)
        return Tile::meta_data;
    else if (data_)
        return Tile::data;
    else
        raise(SanityCheckError, "tile has no metadata or data.  what have you done?");
}

void
Tile::update()
{
    if (tilemap_)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->update();
        )
    }
    else if (data_)
    {
        data_->retrieve(NOT_THREADED);
        float max_log = data_->data()->max_log();
        data_->release(NOT_THREADED);

        set_max_log(max_log, 0);
    }
}

TileIterator::TileIterator(
    Tile** link,
    const TilePtr& parent,
    bool local,
    size_t me,
    const PermutationPtr& p

) :
    tile_(0),
    link_(link),
    parent_(parent),
    stops_(0),
    iters_(0),
    me_(me),
    maxdepth_(parent->get_map()->maxdepth() - 1), //off by one problem
    local_(local),
    depth_(0),
    done_(false),
    perm_(p),
    buffer_(0)
{
    stops_ = new TileMap::iterator[maxdepth_];
    iters_ = new TileMap::iterator[maxdepth_];

#if 0
    usi nbuffers = maxdepth_ + 1;
    buffer_ = new Tile**[nbuffers];
    for (usi i=0; i < nbuffers; ++i)
    {
        buffers_[i] = new Tile*[tile->maxsize()]
    }
#endif
}

TileIterator::~TileIterator()
{
    delete[] stops_;
    delete[] iters_;
    free(buffer_);
}

void
TileIterator::start()
{
    depth_ = 0;
    done_ = false;

    stop_ = parent_->get_map()->end();
    iter_ = parent_->get_map()->begin();
    stops_[0] = stop_;
    iters_[0] = iter_;

    --iter_; //trick the iterator
    next();
}

void
TileIterator::next()
{
    ++iter_;

    if (iter_ == stop_)
    {
        --depth_;
        if (depth_ == -1)
        {
            //off the end
            done_ = true;
            return;
        }
        iter_ = iters_[depth_];
        stop_ = stops_[depth_];
        next();
    }

    tile_ = *iter_;

    while (tile_ == 0) //keep going until nonnull
        next();

    if (done_)
        return;

    if (local_ && tile_->has_owner_process() && tile_->owner_process() != me_)
        next(); //move along

    if (done_)
        return;

    //now we can go ahead and actually iterate
    if (depth_ == maxdepth_)
    {
        (*link_) = tile_;
        return;
    }

    //not yet at the final recursion depth
    iters_[depth_] = iter_;

    ++depth_;

    iter_ = tile_->get_map()->begin();
    stop_ = tile_->get_map()->end();
    stops_[depth_] = stop_;
    next();
}

bool
TileIterator::done()
{
    return done_;
}

TileIteratorWorkspace::TileIteratorWorkspace()
    : buffers(0)
{
    usi nlevels = YetiRuntime::max_depth() + 1;
    buffers = new char*[nlevels];

    const uli* sizes = YetiRuntime::max_block_sizes();

    //allocate pointer size for all metadata levels
    for (usi i=1; i < nlevels; ++i)
        buffers[i] = reinterpret_cast<char*>( ::malloc(sizes[i] * sizeof(void*)) );

    buffers[0] = reinterpret_cast<char*>( ::malloc(sizes[0] * sizeof(double)) );
}

TileIteratorWorkspace::~TileIteratorWorkspace()
{
    usi nlevels = YetiRuntime::max_depth() + 1;
    for (usi i=0; i < nlevels; ++i)
        ::free(buffers[i]);

    delete[] buffers;
}

DefaultTileMapBuilder::DefaultTileMapBuilder(const PermutationGroupPtr &grp)
    : pgrp_(grp)
{

}

TileMapPtr
DefaultTileMapBuilder::build_map(
    const IndexRangeTuplePtr& tuple,
    const TilePtr& parent
)
{
    return new TileMap(tuple, pgrp_, parent);
}

TileRegistry::TileRegistry()
{
}

void
TileRegistry::print(std::ostream &os) const
{
    os << "Tile Registry";
    node_map::const_iterator it(tilemap_.begin());
    node_map::const_iterator stop(tilemap_.end());
    ++Env::indent;
    for ( ; it != stop; ++it)
    {
        os << endl;
        it->first->print(os);
    }
    --Env::indent;
}

void
TileRegistry::register_tile(Tile** tileptr)
{
    Tile* tile = *tileptr;
    IndexRangeLocationPtr loc = new IndexRangeLocation(tile->get_index_ranges());
    tilemap_[loc] = tileptr;
}

Tile**
TileRegistry::get_tile(const IndexRangeTuplePtr &tuple)
{
    IndexRangeLocationPtr loc = new IndexRangeLocation(tuple);
    node_map::const_iterator it(tilemap_.find(loc));
    if (it == tilemap_.end())
    {
        cerr << "Invalid index range: " << endl;
        cerr << loc << endl;
        cerr << "Valid index ranges: " << endl;
        print(cerr); cerr << endl;
        raise(SanityCheckError, "invalid index range location passed to registry");
    }
    return it->second;
}

bool
TileRegistry::has_tile(const IndexRangeTuplePtr& tuple)
{
    IndexRangeLocationPtr loc = new IndexRangeLocation(tuple);
    node_map::iterator it(tilemap_.find(loc));
    return it != tilemap_.end();
}

void
TileRegistry::remove_tile(const IndexRangeTuplePtr &tuple)
{
    IndexRangeLocationPtr loc = new IndexRangeLocation(tuple);
    node_map::iterator it(tilemap_.find(loc));
    if (it != tilemap_.end())
        tilemap_.erase(it);
}

void
TileRegistry::delete_tile(Tile** tileptr)
{
    boost::intrusive_ptr_release(*tileptr);
    *tileptr = 0;
}

void
TileRegistry::delete_tiles()
{
    node_map::const_iterator it = tilemap_.begin();
    node_map::const_iterator stop = tilemap_.end();
    for ( ; it != stop; ++it)
    {
        delete_tile(it->second);
    }
    tilemap_.clear();
}

uli
TileRegistry::size() const
{
    return tilemap_.size();
}
