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
    else
        return false; //no symmetry
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
    tuple_(tuple),
    pgrp_(pgrp),
    nindex_(tuple->nindex()),
    tiles_(0),
    tmptuple_(reinterpret_cast<IndexRange**>(yeti_malloc_indexptr())),
    filler_(0),
    maxdepth_(tuple->maxdepth()),
    mindepth_(tuple->mindepth()),
    parents_(reinterpret_cast<Tile**>(yeti_malloc_indexptr())),
    index_strides_(yeti_malloc_indexset()),
    is_tiled_index_(yeti_malloc_perm()),
    rigorously_zero_tiles_(0),
    nparents_(1),
    estimater_(0)
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

    if (is_aligned())
        iterate_zero_check(0);
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
    nindex_(0),
    tiles_(0),
    tmptuple_(0),
    maxdepth_(tile->depth() + 1),
    mindepth_(tile->depth() + 1),
    filler_(0),
    rigorously_zero_tiles_(new char[1]),
    estimater_(0)
{
    parents_[0] = parent;
    uli size = 1;
    tiles_ = new CountableArray<Tile>(size);
    tiles_->insert(0, tile);
    rigorously_zero_tiles_[0] = 0;
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

void
TileMap::allow_all_tiles()
{
    ::memset(rigorously_zero_tiles_, 0, this->size() * sizeof(char));
}

TileMap::iterator
TileMap::begin() const
{
    return tiles_->begin();
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

void
TileMap::fill(const ThreadedTileElementComputerPtr& filler)
{
    filler_ = filler;
    if (filler_)
        estimater_ = filler_->get_estimater(maxdepth_);

    iterate_fill(0);
}

IndexRangeTuplePtr
TileMap::get_index_ranges() const
{
    return tuple_;
}

void
TileMap::insert(Tile* tile)
{
#if YETI_SANITY_CHECK
    if (!tile) //cannot insert null tile
    {
        raise(SanityCheckError, "cannot insert null tile into tile map");
    }
#endif

    uli idx = indexer_->index(tile->indices());
    insert(idx, tile);
}

void
TileMap::insert(uli index, Tile* tile)
{
#if YETI_SANITY_CHECK
    Tile* test = tiles_->get(index);
    if (test)
        raise(SanityCheckError, "cannot insert on occupied tile location");

    if (rigorously_zero_tiles_[index])
    {
        cerr << "cannot insert rigorously zero tile" << endl;
        abort();
    }
#endif
    tiles_->insert(index, tile);
    tile->set_parent(this); //claim ownership
}

bool
TileMap::is_aligned() const
{
    return maxdepth_ == mindepth_;
}

const usi*
TileMap::is_tiled_index() const
{
    return is_tiled_index_;
}

void
TileMap::iterate_fill(
    usi index
)
{
    if (index == nindex_)
    {
        insert_new_tile();
        return;
    }

    uli idx = offsets_[index];
    uli stop = nindices_[index] + idx;
    IndexRange* range = tuple_->get(index);
    //maybe or maybe not recurse on this
    bool use_subindex = range->depth() >= maxdepth_;
    for ( ; idx != stop; ++idx) //loop all index values within the given range
    {
        //set the index number at the given index position
        tmpindices_[index] = idx;

        //set the IndexRange at the given index position
        tmptuple_[index] = use_subindex ? range->get_subindex(idx) : range;
        iterate_fill(index + 1);
    }
}

void
TileMap::insert_new_tile()
{
    uli index = indexer_->index(tmpindices_);
    if (rigorously_zero_tiles_[index])
        return;

    if (estimater_ && estimater_->max_log(tmpindices_) < YetiRuntime::matrix_multiply_cutoff)
        return;

    Tile* tile = get(index);
    if (tile) //already exists
        return;

    //this corresponds to a permutationally unique block
    IndexRangeTuplePtr subtuple
        = new IndexRangeTuple(nindex_, tmptuple_);

    //create the tile
    uli* indexarr = yeti_malloc_indexset();
    ::memcpy(indexarr, tmpindices_, nindex_ * sizeof(size_t));

    //make(idxset, IndexSet, indexarr, nindex_);
    tile = new Tile(subtuple, indexarr, this);
    insert(index, tile);
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

void
TileMap::tile_zero_check()
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
        }
    }
}

void
TileMap::iterate_zero_check(usi index)
{
    if (index == nindex_)
    {
        tile_zero_check();
        return;
    }

    uli idx = offsets_[index];
    uli stop = nindices_[index] + idx;
    for ( ; idx != stop; ++idx) //loop all index values within the given range
    {
        //set the index number at the given index position
        tmpindices_[index] = idx;
        iterate_zero_check(index + 1);
    }
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

    os << Env::indent << "Tile Map at " << (void*) this << " nref=" << nref();
    for (usi i=0; i < nparents_; ++i)
    {
        os << endl << Env::indent << "Parent " << i << ": " << (void*) parents_[i];
    }
    os << endl << Env::indent << "Sizes: " << ClassOutput<const uli*>::str(nindex_, nindices());
    os << endl << Env::indent << "Offsets: " << ClassOutput<const uli*>::str(nindex_, offsets());
    //os << endl << indexer_;

    usi idx = 0;
    foreach_nonnull(tile, this, Tile,
        os << endl << Env::indent;
        ++Env::indent;
        tile->retrieve(NOT_THREADED);
        os << idx << " -> ";
        tile->print(os, false); //do not indent header
        tile->release(NOT_THREADED);
        --Env::indent;
        ++idx;
    )
}

void
TileMap::retrieve()
{
}

void
TileMap::release()
{
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

    if (idx == nparents_)
        return;

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
        incref();
        tilemap_ = config_->tile_map_builder->build_map(
                ranges_,
                this
                );
        decref();
    }
    init_estimate();
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices,
    //const PermutationGroupPtr& pgrp,
    const TensorConfigurationPtr& config
) : ranges_(ranges),
    //pgrp_(pgrp),
    indices_(indices),
    parent_(0),
    tilemap_(0),
    owner_process_(0),
    has_owner_process_(false), //top level tile, so no owner,
    data_(0),
    config_(config)
{
    init();
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices,
    //const PermutationGroupPtr& pgrp,
    TileMap* parent
) : ranges_(ranges),
    //pgrp_(pgrp),
    indices_(indices),
    parent_(parent),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(parent->config())
{
    init();
}

Tile::Tile(
    const IndexRangeTuplePtr& ranges,
    size_t* indices
) : ranges_(ranges),
    indices_(indices),
    parent_(0),
    owner_process_(0),
    has_owner_process_(false), //no owner... yet
    data_(0),
    config_(0)
{
    init();
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
    indices_(0)
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

    init();

    estimate_.depth = 0;
    estimate_.maxlog = LOG_ZERO;
}

Tile::~Tile()
{
    if (tilemap_)
        tilemap_->remove_parent(this);
    yeti_free_indexset(indices_);
}

void
Tile::add_tiles(const TileIteratorPtr& iter)
{
    if (depth() == 1)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            iter->add_tile(tile);
        )
    }
    else if (tilemap_)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            tile->add_tiles(iter);
        )
    }
    else
    {
        raise(SanityCheckError, "Cannot add tiles in data tile");
    }
}

void
Tile::accumulate(
    const TilePtr& tile,
    const ThreadedSortPtr& sort,
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
            Sort* sorter = sort->get_sorter(0);
            TileIteratorWorkspace* workspace
                = get_workspace<TileIteratorWorkspace>(NOT_THREADED);
            usi depth = tilemap_->maxdepth();
            Tile** sorted_tiles = reinterpret_cast<Tile**>(workspace->buffers[depth]);
            sorter->configure(tile->get_map()->nindices());
            sorter->sort_noscale<Tile*>(unsorted_tiles, sorted_tiles);
            acc = sorted_tiles;
        }

        uli idx = 0;
        for ( ; it != stop; ++it, ++acc, ++idx)
        {
            Tile* tile_acc = *acc;
            if (!tile_acc) //nothing to accumulate
                continue;

            tile_acc->retrieve(NOT_THREADED);

            cout << "accumulating index " << idx << " at depth " << depth() << endl;

            Tile* tile_me = *it;
            if (!tile_me)
            {
                if (tilemap_->is_rigorously_zero(it))
                {
                    cout << "tile rigorously zero" << endl;
                    continue; //nothing to do here
                }
                else
                {
                    uli* newindices = yeti_malloc_indexset();
                    ::memcpy(newindices, tile_acc->indices(), nindex());
                    tile_me = new Tile(
                            tile_acc->get_index_ranges(),
                            newindices,
                            config_
                         );
                    tilemap_->insert(idx, tile_me);
                    cout << "accumulating new tile" << endl;
                }
            }


            tile_me->retrieve(NOT_THREADED);
            tile_me->accumulate(tile_acc, sort, scale);
            tile_me->release(NOT_THREADED);

            tile_acc->release(NOT_THREADED);
        }
    }
    else //data tile
    {

        if (data_->data()->null()) //register my allocation
            config_->data_factory->register_allocation(data_.get());

        AccumulateTask* task = new AccumulateTask(
            data_.get(),
            tile->get_data(),
            tile->get_index_ranges(),
            sort,
            scale
        );

        GlobalQueue::add(this, task);
    }
}

void
Tile::accumulate(
    const MatrixPtr& matrix,
    const ThreadedSortPtr& sort,
    double scale
)
{
    if (tilemap_)
    {
        uli* tmp = yeti_malloc_indexset();

        Tile** it(tilemap_->begin());
        Tile** stop(tilemap_->end());

        Matrix** unsorted_blocks = matrix->get_map()->begin();
        TileMap* srcmap = matrix->get_unique_tile()->get_map();
        Permutation* src_perm = matrix->get_unique_permutation();

        Matrix** acc;
        if (sort->get_permutation()->is_identity())
        {
            acc = unsorted_blocks;
        }
        else
        {
            Sort* sorter = sort->get_sorter(0);
            TileIteratorWorkspace* workspace
                = get_workspace<TileIteratorWorkspace>(NOT_THREADED);
            usi depth = tilemap_->maxdepth();
            Matrix** sorted_blocks = reinterpret_cast<Matrix**>(workspace->buffers[depth]);


            src_perm->permute(srcmap->nindices(), tmp);
            sorter->configure(tmp);
            sorter->sort_noscale<Matrix*>(unsorted_blocks, sorted_blocks);
            acc = sorted_blocks;
        }

        uli idx = 0;
        for ( ; it != stop; ++it, ++acc, ++idx)
        {
            Matrix* src = *acc;
            if (!src) //nothing to accumulate
                continue;

            Tile* src_tile = src->get_unique_tile();
            src->retrieve(NOT_THREADED);

            Tile* dst = *it;
            if (!dst)
            {
                if (tilemap_->is_rigorously_zero(it))
                {
                    continue; //nothing to do here
                }
                else
                {
                    uli* newindices = yeti_malloc_indexset();
                    src->get_unique_permutation()
                            ->permute(src_tile->indices(), newindices);
                    dst = new Tile(
                            src->get_index_ranges(),
                            newindices,
                            config_
                         );
                    tilemap_->insert(idx, dst);
                }
            }


            dst->retrieve(NOT_THREADED);
            dst->accumulate(src, sort, scale);
            dst->release(NOT_THREADED);
            src->release(NOT_THREADED);
        }

        yeti_free_indexset(tmp);
    }
    else //data tile
    {

        if (data_->data()->null()) //register my allocation
            config_->data_factory->register_allocation(data_.get());

        AccumulateTask* task = new AccumulateTask(
            data_.get(),
            matrix->get_data(),
            matrix->get_index_ranges(),
            sort,
            scale
        );

        GlobalQueue::add(this, task);
    }
}

AccumulateTask::AccumulateTask(
    DataBlock *dst,
    DataBlock *src,
    IndexRangeTuple* src_ranges,
    const ThreadedSortPtr& sort,
    double scale
)
    :
  dst_(dst),
  src_(src),
  src_ranges_(src_ranges),
  sort_(sort),
  scale_(scale)
{
}

void
AccumulateTask::print(std::ostream &os) const
{
    os << "Accumulate Task" << endl;
    os << "Destination: " << dst_ << endl;
    os << "Source: " << src_;
}

void
AccumulateTask::run(uli threadnum)
{
    src_->retrieve(threadnum);
    dst_->retrieve(threadnum);

    Data* src = src_->data();
    Data* dst = dst_->data();

    Sort* sorter = sort_->get_sorter(threadnum);
    if (sorter)
        sorter->configure(src_ranges_);
    dst->accumulate(src, scale_, sorter);

    src_->release(threadnum);
    dst_->release(threadnum);
}

TensorConfiguration*
Tile::config() const
{
    return config_.get();
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
    if (DATA_DISTRIBUTION_DEPTH == tilemap_->maxdepth()) //distribute here
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
            bool eq = false;
            if (data_->data()->nonnull())
            {
                eq = data_->data()->equals(*data);
                (*data) += data_->data()->size();
            }
            data_->release(NOT_THREADED);
            return eq;
        }
        else //must have data!
            return false;
    }
}

bool
Tile::equals(const TileElementComputerPtr& filler)
{
    ThreadedTileElementComputerPtr thread_filler = new ThreadedTileElementComputer(filler);
    return this->equals(thread_filler);
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

bool
Tile::exists(const uli *indices) const
{
    return tilemap_->exists(indices);
}

void
Tile::fill()
{
    if (config_->filler->mindepth() > depth())
        return; //we are done here

    if (tilemap_)
    {
        tilemap_->retrieve();
        tilemap_->fill(config_->filler);
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
    else if (storage_type() != Data::recomputed) //data tile
    {

        FillTask* task = new FillTask(this, config_->filler.get());
        GlobalQueue::add(this, task);

        if (data_->data()->null()) //register my allocation
            config_->data_factory->register_allocation(data_.get());

    }
}

FillTask::FillTask(Tile *tile, ThreadedTileElementComputer *filler)
    :
  tile_(tile),
  filler_(filler)
{
}

void
FillTask::run(uli threadnum)
{
    DataBlock* data = tile_->get_data();
    data->retrieve(threadnum);
    filler_->compute(tile_, data->data(), threadnum);
    float maxlog = data->data()->max_log();
    data->release(threadnum);

    //setting log from min depth
    tile_->set_max_log(maxlog, 0);
}

void
FillTask::print(std::ostream &os) const
{
    os << "Fill Task" << endl;
    os << tile_;
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

TileMap*
Tile::get_parent() const

{
    return parent_;
}

TileIteratorPtr
Tile::get_iterator()
{
    uli ntiles = ntiles_nonnull();
    TileIteratorPtr iter = new TileIterator(ntiles);
    add_tiles(iter);
    return iter;
}

IndexRangeTuple*
Tile::get_index_ranges() const
{
    return ranges_.get();
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

uli*
Tile::get_location(const usi* subset, usi nsub) const
{
    uli* location = yeti_malloc_tile_location();
    get_location(location, subset, nsub);
    return location;
}

void
Tile::get_location(uli* location, const usi* subset, usi nsub) const
{
    if (parent_)
    {
        Tile* parent_tile = parent_->get_parent_tile();
        usi mydepth = depth();
        if (subset)
        {
            const uli* sizes = parent_->nindices();
            const uli* offsets = parent_->offsets();
            const usi* subsetptr = subset + nsub - 1;
            uli stride = 1;
            uli index = 0;
            for ( ; subsetptr >= subset; --subsetptr)
            {
                usi subidx = *subsetptr;
                index += (indices_[subidx] - offsets[subidx]) * stride;
                stride *= sizes[subidx];
            }
            location[mydepth] = index;
        }
        else
        {
            location[mydepth] = parent_->index(indices_);
        }
        parent_tile->get_location(location, subset, nsub);
    }
    else
    {
        return; //top level... do nothing
    }
}

TileMap*
Tile::get_tile_map(const uli* location, bool create_location)
{
    if (!tilemap_)
    {
        raise(SanityCheckError, "cannot call get tile on tile without tile map");
    }

    usi seekdepth = depth() - 1;
    if (seekdepth == 0)
    {
        return tilemap_.get();
    }

    uli compidx = location[seekdepth];
    Tile* parent = tilemap_->get(compidx);
    if (!parent && create_location) //create the parent tile
    {
        usi nidx = nindex();
        uli* indices = yeti_malloc_indexset();
        tilemap_->compute_indices(compidx, indices);
        IndexRangeTuplePtr tuple = new IndexRangeTuple(nidx);
        for (usi i=0; i < nidx; ++i)
        {
            IndexRange* range =  ranges_->get(i)->get_subindex(indices[i]);
            tuple->set(i, range);
        }
        parent = new Tile(tuple, indices, tilemap_.get());
        tilemap_->insert(compidx, parent);
    }

    if (parent)
    {
        parent->retrieve(NOT_THREADED);
        TileMap* tilemap = parent->get_tile_map(location, create_location);
        parent->release(NOT_THREADED);
        return tilemap;
    }
    else
    {
        return 0;
    }

}

Tile*
Tile::get_tile(const uli* location)
{
    TileMap* tilemap = get_tile_map(location, false); //do not create the location
    if (tilemap)
        return tilemap->get(location[0]);
    else
        return 0;
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

const size_t*
Tile::indices() const
{
    return indices_;
}

void
Tile::insert(const uli* location, Tile* tile)
{
    TileMap* tilemap = get_tile_map(location, true); //create the location
    tilemap->insert(location[0], tile);
}

bool
Tile::is_aligned() const
{
    return tilemap_->depths_aligned();
}

bool
Tile::is_equivalent(const TilePtr& tile) const
{
    if (type() != tile->type())
        return false;

    if (get_index_ranges()->nindex() != tile->get_index_ranges()->nindex())
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
    if (ranges_->nindex() == 0)
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
    ranges_->nindex();
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
    uli nmax = 1;
    IndexRangeTuple::iterator it(ranges_->begin());
    IndexRangeTuple::iterator stop(ranges_->end());
    for ( ; it != stop; ++it)
    {
        IndexRange* range = *it;
        uli n = range ? range->n() : 1;
        nmax *= n;
    }
    return nmax;
}

uli
Tile::ntiles_nonnull()
{
    if (depth() == 1)
    {
        uli ntiles = 0;
        foreach_nonnull(tile, tilemap_, Tile,
            ++ntiles;
        )
        return ntiles;
    }
    else if (tilemap_)
    {
        uli ntiles = 0;
        foreach_nonnull(tile, tilemap_, Tile,
            ntiles += tile->ntiles_nonnull();
        )
        return ntiles;
    }
    else
    {
        raise(SanityCheckError, "cannot call ntiles on data tile");
    }
}

uli
Tile::ntiles_owned()
{
    return ntiles_owned(DATA_DISTRIBUTION_DEPTH);
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
Tile::parameter_reduce(usi nparams)
{
    if (tilemap_)
    {
        raise(SanityCheckError, "cannot call parameter reduce on a tile with a tile map");
    }


    usi offset = nparams;
    usi nidx = nindex() - nparams;
    for (usi i=0; i < nidx; ++i)
    {
        indices_[i] = indices_[i + offset];
    }

    ranges_->slice_front(nparams);
}

void
Tile::print(std::ostream &os) const
{
    print(os, true);
}

void
Tile::print(ostream& os, bool indent_header) const
{
    if (indent_header)
        os << Env::indent;

    os << ClassOutput<const size_t*>::str(nindex(), indices_) << "  " << ranges_;
    if (has_owner_process_)
        os << " Node " << owner_process_;

    os << "  Max: 10^(" << stream_printf("%8.4f", estimate_.maxlog)
        << ") from depth " << estimate_.depth;

    os << " nref=" << nref();

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
    if (tilemap_) //metadata tile
    {
        tilemap_->retrieve();
        if (config_->data_factory->storage_type() == Data::recomputed)
        {
            ThreadedTileElementComputer* filler =
                config_->data_factory->get_element_computer();
            usi predepth = filler->mindepth();
            filler->set_mindepth(tilemap_->maxdepth());
            tilemap_->fill(filler);
            filler->set_mindepth(predepth);
        }
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
Tile::set_parent(TileMap* tilemap)
{
    parent_ = tilemap;
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
        data_->retrieve(NOT_THREADED);
        data_->sort(sort, buffer);
        data_->release(NOT_THREADED);
    }

    //now sort the index ranges
    IndexRange** rangeptr = ranges_->begin();
    IndexRange** rangebuf = reinterpret_cast<IndexRange**>(buffer);
    sort->get_permutation()->permute(rangeptr, rangebuf);
    ::memcpy(rangeptr, rangebuf, ranges_->nindex() * sizeof(IndexRange*));

    //and the indices
    size_t* indexbuf = reinterpret_cast<size_t*>(buffer);
    sort->get_permutation()->permute(indices_, indexbuf);
    ::memcpy(indices_, indexbuf, nindex() * sizeof(size_t));
}

void
Tile::sort(const PermutationPtr& p)
{
    SortPtr s(new Sort(p, ranges_));

    uli max_blocksize = this->max_blocksize();
    uli index_size = YetiRuntime::max_nindex() * sizeof(void*);
    uli buffer_size = max_blocksize > index_size ? max_blocksize : index_size;

    void* buffer = malloc(buffer_size);

    _sort(s, buffer);
    free(buffer);
}

Data::storage_t
Tile::storage_type() const
{
    return config_->data_factory->storage_type();
}

void
Tile::tally(const TileDistributerPtr &distr)
{
    if (!tilemap_)
        raise(SanityCheckError, "cannot distribute on tile without data map");

    tilemap_->retrieve();
    if (DATA_DISTRIBUTION_DEPTH == tilemap_->maxdepth()) //distribute here
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

TileIterator::TileIterator(uli n) :
    tiles_(n),
    count_(0)
{
}

Tile**
TileIterator::begin() const
{
    return tiles_.begin();
}

Tile**
TileIterator::end() const
{
    return tiles_.end();
}

void
TileIterator::add_tile(Tile *tile)
{
    tiles_.insert(count_, tile);
    ++count_;
}

uli
TileIterator::ntiles() const
{
    return tiles_.n();
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

