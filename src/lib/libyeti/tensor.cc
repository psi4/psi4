#include "tensor.h"
#include "contraction.h"
#include "index.h"
#include "class.h"
#include "exception.h"
#include "permutation.h"
#include "env.h"
#include "data.h"
#include "matrix.h"
#include "runtime.h"
#include "dataimpl.h"

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

TensorConfiguration::TensorConfiguration(
    const PermutationGroupPtr& grp,
    const std::string& tensor_name
)
    :
    tile_distribution_type(Tile::replicated),
    map_storage_type(TileMap::permanent),
    distr_indices(yeti_malloc_perm()),
    nindex_distr(grp->nindex()),
    print_data(true),
    data_factory(0),
    priority(Tensor::alpha_tensor),
    cache(new LayeredDataCache),
    tile_map_builder(0),
    filler(0),
    name(tensor_name),
    data_mode(new DataMode)
{
    filters.push_back(new PermutationalSymmetryFilter(grp));
}

TensorConfiguration::TensorConfiguration(
)
    :
    tile_distribution_type(Tile::replicated),
    map_storage_type(TileMap::permanent),
    distr_indices(yeti_malloc_perm()),
    nindex_distr(0),
    print_data(true),
    data_factory(0),
    priority(Tensor::alpha_tensor),
    cache(new LayeredDataCache),
    tile_map_builder(0),
    filler(0),
    name(""),
    data_mode(new DataMode)
{
}

TensorConfiguration::~TensorConfiguration()
{
    delete data_mode;
    yeti_free_perm(distr_indices);
}

TensorConfigurationPtr
TensorConfiguration::copy(const std::string& name) const
{
    TensorConfiguration* copy = new TensorConfiguration();

    std::list<TileFilterPtr>::const_iterator it(filters.begin());
    std::list<TileFilterPtr>::const_iterator stop(filters.end());
    for ( ; it != stop; ++it)
    {
        copy->filters.push_back((*it)->copy());
    }

    copy->name = name;
    copy->tile_distribution_type = tile_distribution_type;
    copy->nindex_distr = nindex_distr;
    copy->print_data = print_data;
    copy->priority = priority;
    copy->filler = filler;
    copy->tile_map_builder = tile_map_builder;
    copy->data_mode->flag = data_mode->flag;

    //This should be malloc'd by the tensor configuration constructo
    //copy->data_mode = new DataMode;
    //copy->distr_indices = yeti_malloc_perm();
    ::memcpy(copy->distr_indices, distr_indices, nindex_distr * sizeof(usi));

    /** CHANGE! For now just use the same cache */
    copy->cache = cache;

    if (data_factory)
        copy->data_factory = data_factory->copy();

    return copy;
}

const DataMode*
TensorConfiguration::get_data_mode() const
{
    return data_mode;
}

Tensor::Tensor(
    const std::string& name,
    const IndexRangeTuplePtr& tuple,
    const PermutationGroupPtr& pgrp
) :
    pgrp_(pgrp),
    Tile(
        tuple,
        IndexSet::get_zero_set(),
        new TensorConfiguration(pgrp, name)
    ),
    alignment_depth_(0)
{
    config_->tile_map_builder
        = new DefaultTileMapBuilder(pgrp_);

    init();

    for (usi i=0; i < nindex(); ++i)
        config_->distr_indices[i] = i;
}

Tensor::Tensor(
    const std::string& name,
    const IndexRangeTuplePtr& tuple,
    const PermutationGroupPtr& pgrp,
    const TensorConfigurationPtr& config
) :
    pgrp_(pgrp),
    Tile(
        tuple,
        IndexSet::get_zero_set(),
        config
    ),
    alignment_depth_(0)
{
    init();
}

Tensor::Tensor(const TensorPtr& t)
    :
    pgrp_(t->pgrp_),
    Tile(
        IndexRangeTuple::get_unit_range(t->get_index_ranges()),
        IndexSet::get_zero_set()
      )
{
    config_ = t->config_;
    ::memcpy(
        config_->distr_indices,
        t->config_->distr_indices,
        sizeof(usi) * t->config_->nindex_distr
    );

    init(); //no depth restrictions

    if (t->nindex() == 0) //create a null tilemap
    {
        incref();
        tilemap_ = new TileMap(t.get(), this);
        t->parent_ = tilemap_.get();
        decref();
    }
    else
    {
        tilemap_->insert(t.get());
    }
}

Tensor::Tensor(const TensorConfigurationPtr& config)
    : Tile(config)
{
}

Tensor::~Tensor()
{
}

void*
Tensor::operator new(uli n)
{
    return ::malloc(sizeof(Tensor));
}

void
Tensor::operator delete(void* ptr)
{
    ::free(ptr);
}

void
Tensor::accumulate(
    const TensorPtr &tile,
    const PermutationPtr &p,
    double scale
)
{
#if YETI_SANITY_CHECK
//    if (!tile->get_permutation_grp()->equals(this->pgrp_))
//        raise(SanityCheckError, "cannot accumulate tensors with incompatible symmetries");
#endif
    SortPtr sort = new Sort(p, ranges_);
    Tile::accumulate(tile, sort, scale);
}

void
Tensor::allocate()
{
    if (is_allocated())
        return;

    //if we have a registry of tiles
    if (registry_) //make sure we only build non-null tiles
    {

        std::list<IndexRangeTuplePtr>::const_iterator it(nonnull_tuples_.begin());
        std::list<IndexRangeTuplePtr>::const_iterator stop(nonnull_tuples_.end());
        for ( ; it != stop; ++it)
        {
            if (!registry_->has_tile(*it))
            {

            }
        }

        it = nonnull_tuples_.begin();
        for ( ; it != stop; ++it)
        {
            registry_->remove_tile(*it);
        }

        //clear all remaining tiles in the registry
        registry_->delete_tiles();
    }

    pgrp_->close();

    retrieve(NOT_THREADED);

    if (config_->map_storage_type == TileMap::temporary)
        return; //nothing to allocate

    distribute();

    Tile::fill(); //construct the metadata
    config_->data_factory->configure(this);
    Tile::allocate();

    release(NOT_THREADED);
}

Tensor*
Tensor::copy(const std::string& name)
{
    retrieve(NOT_THREADED);
    TensorConfigurationPtr newconfig = config_->copy(name);
    Tensor* copy = new Tensor(newconfig);
    copy->alignment_depth_ = alignment_depth_;
    copy->pgrp_ = pgrp_;

    Tile::copy_to(copy);

    release(NOT_THREADED);
    return copy;
}

void
Tensor::distribute()
{
    pgrp_->close();

    //distribution will require us to have a data factory
    if (!config_->data_factory) //this must exist in memory
        config_->data_factory = new MemoryBlockFactory<double>;

    if (YetiRuntime::nproc() == 1)
    {
        return;
    }

    return;

    TileDistributerPtr distr(
        new TileDistributer(config_->distr_indices, config_->nindex_distr)
    );

    //now that I have the distribution object, tally all of the tiles
    //at a given level
    foreach_nonnull(tile, tilemap_, Tile,
        tile->tally(distr);
    )

    distr->sort_tasks();

    //now that all of the tiles have been tallied, sort the distribution
    foreach_nonnull(tile, tilemap_, Tile,
        tile->distribute(distr);
    )

}

void
Tensor::init()
{
    parent_ = 0;

    if (ranges_->size() > 0)
        alignment_depth_ = ranges_->mindepth();
}

usi
Tensor::alignment_depth() const
{
    return alignment_depth_;
}

PermutationGroup*
Tensor::antisymmetric_permutation_group(
    const IndexRangeTuplePtr &bratuple,
    const IndexRangeTuplePtr &kettuple
)
{
    usi nbra = bratuple->size();
    usi nket = bratuple->size();

    if (nbra != nket)
        raise(SanityCheckError, "number of bra and ket indices must be the same");

    usi nidx = nbra + nket;

    short minus = -1;
    short plus = 1;

    PermutationGroup*  grp = new PermutationGroup(nidx);

    //determine all transpositions in the bra
    usi offset = 0; //no offset for bra indices
    for (usi i=0; i < nbra; ++i)
    {
        IndexRange* idx_i(bratuple->get(i));
        for (usi j=i+1; j < nbra; ++j)
        {
            IndexRange* idx_j(bratuple->get(j));
            if (idx_i->equals(idx_j))
            {
                //add the transposition
                make(p, Permutation, nidx, i + offset, j + offset, minus);
                grp->add(p);
            }
        }
    }

    offset = nbra; //start ket indices after bra indices
    for (usi i=0; i < nket; ++i)
    {
        IndexRange* idx_i(kettuple->get(i));
        for (usi j=i+1; j < nket; ++j)
        {
            IndexRange* idx_j(kettuple->get(j));
            if (idx_i->equals(idx_j))
            {
                //add the transposition
                make(p, Permutation, nidx, i + offset, j + offset, minus);
                grp->add(p);
            }
        }
    }

    grp->close();

    return grp;
}

Tensor*
Tensor::build(
    const std::string& name,
    usi nidx,
    IndexRange* range
)
{
    PermutationGroupPtr grp = new PermutationGroup(nidx);
    IndexRangeTuplePtr tuple = new IndexRangeTuple(nidx, range);

    return new Tensor(name, tuple, grp);
}

template <class data_t>
data_t
Tensor::dot_product(const TensorPtr &tensor)
{
    if (tensor->nindex() != nindex())
        raise(SanityCheckError, "different number of indices in tensor dot product");

    //create two matrices
    usi nrows = 0;
    usi ncols = 0;
    usi nlink = nindex();

    MatrixIndexPtr lindex = new MatrixIndex(MatrixIndex::front, nlink, MatrixIndex::back, nrows);
    MatrixIndexPtr rindex = new MatrixIndex(MatrixIndex::front, nlink, MatrixIndex::back, ncols);

    PermutationGroupPtr scalar_grp = new PermutationGroup(0);
    IndexRangeTuplePtr tuple = new IndexRangeTuple(0);

    TensorPtr product_tensor = new Tensor("dot product", tuple, scalar_grp);
    product_tensor->config_->data_factory = new MemoryBlockFactory<data_t>;

    ContractionPtr cxn = new Contraction(
        1.0,
        this, tensor,
        lindex, rindex,
        scalar_grp,
        pgrp_->get_identity(),
        tensor->get_permutation_grp()->get_identity(),
        product_tensor
    );

    TensorPtr product = cxn->get_product_tensor();

    Contraction::configure();
    Contraction::run();



    //get the data block
    uli index = 0;
    Data* data = product->get_map()->get(index)->get_data()->data();
    data_t* d = data->get<data_t>();
    return *d;
}

int
Tensor::dot_product_int(const TensorPtr& tensor)
{
    return dot_product<int>(tensor);
}

quad
Tensor::dot_product_quad(const TensorPtr& tensor)
{
    return dot_product<quad>(tensor);
}

double
Tensor::dot_product_double(const TensorPtr& tensor)
{
    return dot_product<double>(tensor);
}

float
Tensor::dot_product_float(const TensorPtr& tensor)
{
    return dot_product<float>(tensor);
}

bool
Tensor::equals(const void *data)
{
    return Tile::equals(data);
}


bool
Tensor::equals(const ThreadedTileElementComputerPtr &filler)
{
    return Tile::equals(filler);
}

bool
Tensor::equals(const TensorPtr& tensor)
{
    Env::err0() << "not yet implemented at"
        << __FILE__ << " " << __LINE__ << endl;
    abort();
}

void
Tensor::fill()
{
    fill(config_->filler);
}

void
Tensor::fill(const ThreadedTileElementComputerPtr& new_filler)
{
    config_->filler = new_filler;

    allocate();

    set_write_mode();
    Tile::fill();
    set_read_mode();
}

void
Tensor::fill(const TileElementComputerPtr& val)
{
    if (val)
    {
        fill(new ThreadedTileElementComputer(val));
    }
    else
    {
        ThreadedTileElementComputer* nullptr = 0;
        fill(nullptr);
    }

    //ThreadedTileElementComputerPtr filler = new TileFiller(val);
    //fill(filler);
}

PermutationGroup*
Tensor::get_permutation_grp() const
{
    return pgrp_.get();
}

const std::string&
Tensor::name() const
{
    return config_->name;
}

void
Tensor::print(ostream& os) const
{
    os << Env::indent << "Tensor: " << config_->name << endl;
    ++Env::indent;
    os << Env::indent << "Total symmetry: " << endl
       << pgrp_ << endl;
    --Env::indent;

    Tile::print(os);
}

Tensor::contraction_priority_t
Tensor::priority() const
{
    return config_->priority;
}

void
Tensor::push_location(const TileLocationPtr& loc) const
{
    //do nothing... this is the end
}

void
Tensor::register_nonnull_tiles(const IndexRangeTuplePtr &tuple)
{
    usi depth = tuple->mindepth();
    if (!registry_)
        registry_ = new TileRegistry;

    retrieve(NOT_THREADED);
    //register tiles at the given depth
    register_tiles(registry_, depth);

    //register all the tuples at the depth
    if (tuple->is_aligned())
    {
        nonnull_tuples_.push_back(tuple);
    }
    else
    {
        std::list<IndexRangeTuplePtr> tuples;
        tuple->get_index_tuples(nonnull_tuples_, depth);
    }

    release(NOT_THREADED);
}

void
Tensor::finalize_mode()
{
    switch(config_->data_mode->flag)
    {
        case DataMode::read: break; //do nothing
        case DataMode::write: config_->cache->finalize_write(); break;
        case DataMode::accumulate: config_->cache->finalize_accumulate(); break;
    }

}

void
Tensor::set_read_mode()
{
    finalize_mode();
    config_->data_mode->flag = DataMode::read;
}

void
Tensor::set_write_mode()
{
    finalize_mode();
    config_->data_mode->flag = DataMode::write;
}

void
Tensor::set_accumulate_mode()
{
    finalize_mode();
    config_->data_mode->flag = DataMode::accumulate;
}

void
Tensor::sizes(std::map<size_t, size_t>& final_size_counts) const
{
    std::map<size_t, size_t> total_size_counts;
    std::map<size_t, size_t> range_size_counts;

    //IndexRangeTuple* tuple = tensor->get_index_ranges();

    final_size_counts[sizeof(double)] = 1; //start the map off

    foreach(range, ranges_, IndexRange,

        range->sizes(range_size_counts);


        total_size_counts = final_size_counts;
        final_size_counts.clear();
        std::map<size_t, size_t>::iterator itmap = total_size_counts.begin();
        std::map<size_t, size_t>::iterator stopmap = total_size_counts.end();
        for ( ; itmap != stopmap; ++itmap)
        {
            size_t total_size = itmap->first;
            size_t total_count = itmap->second;
            std::map<size_t, size_t>::iterator itlist(range_size_counts.begin());
            std::map<size_t, size_t>::iterator stoplist(range_size_counts.end());
            for ( ; itlist != stoplist; ++itlist)
            {
                size_t range_size = itlist->first;
                size_t range_count = itlist->second;
                final_size_counts[total_size*range_size] += total_count*range_count;
            }
        }

        total_size_counts.clear();
        range_size_counts.clear();
    )
}

void
Tensor::sort(const PermutationPtr &p)
{
    retrieve(NOT_THREADED);
    Tile::sort(p);
    release(NOT_THREADED);
    //get the permuted group
    PermutationGroupPtr newgrp = pgrp_->conjugate(p);
    pgrp_->clear();
    pgrp_->add(newgrp);
    pgrp_->close();
}

bool
yeti::tensor_less(const TensorPtr &l, const TensorPtr &r)
{
    //distributed tensors are always greater than replicted tensors
    if (l->config()->tile_distribution_type != r->config()->tile_distribution_type)
        return l->config()->tile_distribution_type < r->config()->tile_distribution_type;

    //in core is always faster than recomputing which is probably faster than going to disk
    Data::storage_t lstorage = l->config()->data_factory->storage_type();
    Data::storage_t rstorage = r->config()->data_factory->storage_type();
    if (lstorage != rstorage)
        return lstorage < rstorage;

    //everything is the same now, select based on size
    size_t lsize = l->total_size();
    size_t rsize = r->total_size();
    if (lsize < rsize)
        return lsize < rsize;

    return false; //exactly equal
}

size_t
Tensor::total_size() const
{
    size_t ntot = 1;
    IndexRangeTuple::iterator it(ranges_->begin());
    IndexRangeTuple::iterator stop(ranges_->end());
    for ( ; it != stop; ++it)
    {
        IndexRange* range(*it);
        ntot *= range->ntot();
    }
    return ntot;
}


SubtensorTileMapBuilder::SubtensorTileMapBuilder(
    Tensor* tensor,
    const IndexRangeTuplePtr& subtuple
)
    :
    tensor_(tensor),
    subtuple_(subtuple),
    tile_registry_(new TileRegistry),
    subdepth_(MAX_DEPTH)
{
    //figure out the recursion depth for the tile map registry
    IndexRangeTuplePtr tuple = tensor->get_index_ranges();
    for (usi i=0; i < subtuple_->size(); ++i)
    {
        IndexRange* parent = tuple->get(i);
        IndexRange* child = subtuple->get(i);
        usi depth = parent->get_subdepth_alignment(child);
        if (depth < subdepth_)
            subdepth_ = depth;
    }

    /**
        Now that we have the alignment depth check the tensor to see what type of "subtensor" we have.
        If the subtensor alignment depth is equal to the parent tensor alignment depth,
        then we have to log all the tile maps from the alignment depth.
    */
    tensor_->register_tiles(tile_registry_, subdepth_);
}

SubtensorTileMapBuilder::~SubtensorTileMapBuilder()
{
}

TileMapPtr
SubtensorTileMapBuilder::build_map(
    const IndexRangeTuplePtr &tuple,
    const TilePtr &parent
)
{
    bool reuse = true;
    for (usi i=0; i < tuple->size(); ++i)
    {
        if (tuple->get(i)->depth() != subdepth_)
        {
            reuse = false;
            break;
        }
    }

    if (reuse)
    {
        Tile* tile = *(tile_registry_->get_tile(tuple));
        if (!tile)
        {
            cerr << "Subtensor built from "
                << tensor_->name() << " is not a proper subtensor"
                << endl;
            abort();
        }
        parent->set_max_log(tile->get_max_log());
        tile->get_map()->add_parent(parent.get());
        return tile->get_map();
    }
    else
    {
        return new TileMap(tuple, tensor_->get_permutation_grp(), parent);
    }

}



