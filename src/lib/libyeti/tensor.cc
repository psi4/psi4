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
    data_factory(new MemoryBlockFactory<double>),
    priority(Tensor::alpha_tensor),
    cache(new LayeredDataCache),
    tile_map_builder(new DefaultTileMapBuilder(grp)),
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

DataBlockFactory*
TensorConfiguration::get_data_factory() const
{
    return data_factory.get();
}

Tensor::Tensor(
    const std::string& name,
    const IndexRangeTuplePtr& tuple,
    const PermutationGroupPtr& pgrp
) :
    pgrp_(pgrp),
    Tile(
        tuple,
        IndexRange::get_zero_set(),
        new TensorConfiguration(pgrp, name)
    ),
    nparams_(0),
    params_(0),
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
        IndexRange::get_zero_set(),
        config
    ),
    nparams_(0),
    params_(0),
    alignment_depth_(0)
{
    init();
}

Tensor::Tensor(
    const std::string& name,
    const IndexRangeTuplePtr& tuple,
    const PermutationGroupPtr& pgrp,
    const TensorConfigurationPtr& config,
    uli* params,
    usi nparams
) :
    pgrp_(pgrp),
    Tile(
        tuple,
        IndexRange::get_zero_set(),
        config
    ),
    nparams_(nparams),
    params_(params),
    alignment_depth_(0)
{
    init();
}

Tensor::Tensor(const TensorPtr& t)
    :
    pgrp_(t->pgrp_),
    Tile(
        IndexRangeTuple::get_unit_range(t->get_index_ranges()),
        IndexRange::get_zero_set(),
        t->config_
      ),
    nparams_(0),
    params_(0)
{
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

Tensor::~Tensor()
{
    if (params_)
        yeti_free_indexset(params_);
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
    const TensorPtr& tile,
    const PermutationPtr& p,
    double scale
)
{
    ThreadedSortPtr sort = new ThreadedSort(p);
    MatrixPtr matrix = 0;
    if (!pgrp_->contains(tile->get_permutation_grp()))
    {
        //we must turn the tile into an unpacked matrix
        usi nrows = nindex();
        usi ncols = 0;
        MatrixIndex::matrix_index_t rowtype = MatrixIndex::front;
        MatrixIndex::matrix_index_t coltype = MatrixIndex::back;
        MatrixIndex* index = new MatrixIndex(rowtype, nrows, coltype, ncols);
        MatrixConfiguration* config
            = new MatrixConfiguration(
                    index,
                    tile->get_permutation_grp(),
                    tile->get_permutation_grp()
                        ->get_identity()
               );

        PermutationGroupPtr cxngrp =get_permutation_grp();
        config->configure_quotient_set(cxngrp);
        config->configure_isotropy(cxngrp);
        matrix = new Matrix(tile, config);
        matrix->set_as_multiplicand();
        matrix->retrieve(NOT_THREADED);
        Tile::accumulate(matrix, sort, scale);
        matrix->release(NOT_THREADED);
    }
    else
    {
        if (!p->is_identity())
        {
            //validate the permutation passed in
            uli* indices = yeti_malloc_indexset();
            uli nindex = p->nindex();
            for (uli i=0; i < nindex; ++i)
            {
                indices[i] = i;
            }
            uli* permuted = yeti_malloc_indexset();
            p->permute(indices, permuted);
            if (tile->get_permutation_grp()->improves_sort(permuted))
            {
                cerr << "Invalid accumulate specified for tensor.  Permutation"
                        " violates canonical indexing." << endl;
                abort();
            }
            yeti_free_indexset(indices);
            yeti_free_indexset(permuted);
        }

        Tile::accumulate(tile, sort, scale);
    }

    config_->data_factory->allocate_blocks();
    GlobalQueue::run();

    //only clear when you reach here to ensure matrix
    //is persistent in memory for the entire contraction
    matrix = 0;
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

    if (ranges_->nindex() > 0)
        alignment_depth_ = ranges_->mindepth();
}

void
Tensor::configure(const DataBlockFactoryPtr& factory)
{
    factory->configure(this);
    config_->data_factory = factory;
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
    usi nbra = bratuple->nindex();
    usi nket = bratuple->nindex();

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
Tensor::equals(const TileElementComputerPtr &filler)
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
    if (config_->filler)
    {
        fill(config_->filler);
    }
    else
    {
        TileElementComputerPtr filler
            = new MemsetElementComputer;
        fill(filler);
    }
}

void
Tensor::fill(const ThreadedTileElementComputerPtr& new_filler)
{
    config_->filler = new_filler;

    //this will build the queue of fill tasks
    Tile::fill();

    //allocate any new blocks in the filler
    config_->data_factory->allocate_blocks();

    set_write_mode();
    GlobalQueue::run();
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
        TileElementComputerPtr filler
            = new MemsetElementComputer;
        fill(new ThreadedTileElementComputer(filler));
    }
}

IndexRange*
Tensor::get_index(usi index) const
{
    return ranges_->get(index);
}

Tensor*
Tensor::get_parameterized_tensor(usi nindex)
{
    TensorPtr split_tensor = split(nindex);

    usi ntop_index = nindex;
    IndexRangeTuplePtr top_tuple = new IndexRangeTuple(ntop_index);
    for (usi i=0; i < ntop_index; ++i)
        top_tuple->set(i, split_tensor->get_index_ranges()->get(i));

    usi nbottom_index = this->nindex() - nindex;
    IndexRangeTuplePtr bottom_tuple = new IndexRangeTuple(nbottom_index);
    for (usi i=0; i < nbottom_index; ++i)
        bottom_tuple->set(i, ranges_->get(i + ntop_index));

    //figure out the permutation groups for the "row" and "column" indices
    MatrixIndexPtr matrix_index = new MatrixIndex(
        MatrixIndex::front,
        ntop_index,
        MatrixIndex::back,
        nbottom_index
    );

    MatrixConfigurationPtr matrix_config = new MatrixConfiguration(
        matrix_index,
        pgrp_,
        pgrp_->get_identity()
    );

    PermutationGroupPtr top_pgrp = matrix_config->get_row_permutation_grp();
    PermutationGroupPtr bottom_pgrp = matrix_config->get_col_permutation_grp();

    std::string topname = config_->name + " top level";
    TensorConfigurationPtr top_config = config_->copy(topname);
    Tensor* top_tensor = new Tensor(
                            topname,
                            top_tuple,
                            top_pgrp,
                            top_config
                           );

    std::string bottom_name = config_->name + " bottom level";
    TensorConfigurationPtr bottom_config = config_->copy(bottom_name);

    usi* top_subset = yeti_malloc_perm();
    for (usi i=0; i < nindex; ++i)
        top_subset[i] = i;

    usi* bottom_subset = yeti_malloc_perm();
    for (usi i=0; i < nindex; ++i)
    {
        bottom_subset[i] = i + ntop_index;
    }

    TileIteratorPtr blockiter = split_tensor->get_iterator();
    Tile** it = blockiter->begin();
    Tile** stop = blockiter->end();
    usi split_depth = split_tensor->depth();

    uli* top_location = yeti_malloc_indexset();
    uli* bottom_location = yeti_malloc_indexset();
    top_tensor->retrieve(NOT_THREADED);

    for ( ; it != stop; ++it)
    {
        Tile* tile = *it;
        tile->get_location(top_location, top_subset, ntop_index);
        //an extra metadata layer was padded on top of this tile
        //the parent will have the exact same bottom level indices
        //get the parent tile and grab a location that removes this
        //extra layer of metadata padding
        //start the location writing at bottom_location - 1
        tile->get_parent()->get_parent_tile()
            ->get_location(bottom_location - 1, bottom_subset, nbottom_index);

        uli tile_index = top_location[0];

        TileMap* tilemap = top_tensor->get_tile_map(top_location, true); //create recursively
        Tile* child = tilemap->get(tile_index);

        if (!child) //create the child and insert it
        {
            uli* params = yeti_malloc_indexset();
            ::memcpy(params, tile->indices(), ntop_index * sizeof(uli));
            child = new Tensor(
                            bottom_name,
                            bottom_tuple,
                            bottom_pgrp,
                            bottom_config,
                            params,
                            ntop_index
                        );
            tilemap->insert(tile_index, child);
            child->set_parent(0); //clear ownership for this child
        }
        child->retrieve(NOT_THREADED);
        tilemap = child->get_tile_map(bottom_location, true); //create recursively
        tilemap->insert(bottom_location[0], tile);
        child->release(NOT_THREADED);

        //have the tile reorder itself so that the parameters no longer affect metadata structure
        tile->parameter_reduce(ntop_index);
    }
    top_tensor->release(NOT_THREADED);

    yeti_free_perm(bottom_subset);
    yeti_free_perm(top_subset);
    yeti_free_indexset(top_location);
    yeti_free_indexset(bottom_location);

    blockiter = 0;
    split_tensor = 0;

    return top_tensor;
}

void
Tensor::get_params(uli &i, uli &j) const
{
    i = params_[0];
    j = params_[1];
}

PermutationGroup*
Tensor::get_permutation_grp() const
{
    return pgrp_.get();
}

bool
Tensor::is_empty() const
{
    if (tilemap_)
    {
        foreach_nonnull(tile, tilemap_, Tile,
            return false;
        )
    }
    return true;
}

const std::string&
Tensor::name() const
{
    return config_->name;
}
usi
Tensor::nparams() const
{
    return nparams_;
}

void
Tensor::print(ostream& os) const
{
    os << Env::indent << "Tensor: " << config_->name << endl;
    if (params_)
    {
        os << Env::indent << "Parameters: ";
        print_indices(nparams_, params_);
    }

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
Tensor::_release(uli threadnum)
{
    Tile::_release(threadnum);
}

void
Tensor::_retrieve(uli threadnum)
{
    if (!config_->data_factory)
    {
        //use default memory block factory
        config_->data_factory = new MemoryBlockFactory<double>;
    }

    Tile::_retrieve(threadnum);
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

Tensor*
Tensor::split(usi nindex_split)
{
    usi nindex = ranges_->nindex();
    IndexRangeTuplePtr split_tuple = new IndexRangeTuple(nindex);
    for (usi i=0; i < nindex_split; ++i)
    {
        IndexRange* newrange = ranges_->get(i)->split_bottom_range();
        split_tuple->set(i, newrange);
    }

    for (usi i=nindex_split; i < nindex; ++i)
    {
        IndexRange* newrange = ranges_->get(i)->shift_bottom_range();
        split_tuple->set(i, newrange);
    }

    string newname = config_->name + " split index";
    TensorConfigurationPtr newconfig = config_->copy(newname);
    newconfig->data_factory = new SubsetBlockFactory(newconfig->data_factory);

    Tensor* split_tensor = new Tensor(newname, split_tuple, pgrp_, newconfig);

    Tile* nonconst_me = const_cast<Tensor*>(this);
    recurse_split(
        nonconst_me,
        split_tensor
    );

    return split_tensor;
}

void
Tensor::recurse_split(
    Tile* src,
    Tile* dst
)
{
    src->retrieve(NOT_THREADED);
    dst->retrieve(NOT_THREADED);


    usi depth = src->depth();
    if (depth == 0)
    {
        finalize_split(src, dst);
        return;
    }

    if (src->get_map()->size() != dst->get_map()->size())
    {
        cerr << "destination map must match src map" << endl;
        cerr << "Source" << endl;
        cerr << src->get_map() << endl;
        cerr << "Destination" << endl;
        cerr << dst->get_map() << endl;
        abort();
    }

    IndexRangeTuple* ranges = dst->get_index_ranges();
    uli nindex = dst->nindex();

    Tile** itsrc = src->get_map()->begin();
    Tile** itdst = dst->get_map()->begin();
    Tile** stop = src->get_map()->end();
    usi idx = 0;
    for ( ; itsrc != stop; ++itdst, ++itsrc, ++idx)
    {
        Tile* src_tile = *itsrc;
        if (!src_tile)
            continue;

        //create the destination tile
        IndexRangeTuplePtr tuple
            = new IndexRangeTuple(nindex);
        uli* newindices = yeti_malloc_indexset();
        ::memcpy(newindices, src_tile->indices(), nindex * sizeof(uli));
        for (uli i=0; i < nindex; ++i)
        {
            tuple->set(i, ranges->get(i)->get_subindex(newindices[i]));
        }

        Tile* dst_tile = new Tile(tuple, newindices, dst->config());
        dst->get_map()->insert(idx, dst_tile);

        recurse_split(src_tile, dst_tile);
    }

    src->release(NOT_THREADED);
    dst->release(NOT_THREADED);
}

void
Tensor::finalize_split(
    Tile* src,
    Tile* dst
)
{

    DataBlock* parent_data = src->get_data();
    Tile** it = dst->get_map()->begin();
    Tile** stop = dst->get_map()->end();

    //have the tile map create a dense set of tiles
    //do not filter anything based on permutational symmetry or otherwise
    dst->get_map()->allow_all_tiles();
    dst->get_map()->fill(0);

    uli offset = 0;
    for ( ; it != stop; ++it)
    {
        Tile* dst_tile = *it;
        if (!dst_tile)
            continue;

        dst_tile->retrieve(NOT_THREADED);
        SubsetDataBlock* subset = static_cast<SubsetDataBlock*>(dst_tile->get_data());
        subset->configure(parent_data, offset);
        offset += subset->data()->size();
        dst_tile->release(NOT_THREADED);
    }

}

bool
yeti::tensor_less(const TensorPtr &l, const TensorPtr &r)
{
    //distributed tensors are always greater than replicted tensors
    if (l->config()->tile_distribution_type != r->config()->tile_distribution_type)
        return l->config()->tile_distribution_type < r->config()->tile_distribution_type;

    //in core is always faster than recomputing which is probably faster than going to disk
    Data::storage_t lstorage = l->storage_type();
    Data::storage_t rstorage = r->storage_type();
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
    subdepth_(MAX_DEPTH)
{
    //figure out the recursion depth for the tile map registry
    IndexRangeTuplePtr tuple = tensor->get_index_ranges();
    for (usi i=0; i < subtuple_->nindex(); ++i)
    {
        IndexRange* parent = tuple->get(i);
        IndexRange* child = subtuple->get(i);
        usi depth = parent->get_subdepth_alignment(child);
        if (depth < subdepth_)
            subdepth_ = depth;
    }
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
    for (usi i=0; i < tuple->nindex(); ++i)
    {
        if (tuple->get(i)->depth() != subdepth_)
        {
            reuse = false;
            break;
        }
    }

    return new TileMap(tuple, tensor_->get_permutation_grp(), parent);
}

