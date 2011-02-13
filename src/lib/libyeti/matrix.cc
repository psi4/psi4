#include "matrix.h"
#include "tile.h"
#include "index.h"
#include "data.h"
#include "class.h"
#include "permutation.h"
#include "exception.h"
#include "env.h"
#include "tensor.h"
#include "contraction.h"
#include "contractionimpl.h"
#include "sort.h"
#include "cache.h"
#include "malloc.h"
#include "runtime.h"

using namespace yeti;
using namespace std;

DECLARE_MALLOC(Matrix);
DECLARE_MALLOC(MatrixMap);

Matrix::Matrix(
   const MatrixConfigurationPtr& config,
   const MatrixPtr& parent,
   const size_t rowindex,
   const size_t colindex,
   usi nisotropy
)
   :
    parent_(parent.get()),
    config_(config),
    colindex_(colindex),
    rowindex_(rowindex),
    nrows_(0),
    ncols_(0),
    data_(0),
    tiles_(0),
    perms_(0),
    ngenerator_(0),
    nisotropy_(nisotropy),
    mtype_(none),
    unique_tile_(false),
    tmpindices_(yeti_malloc_indexset()),
    constructed_(false),
    maxlog_(0),
    tuple_(0)
{
    set_as_multiplicand();
}

Matrix::Matrix(
   const TensorPtr& tensor,
   const MatrixConfigurationPtr& config
)
   :
    parent_(0), //null parent, top-level tile
    Imap_(0),
    config_(config),
    colindex_(0),
    rowindex_(0),
    nrows_(0),
    ncols_(0),
    data_(0),
    nisotropy_(1),
    mtype_(none),
    ngenerator_(0),
    perms_(0),
    tiles_(0),
    unique_tile_(false),
    tmpindices_(yeti_malloc_indexset()),
    constructed_(false),
    maxlog_(tensor->max_log()),
    tuple_(0)
{
    //this is a top level matrix so we must allocate the generators in construction
    allocate_generators(config->get_quotient_set()->order());

    //loop the quotient set and add the parent tile with every permutation
    PermutationSet::iterator it(config->get_quotient_set()->begin());
    PermutationSet::iterator stop(config->get_quotient_set()->end());
    for ( ; it != stop; ++it)
        add_generator(tensor.get(), (*it).get());

    set_as_multiplicand();
}

Matrix::Matrix(
    const MatrixPtr& lmatrix,
    const MatrixPtr& rmatrix,
    const ContractionPtr& cxn
)
   :
    parent_(0), //top level matrix
    config_(cxn->get_product_matrix_configuration()),
    colindex_(0),
    rowindex_(0),
    nrows_(lmatrix->ncols()),
    ncols_(rmatrix->ncols()),
    data_(0),
    tiles_(0),
    perms_(0),
    ngenerator_(0),
    mtype_(none),
    nisotropy_(1), //top level, so prefactor is always 1
    unique_tile_(false),
    tmpindices_(yeti_malloc_indexset()),
    Imap_(0),
    constructed_(false),
    maxlog_(0),
    tuple_(0)
{
    //add generators for each permutation in the product set
    //these are the permutations that must be applied to the matrix blocks
    //in order to give the tensor the correct permutational symmetry
    PermutationSet* set = cxn->get_product_set();

    //very important here to get the cxn product tensor which
    //has the appropriate extra metadata "padding"
    Tensor* tensor = cxn->get_product_tensor();

    //this is a top level matrix so we must allocate the generators in construction
    allocate_generators(set->order());

    PermutationSet::iterator it(set->begin());
    PermutationSet::iterator stop(set->end());
    for ( ; it != stop; ++it)
    {
        add_generator(tensor, (*it).get());
    }

    accumulate_product(
        lmatrix.get(),
        rmatrix.get(),
        cxn.get(),
        NOT_THREADED
    );

    //register any new allocations
    tensor->config()->get_data_factory()->allocate_blocks();
}

Matrix::~Matrix()
{
    ++Env::indent;
    //it's very important that we delete the map and the data in the corre
    Imap_ = 0;
    data_ = 0;

    --Env::indent;

    if (unique_tile_) //data tile
    {
        yeti_free_matrix_generator(perms_);
        yeti_free_matrix_generator(tiles_);
    }
    else
    {
        yeti_free_matrix_generator_set(perms_);
        yeti_free_matrix_generator_set(tiles_);
    }

    yeti_free_indexset(tmpindices_);


}

#if 0
void
Matrix::accumulate(
    const MatrixPtr &matrix,
    double scale
)
{
#if YETI_SANITY_CHECK
    if (matrix->nrows() != nrows() ||
        matrix->ncols() != ncols() ||
        matrix->depth() != depth())
        raise(SanityCheckError, "matrices are not aligned for accumulation");
#endif
    if (Imap_) //metadata
    {
        accumulate_metadata(matrix, scale);
    }
    else
    {
        accumulate_data(matrix, scale);
    }
}
#endif

#if 0
void
Matrix::accumulate_metadata(
    const MatrixPtr &matrix,
    double scale
)
{
    MatrixBlockMap::iterator itl(Imap_->begin());
    MatrixBlockMap::iterator itr(matrix->get_row_map()->begin());
    MatrixBlockMap::iterator lstop(Imap_->end());
    MatrixBlockMap::iterator rstop(matrix->get_row_map()->end());

    while(     itl != lstop
            && itr != rstop )

    for (; itl != lstop; ++itl, ++itr)
    {
        MatrixRow* r_row(*itr);
        MatrixRow* l_row(*itl);
        if (!r_row) //nothing to accumulate
            continue;

        //now that we have compatible rows, we can accumulate
        accumulate_row(
            l_row,
            r_row,
            scale
        );
    }
}

void
Matrix::accumulate_data(
    const MatrixPtr &matrix,
    double scale
)
{
    DataBlock* lblock = matrix->get_data();
    DataBlock* rblock = matrix->get_data();

    lblock->retrieve();
    rblock->retrieve();

    Data* ldata = lblock->data();
    Data* rdata = rblock->data();
    ldata->accumulate(rdata, scale, 0);


    lblock->release();
    rblock->release();
}


void
Matrix::accumulate_row(
    MatrixRow* l_row,
    MatrixRow* r_row,
    double scale
)
{
    Matrix** itl(l_row->begin());
    Matrix** lstop(l_row->end());
    Matrix** itr(r_row->begin());
    Matrix** rstop(r_row->end());
    for ( ; itl != lstop; ++itl, ++itr)
    {
        Matrix* l(*itl);
        Matrix* r(*itr);
        if (!r) //nothing to accumulate
            continue;

        l->retrieve();
        r->retrieve();
        l->accumulate(r, scale);
    }
}
#endif

void
Matrix::accumulate_data_product(
    Matrix* lmatrix,
    Matrix* rmatrix,
    Contraction* cxn,
    uli threadnum
)
{
    DataBlock* ldata = lmatrix->get_data();
    DataBlock* rdata = rmatrix->get_data();

    ldata->retrieve(threadnum);
    rdata->retrieve(threadnum);
    data_->retrieve(threadnum);

    data_->lock(); //only accumulate to one block at a time

    //determine the type for the contraction
    TemplateInfo::type_t ltype = ldata->data()->type();
    TemplateInfo::type_t rtype = rdata->data()->type();

    //select the larger of the two values
    //integers need to go up in precision
    //doubles can go down in precision
    TemplateInfo::type_t cxn_data_type =
        ltype > rtype ? ltype : rtype;

    uli nrows = lmatrix->ncols();
    uli ncols = rmatrix->ncols();
    uli nlink = rmatrix->nrows();

#if YETI_SANITY_CHECK
    uli lmat_nrows = lmatrix->nrows();
    if (lmat_nrows != nlink)
    {
        cerr << lmatrix->get_unique_tile() << endl;
        cerr << rmatrix->get_unique_tile() << endl;
        raise(SanityCheckError, "matrices not aligned for multiplication");
    }
#endif

    usi l_niso = lmatrix->nisotropy();
    usi r_niso = rmatrix->nisotropy();

#if YETI_SANITY_CHECK
    if (l_niso != r_niso)
    {
        cout << lmatrix << endl;
        cout << rmatrix << endl;
        throw ValuesNotEqual<usi>("nisotropy",l_niso, r_niso,__FILE__, __LINE__);
    }
#endif

    ContractionEngine* engine = cxn->get_engine();

    engine->contract(
        ldata->data(),
        rdata->data(),
        data_->data(),
        lmatrix->ncols(),
        rmatrix->ncols(),
        lmatrix->nrows(),
        cxn->left_tmp_space(),
        cxn->right_tmp_space(),
        cxn->product_tmp_space(),
        cxn->scale_factor() * l_niso,
        cxn_data_type
    );

    data_->unlock();

    ldata->release(threadnum);
    rdata->release(threadnum);
    data_->release(threadnum);
}

void
Matrix::accumulate_product(
    Matrix* lmatrix,
    Matrix* rmatrix,
    Contraction* cxn,
    uli threadnum
)
{
    this->set_as_product();
    lmatrix->set_as_multiplicand();
    rmatrix->set_as_multiplicand();

    this->retrieve(threadnum);
    lmatrix->retrieve(threadnum);
    rmatrix->retrieve(threadnum);

    MatrixMap* lmap(lmatrix->get_map());
    MatrixMap* rmap(rmatrix->get_map());

    bool l_is_metadata = lmap;

#if YETI_SANITY_CHECK
    bool r_is_metadata = rmap;
    bool p_is_metadata = Imap_;

    if (l_is_metadata != r_is_metadata)
        raise(SanityCheckError, "meta data depths not aligned for multiplication");

    if (lmatrix->get_unique_tile()->config()->get_data_mode()->flag != DataMode::read)
        raise(SanityCheckerror, "left matrix not configured for read");
    if (rmatrix->get_unique_tile()->config()->get_data_mode()->flag != DataMode::read)
        raise(SanityCheckerror, "right matrix not configured for read");
    if (get_unique_tile()->config()->get_data_mode()->flag != DataMode::accumulate)
    {
        cerr << "Data mode = " << (void*) get_unique_tile()->config()->get_data_mode() << endl;
        raise(SanityCheckerror, "product matrix not configured for accumulate");
    }
#endif

    if (!l_is_metadata) //no metadata work to do
        accumulate_data_product(lmatrix, rmatrix, cxn, threadnum);
    else
        accumulate_metadata_product(lmatrix, rmatrix, cxn, threadnum);

    this->release(threadnum);
    rmatrix->release(threadnum);
    lmatrix->release(threadnum);
}

void
Matrix::accumulate_metadata_product(
    Matrix *lmatrix,
    Matrix *rmatrix,
    Contraction *cxn,
    uli threadnum
)
{
    uli nlink = lmatrix->ncxn_rows();
    uli nrows = lmatrix->ncxn_cols();
    uli ncols = rmatrix->ncxn_cols();

    Matrix** ldata = lmatrix->get_map()->begin();
    Matrix** rdata = rmatrix->get_map()->begin();

    usi ldepth = lmatrix->depth();
    usi rdepth = rmatrix->depth();
    usi mydepth = depth();

    if (mydepth < ldepth && mydepth < rdepth) //don't iterate this matrix
    {
        //the number of rows and columns should only be 1
        if (nrows != 1 && ncols != 1)
        {
            raise(SanityCheckError, "matrix of lower depth should have only 1 row and 1 col");
        }

        //this matrix should accumulate the contributions
        for (uli link=0; link < nlink; ++link, ++ldata, ++rdata)
        {
            Matrix* l = *ldata;
            Matrix* r = *rdata;

            if (!l || !r) //null matrices, nothing to accumulate
                continue;

            float log_max_product = l->max_log() + r->max_log();
            if (log_max_product < YetiRuntime::matrix_multiply_cutoff)
                continue;

            if (mydepth == 0) //data task
            {
                cxn->add_task(l, r, this);
                Tile* tile = get_unique_tile();
                if (tile->get_data()->data()->null())
                {
                    tile->config()->get_data_factory()
                        ->register_allocation(tile->get_data());
                }
            }
            else
            {
                accumulate_product(l, r, cxn, threadnum);
            }
        }
        return;
    }


#if YETI_SANITY_CHECK
    uli nrows_me = this->nrows();
    uli ncols_me = this->ncols();
    Tile* mytile = get_unique_tile();
    if (nrows_me != nrows || ncols_me != ncols)
    {
        Tile* ltile = lmatrix->get_unique_tile();
        cout << ClassOutput<const uli*>::str(ltile->nindex(), ltile->index_sizes()) << endl;
        ltile->config()->print_data = false;
        cout << ltile << endl;

        Tile* rtile = rmatrix->get_unique_tile();
        cout << ClassOutput<const uli*>::str(rtile->nindex(), rtile->index_sizes()) << endl;
        rtile->config()->print_data = false;
        cout << rtile << endl;

        cout << ClassOutput<const uli*>::str(mytile->nindex(), mytile->index_sizes()) << endl;
        mytile->config()->print_data = true;
        cout << mytile << endl;

        abort();
    }
#endif

    incref();
    //now accumulate the contraction

    TileMap* tilemap = get_unique_tile()->get_map();
    Matrix** pdata = get_map()->begin();
    for (uli link=0; link < nlink; ++link, ldata += nrows, rdata += ncols)
    {
        Matrix** lptr = ldata;
        Matrix** pptr = pdata;
        uli compidx = 0;
        for (uli row=0; row < nrows; ++row, ++lptr)
        {
            Matrix** rptr = rdata;
            for (uli col=0; col < ncols; ++col, ++rptr, ++pptr, ++compidx)
            {
                Matrix* l = *lptr;
                Matrix* r = *rptr;
                Matrix* p = *pptr;

                if (!l || !r) //null matrices, nothing to accumulate
                    continue;

                if (!p) //block doesnt exist, attempt to construct it
                {

                    //not rigorously zero, so attempt to build the subblock
                    p = build_subproduct(
                            compidx,
                            row,
                            col,
                            l,
                            r,
                            cxn
                        );

                    //this product matrix is not meant to be built
                    if (!p)
                        continue;

                    p->incref(); //make sure to increment the refcount to ensure non-deletion
                    (*pptr) = p;
                }

                //we have non-null product and non-null multiplicands
                if (ldepth == DATA_DISTRIBUTION_DEPTH)
                {
                    cxn->add_task(l, r, p);
                }
                else
                {
                    p->accumulate_product(l, r, cxn, threadnum);
                }
            }
        }
    }
    decref();
}

Matrix*
Matrix::build_subproduct(
    uli composite_index,
    uli rowindex,
    uli colindex,
    Matrix* lmatrix,
    Matrix* rmatrix,
    Contraction* cxn
)
{
    float log_max_product = lmatrix->max_log() + rmatrix->max_log();

    if (log_max_product < YetiRuntime::matrix_multiply_cutoff)
    {
        return 0;
    }

    Tile* tile = *tiles_;
    TileMap* tilemap = tile->get_map();
    Permutation* p = *perms_;

    size_t* indexset = yeti_malloc_indexset();
    tilemap->compute_indices(composite_index, indexset);
    usi nindex = tile->nindex();

    //from the permutation group, determine if we need to make this block
    if (tilemap->depths_aligned())
    {
        bool unique =
            is_unique_block(
              indexset,
              0,
              cxn->get_default_grp(),
              tile->depth(),
              0 //existence of matrix doesn't matter
           );
        if (!unique)
            return 0;
    }

    Permutation* plowest = PermutationSet::lowest_image(
                perms_, ngenerator_,
                indexset, tmpindices_,
                tile->nindex()
           );

    uli unique_index = tilemap->index(tmpindices_);
    bool is_zero = tilemap->is_rigorously_zero(unique_index);
    if (is_zero) //not to be built
        return 0;

    usi nisotropy = 1; //always 1 for a product matrix
    Matrix* pmatrix = new Matrix(
                        config_,
                        this,
                        rowindex,
                        colindex,
                        nisotropy
                      );

    Tile* newtile = tilemap->get(unique_index);
    if (!newtile)
    {
        IndexRangeTuplePtr tuple = new IndexRangeTuple(nindex);
        for (usi i=0; i < nindex; ++i)
        {
            IndexRange* subindex = tilemap->get_subrange(i, indexset[i]);
            tuple->set(i, subindex);
        }

        //create the tile
        tuple->permute(plowest);
        uli* indexarr = yeti_malloc_indexset();
        ::memcpy(indexarr, tmpindices_, nindex * sizeof(uli));
        newtile = new Tile(
            tuple,
            indexarr,
            tile->get_map()
        );
        newtile->retrieve(NOT_THREADED);
        newtile->release(NOT_THREADED);
        tilemap->insert(newtile);

        //if this is a data tile, register its data block for allocation
        if (tuple->maxdepth() == 0) //data tile
        {
            tile->config()->get_data_factory()
                ->register_allocation(newtile->get_data());
        }
    }

    for (usi i=0; i < ngenerator_; ++i)
    {
        Permutation* pgen = perms_[i];

        //check to see if the permutation maps onto the unique tile
        bool maps_to_tile = pgen->maps_to(indexset, tmpindices_);
        if (maps_to_tile)
        {
            pmatrix->add_generator(newtile, pgen);
        }
    }

    yeti_free_indexset(indexset);

    return pmatrix;
}

void
Matrix::allocate_generators(usi ngen)
{
    if (ngen == 1)
    {
        tiles_ = ptrcast(Tile*, yeti_malloc_matrix_generator());
        perms_ = ptrcast(Permutation*, yeti_malloc_matrix_generator());
        unique_tile_ = true;
    }
    else
    {
        tiles_ = ptrcast(Tile*, yeti_malloc_matrix_generator_set());
        perms_ = ptrcast(Permutation*, yeti_malloc_matrix_generator_set());
        unique_tile_ = false;
    }

}

void
Matrix::add_generator(
    Tile* tile,
    Permutation* p
)
{
    if (tiles_ == 0) //not yet allocated
    {
        allocate_generators(parent_->ngenerator());
        maxlog_ = tile->max_log();
    }
    else
    {
        float max = tile->max_log();
        if (max > maxlog_)
            maxlog_ = max;
    }

    if (!tuple_) set_tuple(tile, p);

    tiles_[ngenerator_] = tile;
    perms_[ngenerator_] = p;

    if (ngenerator_ == 0) //dimensions are not yet computed
        compute_dimensions();

    ++ngenerator_;
}

size_t
Matrix::colindex() const
{
    return colindex_;
}

void
Matrix::compute_dimensions()
{

    Tile* tile = *tiles_;
    tile->retrieve(NOT_THREADED);

    IndexRangeTuple* tuple = tile->get_index_ranges();

    const usi* pmap = (*perms_)->indexmap();
    IndexRange** ranges = tile->get_index_ranges()->begin();

    const usi* iter = config_->get_index()->cols();
    usi n = config_->get_index()->ncolindex();

    uli ndim = 1;
    for (usi i=0; i < n; ++i, ++iter)
    {
        //ndim *= ranges[pmap[*iter]]->n();
        ndim *= tile->nrange(pmap[*iter]);
    }
    ncols_ = ndim;

    iter = config_->get_index()->rows();
    n = config_->get_index()->nrowindex();
    ndim = 1;
    for (usi i=0; i < n; ++i, ++iter)
    {
        //ndim *= ranges[pmap[*iter]]->n();
        ndim *= tile->nrange(pmap[*iter]);
    }
    nrows_ = ndim;

    tile->release(NOT_THREADED);
}

bool
Matrix::data_equals(const MatrixPtr& matrix)
{
    DataBlock* ldata(this->get_data());
    DataBlock* rdata(matrix->get_data());

    ldata->retrieve(NOT_THREADED);
    rdata->retrieve(NOT_THREADED);

    if (ldata->data()->size() != rdata->data()->size())
    {
        cerr << "data blocks not of same size" << endl;
        cerr << ldata << endl;
        cerr << rdata << endl;
        ldata->release(NOT_THREADED);
        rdata->release(NOT_THREADED);
        return false;
    }

    bool eq = ldata->data()->equals(rdata->data()->pointer());

    ldata->release(NOT_THREADED);
    rdata->release(NOT_THREADED);

    return eq;
}

usi
Matrix::depth() const
{
    return get_unique_tile()->depth();
}

bool
Matrix::equals(const MatrixPtr& matrix)
{
    if (this->rowindex() != matrix->rowindex())
    {
        cerr << "rowindex not equals" << endl;
        print(cerr); cerr << endl;
        cerr << matrix << endl;
        return false;
    }

    if (this->colindex() != matrix->colindex())
    {
        cerr << "colindex not equals" << endl;
        print(cerr); cerr << endl;
        cerr << matrix << endl;
        return false;
    }

    if (this->nrows() != matrix->nrows())
    {
        cerr << "ncols not equal" << endl;
        print(cerr); cerr << endl;
        cerr << matrix << endl;
        return false;
    }

    if (this->ncols() != matrix->ncols())
    {
        cerr << "ncols not equal" << endl;
        print(cerr); cerr << endl;
        cerr << matrix << endl;
        return false;
    }

    this->retrieve_as_multiplicand(0); //thread 0
    matrix->retrieve_as_multiplicand(0); //thread 0

    MatrixMap* lmap(get_map());
    MatrixMap* rmap(matrix->get_map());

    bool l_is_metadata = lmap;
    bool r_is_metadata = rmap;
    if (l_is_metadata != r_is_metadata)
        raise(SanityCheckError, "meta data depths not aligned for multiplication");

    if (!l_is_metadata) //no metadata work to do
    {
        return data_equals(matrix);
    }
    else
    {
        return metadata_equals(matrix);
    }

    this->release_as_multiplicand(0); //thread 0
    matrix->release_as_multiplicand(0);
}

Matrix*
Matrix::get(size_t row, size_t col) const
{
    return Imap_->get(row, col);
}

DataBlock*
Matrix::get_data() const
{
    return data_.get();
}

IndexRangeTuple*
Matrix::get_index_ranges() const
{
    return tuple_.get();
}

MatrixConfiguration*
Matrix::get_matrix_config() const
{
    return config_.get();
}

MatrixMap*
Matrix::get_map() const
{
    return Imap_.get();
}

Permutation*
Matrix::get_unique_permutation() const
{
    return *perms_;
}

Tile*
Matrix::get_unique_tile() const
{
    return *tiles_;
}

bool
Matrix::is_unique_block(
    const size_t *indices,
    Permutation* p,
    PermutationGroup *grp,
    usi depth,
    Matrix* matrix
)
{
    if (depth == 1 && matrix)
        return false; //only one data block should be generated

    if (p && !p->is_identity())
    {
        p->permute(indices, tmpindices_);
        return Tile::is_unique(tmpindices_, grp, depth);
    }
    else
    {
        return true; //is unique
    }
}

void
Matrix::make_blocks()
{
    if (Imap_)
        return; //already built

    for (usi idx=0; idx < ngenerator_; ++idx)
        make_blocks(tiles_[idx], perms_[idx]);
}

void
Matrix::make_blocks(
    Tile* parent,
    Permutation* p
)
{
    TileMap* tilemap(parent->get_map());
    if (!tilemap)
        raise(SanityCheckError, "cannot make blocks without a tile map");

    if (!Imap_)
        Imap_ = new MatrixMap(ncxn_rows(), ncxn_cols());

    usi depth = tilemap->maxdepth();

    incref(); //passing this pointer... incref to ensure non-deletion

    uli nrows = ncxn_rows();
    uli ncols = ncxn_cols();

    Permutation* rowcolperm = config_->get_index()->get_row_column_permutation();
    PermutationPtr product = rowcolperm->product(p);

    Tile** tiles = new Tile*[nrows * ncols];
    SortPtr sort = new Sort(product, tilemap->nindices());
    sort->sort_noscale<Tile*>(tilemap->begin(), tiles);

    Tile** tileptr = tilemap->begin();
    for (uli row=0; row < nrows; ++row)
    {
        for (uli col=0; col < ncols; ++col, ++tileptr)
        {
            Tile* tile = *tileptr;
            if (!tile)
                continue;

        }
    }

    if (nrows * ncols != sort->ntot())
    {
        cerr << "sort is not large enough" << endl;
        cerr << "nr=" << nrows << endl;
        cerr << "nc=" << ncols << endl;
        //cerr << get_unique_tile()->get_index_ranges() << endl;
        cerr << get_unique_tile()->get_map() << endl;
        abort();
    }

    tileptr = tiles;
    for (uli row=0; row < nrows; ++row)
    {
        for (uli col=0; col < ncols; ++col, ++tileptr)
        {
            Tile* tile = *tileptr;
            if (!tile) //no work to be done for a null tile
                continue;

            uli rowindex = row;
            uli colindex = col;

            const uli* indices = tile->indices();

            bool unique = false;

            Matrix* matrix(get(rowindex, colindex));
            unique = is_unique_block(
                    indices,
                    p,
                    config_->get_matrix_permutation_grp().get(),
                    depth,
                    matrix
                  );

            if (!unique)
                continue;

            //if matrix doesn't exist yet, build it
            if (!matrix)
            {
                matrix = new Matrix(
                        config_,
                        this,
                        rowindex,
                        colindex,
                        config_->nisotropy(tile, p)
                        );

                //must add generator first to configure matrix properly
                if (depth == 1) //unique generator for data blocks
                    matrix->set_generator(tile, p);
                else
                    matrix->add_generator(tile, p);

                //then call append once configured
                Imap_->append(matrix);
            }
            else
            {
                if (depth == 1)
                    matrix->set_generator(tile, p);
                else
                    matrix->add_generator(tile, p);
            }
        }
    }
    delete[] tiles;

    decref();
}

float
Matrix::max_log() const
{
    return maxlog_;
}

bool
Matrix::metadata_equals(const MatrixPtr& matrix)
{
    MatrixMap* lmap(this->get_map());
    MatrixMap* rmap(matrix->get_map());
    Matrix** itl(lmap->begin());
    Matrix** itr(rmap->begin());
    Matrix** stop(lmap->end());
    for ( ; itl != stop; ++itl, ++itr)
    {
        Matrix* l(*itl);
        Matrix* r(*itr);

        if ( (!l && r) || (!r && l) )
            return false; //mismatched null structure

        if (!l && !r) //both null
            continue;

        l->retrieve(NOT_THREADED);
        r->retrieve(NOT_THREADED);

        bool eq = l->equals(r);

        l->release(NOT_THREADED);
        r->release(NOT_THREADED);

        if (!eq)
            return false;
    }

    //all blocks equal
    return true;
}

void
Matrix::metadata_retrieve(usi depth)
{

    make_blocks();
    foreach_nonnull(m, Imap_, Matrix,
        m->retrieve_to_depth(depth);
    )
}

uli
Matrix::ncols() const
{
    return ncols_;
}

uli
Matrix::ncol_data_blocks()
{
    if (Imap_)
    {
        uli ntot = 0;
        Matrix** mptr = Imap_->begin();
        uli nr = nrows();
        uli nc = ncols();
        for (uli col=0; col < nc; ++col, ++mptr)
        {
            Matrix* m = *mptr;
            m->retrieve(NOT_THREADED);
            ntot += m->ncol_data_blocks();
            m->release(NOT_THREADED);
        }
        return ntot;
    }
    else
    {
        return 1; //this is a data block
    }
}

uli
Matrix::nrow_data_blocks()
{
    if (Imap_)
    {
        uli ntot = 0;
        uli nr = nrows();
        uli nc = ncols();
        Matrix** mptr = Imap_->begin();
        for (uli row=0; row < nr; ++row, mptr += nc)
        {
            Matrix* m = *mptr;
            m->retrieve(NOT_THREADED);
            ntot += m->ncol_data_blocks();
            m->release(NOT_THREADED);
        }
        return ntot;
    }
    else
    {
        return 1; //this is a data block
    }
}

uli
Matrix::ncxn_cols() const
{
    if (ncols_ == 0)
    {
        cerr << "no columns!" << endl;
        abort();
    }
    return ncols_ == 0 ? 1 : ncols_;
}

uli
Matrix::ncxn_rows() const
{
    return nrows_ == 0 ? 1 : nrows_;
}

usi
Matrix::ngenerator() const
{
    return ngenerator_;
}

usi
Matrix::nisotropy() const
{
    return nisotropy_;
}

size_t
Matrix::nrows() const
{
    return nrows_;
}

void
Matrix::print(ostream& os) const
{
    if (!parent_) //this is a top-level tile
        os << Env::indent << "Matrix" << endl;

    os << Env::indent << "M(" << this->rowindex() << "," << this->colindex() << ")"
       << " n=" << nisotropy_
       << " max=10^" << maxlog_;

    uli* permuted = yeti_malloc_indexset();
    for (usi idx=0; idx < ngenerator_; ++idx)
    {
        Tile* tile = tiles_[idx];
        perms_[idx]->permute(tile->indices(), permuted);
        os << endl << Env::indent << ClassOutput<const size_t*>::str(tile->nindex(), permuted);
        os << "  ";
        perms_[idx]->print(os);
    }
    yeti_free_indexset(permuted);

    if (!Imap_)
    {
        TilePtr tile(get_unique_tile());
        if (!tile->config()->print_data || !data_)
            return; //don't print data

        //no submatrices, print the data, maybe
        data_->retrieve(NOT_THREADED);
        os << endl;
        data_->print(os);
        data_->release(NOT_THREADED);
        return;
    }

    ++Env::indent;
    foreach_nonnull(m, Imap_, Matrix,
        os << endl;
        m->retrieve(NOT_THREADED);
        m->print(os);
        m->release(NOT_THREADED);
    )
    --Env::indent;
}

#if 0
Matrix*
Matrix::recursive_get(uli rowindex, uli colindex) const
{
    if (Imap_)
        return Imap_->get(rowindex, colindex);
    else
        return parent_->recursive_get(rowindex, colindex);
}
#endif

void
Matrix::_release(uli threadnum)
{
    //let the tiles know we no longer need them
    for (usi idx=0; idx < ngenerator_; ++idx)
        tiles_[idx]->release(threadnum);

    if (mtype_ == multiplicand)
        release_as_multiplicand(threadnum);
    else if (mtype_ == product)
        release_as_product(threadnum);
    else
    {
        cerr << "cannot release matrix" << endl;
        abort();
    }
}

void
Matrix::release_as_multiplicand(uli threadnum)
{
}

void
Matrix::release_as_product(uli threadnum)
{
}

void
Matrix::retrieve_to_depth(usi depth)
{
    Tile* tile(get_unique_tile());
    TileMap* tilemap(tile->get_map());
    if (tilemap && tilemap->maxdepth() == depth)
        return; //don't build further

    DataBlock* dblock(tile->get_data());
    if (dblock) //meta data tile
        retrieve_as_data_multiplicand(dblock, 0); //thread 0
    else
        metadata_retrieve(depth);
}

void
Matrix::_retrieve(uli threadnum)
{
    //make sure all tiles are available throughout use
    for (usi i=0; i < ngenerator_; ++i)
        tiles_[i]->retrieve(threadnum);

    if (constructed_)
        return;

    if (mtype_ == multiplicand)
        retrieve_as_multiplicand(threadnum);
    else if (mtype_ == product)
        retrieve_as_product(threadnum);
    else
    {
        cerr << "cannot retrieve matrix" << endl;
        abort();
    }

    constructed_ = true;
}

void
Matrix::retrieve_as_data_multiplicand(DataBlock* dblock, uli threadnum)
{
    if (data_)
        return; //already retrieved

#if YETI_SANITY_CHECK
    verify_unique_generator();
#endif

    Permutation* p = *perms_; //(get_unique_permutation());
    Tile* tile = *tiles_; //get_unique_tile());

    if (p->is_identity()) //no work to do on id permutation
    {
        data_ = dblock;
    }
    else
    {
        LayeredDataCachePtr cache(tile->config()->cache);
        data_ = new SortedBlock(
                        tile->config()->get_data_mode(),
                        dblock,
                        cache,
                        tile->get_index_ranges(),
                        perms_,
                        1 //only a single permutation matters
                     );
    }
}

void
Matrix::retrieve_as_metadata_multiplicand(uli threadnum)
{
    make_blocks();
}

void
Matrix::retrieve_as_multiplicand(uli threadnum)
{
    Tile* tile = *tiles_;

    DataBlock* dblock(tile->get_data());
    if (dblock) //meta data tile
        retrieve_as_data_multiplicand(dblock, threadnum);
    else
        retrieve_as_metadata_multiplicand(threadnum);
}

void
Matrix::retrieve_as_product(uli threadnum)
{
    Tile* tile = *tiles_;

    DataBlock* dblock(tile->get_data());
    if (dblock)
        retrieve_as_data_product(dblock, threadnum);
    else
        retrieve_as_metadata_product(threadnum);
}

void
Matrix::retrieve_as_data_product(DataBlock* dblock, uli threadnum)
{
    //we can accumulate directly to the tile
    if (ngenerator_ == 1 && (*perms_)->is_identity())
    {
        data_ = dblock;
        return;
    }

    //we must create sorted accumulate blocks
    Tile* tile = *tiles_;
    data_ = new SortedBlock(
        tile->config()->get_data_mode(),
        dblock,
        tile->config()->cache,
        tile->get_index_ranges(),
        perms_,
        ngenerator_
    );
}

void
Matrix::retrieve_as_metadata_product(uli threadnum)
{
    if (!Imap_)
    {
        Imap_ = new MatrixMap(ncxn_rows(), ncxn_cols());
    }
}

size_t
Matrix::rowindex() const
{
    return rowindex_;
}

void
Matrix::set_as_product()
{
    mtype_ = product;
}

void
Matrix::set_as_multiplicand()
{
    mtype_ = multiplicand;
}

void
Matrix::set_tuple(
    Tile* tile,
    Permutation* p
)
{
    if (p->is_identity())
    {
        tuple_ = tile->get_index_ranges();
    }
    else
    {
        tuple_ = tile->get_index_ranges()->copy(p);
    }
}

void
Matrix::set_generator(
    Tile* tile,
    Permutation* p
)
{
    if (tiles_ == 0)
    {
        allocate_generators(1);
        maxlog_ = tile->max_log();
    }

    tiles_[0] = tile;
    perms_[0] = p;

    if (!tuple_) set_tuple(tile, p);

    if (ngenerator_ == 0) //dimensions are not yet computed
        compute_dimensions();
    else
        raise(SanityCheckError, "generator already set");

    ngenerator_ = 1;
}

void
Matrix::verify_unique_generator()
{
    //data tile!
    if (ngenerator_ > 1)
    {
        for (usi idx=0; idx < ngenerator_; ++idx)
        {
            perms_[idx]->print(cerr); cerr << endl;
        }
        raise(SanityCheckError, "quotient set for data matrix must have only 1 element");
    }
}

MatrixConfiguration::MatrixConfiguration(
    const MatrixIndexPtr& index,
    const PermutationGroupPtr& pgrp,
    const PermutationPtr& indexperm
) : index_(index),
    pgrp_(pgrp),
    rowgrp_(0),
    colgrp_(0),
    fullrowgrp_(0),
    fullcolgrp_(0),
    quotientset_(0),
    cxngrp_(0),
    matrix_grp_(0),
    workspace_(yeti_malloc_indexset()),
    tilegrp_(pgrp),
    indexperm_(indexperm)
{
    //configure everything here...
    if (indexperm->is_identity() || pgrp_->order() == 1)
        pgrp_ = pgrp;
    else
        pgrp_ = pgrp->conjugate(indexperm);

    fullrowgrp_ = pgrp_->subgrp(index->rows(), index->nrowindex());
    fullcolgrp_ = pgrp_->subgrp(index->cols(), index->ncolindex());

    rowgrp_ = pgrp_->compressed_subgrp(index->rows(), index->nrowindex());
    colgrp_ = pgrp_->compressed_subgrp(index->cols(), index->ncolindex());

    if (!indexperm_) //null passed in
        indexperm_ = pgrp_->get_identity();
}

MatrixConfiguration::~MatrixConfiguration()
{
    yeti_free_indexset(workspace_);
}

Indexer*
MatrixConfiguration::colindex(
    const IndexRangeTuplePtr& tuple,
    const PermutationPtr& p
) const
{
    return new Indexer(tuple, index_->cols(), index_->ncolindex(), p);
}

Indexer*
MatrixConfiguration::colindex(
    const uli* sizes,
    const uli* offsets,
    const PermutationPtr& p
) const
{
    return new Indexer(sizes, offsets, index_->cols(), index_->ncolindex(), p);
}

void
MatrixConfiguration::configure_isotropy(
    const PermutationGroupPtr &grp
)
{
    cxngrp_ = grp;
    matrix_grp_ = cxngrp_->union_grp(this->fullcolgrp_);
}

void
MatrixConfiguration::configure_quotient_set(
    const PermutationGroupPtr& cxn_grp
)
{
    PermutationGroupPtr target_grp(cxn_grp->union_grp(fullcolgrp_));
    quotientset_ = tilegrp_->quotient_set(target_grp);
    if (!indexperm_->is_identity())
    {
        quotientset_ = quotientset_->orbit(indexperm_);
    }
}

PermutationGroupPtr
MatrixConfiguration::get_col_permutation_grp() const
{
    return colgrp_;
}

PermutationGroupPtr
MatrixConfiguration::get_full_col_permutation_grp() const
{
    return fullcolgrp_;
}

PermutationGroupPtr
MatrixConfiguration::get_full_permutation_grp() const
{
    return pgrp_;
}

PermutationGroupPtr
MatrixConfiguration::get_full_row_permutation_grp() const
{
    return fullrowgrp_;
}

MatrixIndex*
MatrixConfiguration::get_index() const
{
    return index_.get();
}

PermutationGroupPtr
MatrixConfiguration::get_matrix_permutation_grp() const
{
    return matrix_grp_;
}

PermutationSetPtr
MatrixConfiguration::get_quotient_set() const
{
    return quotientset_;
}

PermutationGroupPtr
MatrixConfiguration::get_row_permutation_grp() const
{
    return rowgrp_;
}

usi
MatrixConfiguration::nisotropy(
    const TilePtr &tile,
    const PermutationPtr& p
)
{
    if (cxngrp_->order() == 1)
        return 1;

    const size_t* indices = tile->indices();//->data();
    p->permute(indices, workspace_);
    int niso = 0;
    PermutationGroup::iterator it(cxngrp_->begin());
    PermutationGroup::iterator stop(cxngrp_->end());
    for ( ; it != stop; ++it)
    {
        if ( (*it)->fixes_arr(workspace_) )
            ++niso;
    }
    niso = cxngrp_->order() / niso;
    return niso;
}

void
MatrixConfiguration::print(std::ostream &os) const
{
    Serializable::print(os);
}

void
MatrixConfiguration::product(
    const MatrixIndex* lindex,
    const MatrixIndex* rindex,
    const IndexRangeTuple* ltuple,
    const IndexRangeTuple* rtuple,
    const size_t* lindexset,
    const size_t* rindexset,
    const Permutation* lperm,
    const Permutation* rperm,
    size_t* indexset,
    IndexRangeTuplePtr& tuple
)
{
    //build the composite index set
    usi nrowidx = lindex->ncolindex(); //target tensor derives from columns
    usi ncolidx = rindex->ncolindex(); //rows are contraction indices
    usi nindex = nrowidx + ncolidx;

    //build the composite index range tuple
    tuple = new IndexRangeTuple(nindex);

    const usi* rmap = rperm->indexmap();
    const usi* lmap = lperm->indexmap();

    const usi* colptr = rindex->cols();
    const usi* rowptr = lindex->cols();

    IndexRange** rangeptr = tuple->begin();

    usi itot = 0;
    usi idx = 0;
    for (usi i=0; i < nrowidx; ++i, ++itot, ++rowptr)
    {
        idx = lmap[*rowptr];
        indexset[itot] = lindexset[idx];//lindexset->index(idx);
        tuple->set(itot, ltuple->get(idx));
    }

    for (usi i=0; i < ncolidx; ++i, ++itot, ++colptr)
    {
        idx = rmap[*colptr];
        indexset[itot] = rindexset[idx];//->index(idx);
        tuple->set(itot, rtuple->get(idx));
    }
}

Indexer*
MatrixConfiguration::rowindex(
    const IndexRangeTuplePtr& tuple,
    const PermutationPtr& p
) const
{
    return new Indexer(tuple, index_->rows(), index_->nrowindex(), p);
}

Indexer*
MatrixConfiguration::rowindex(
    const uli* sizes,
    const uli* offsets,
    const PermutationPtr& p
) const
{
    return new Indexer(sizes, offsets, index_->rows(), index_->nrowindex(), p);
}

MatrixIndex::MatrixIndex(
    usi nrows,
    usi ncols
) : rowindices_(0),
    colindices_(0),
    nrowindex_(nrows),
    ncolindex_(ncols),
    rowtype_(MatrixIndex::front),
    coltype_(MatrixIndex::back)
{
    init();
}

MatrixIndex::MatrixIndex(
    matrix_index_t row_type,
    usi nrows,
    matrix_index_t col_type,
    usi ncols
) : rowindices_(0),
    colindices_(0),
    nrowindex_(nrows),
    ncolindex_(ncols),
    rowtype_(row_type),
    coltype_(col_type)
{
    init();
}

MatrixIndex::~MatrixIndex()
{
    if (rowindices_)
        yeti_free_perm(rowindices_);
    if (colindices_)
        yeti_free_perm(colindices_);
}

void
MatrixIndex::init()
{
#if YETI_SANITY_CHECK
    if (rowtype_ == coltype_)
    {
        raise(SanityCheckError, "cannot have row and column types the same for matrix");
    }
#endif

    if (nrowindex_)
    {
        rowindices_ = yeti_malloc_perm();
    }
    if (ncolindex_)
    {
        colindices_ = yeti_malloc_perm();
    }

    usi* pmap = yeti_malloc_perm();
    if (rowtype_ == front)
    {
        for (usi i=0; i < nrowindex_; ++i)
        {
            rowindices_[i] = i; //rowindices are offset by the number of params
            pmap[i] = i; //parameters are only a single index and do not affect permutations
        }
        for (usi i=0; i < ncolindex_; ++i)
        {
            colindices_[i] = i + nrowindex_;
            pmap[i + nrowindex_] = i + nrowindex_;
        }
    }
    else
    {
        for (usi i=0; i < ncolindex_; ++i)
        {
            colindices_[i] = i;
            pmap[i + nrowindex_] = i;
        }
        for (usi i=0; i < nrowindex_; ++i)
        {
            rowindices_[i] = i + ncolindex_;
            pmap[i] = i + ncolindex_;
        }
    }

    short plus = 1;
    row_to_col_perm_ = new Permutation(pmap, ncolindex_ + nrowindex_, plus);
    //do not free! this assigns the pointer! yeti_free_perm(pmap);
}

Permutation*
MatrixIndex::get_row_column_permutation() const
{
    return row_to_col_perm_.get();
}

MatrixIndex::matrix_index_t
MatrixIndex::col_index_type() const
{
    return coltype_;
}

MatrixIndex::matrix_index_t
MatrixIndex::row_index_type() const
{
    return rowtype_;
}

usi
MatrixIndex::colindex(usi i) const
{
    return colindices_[i];
}

const usi*
MatrixIndex::cols() const
{
    return colindices_;
}

PermutationGroupPtr
MatrixIndex::expand_on_cols(const PermutationGroupPtr &pgrp) const
{
    make(newgrp, PermutationGroup, nindex(), pgrp, colindices_);
    return newgrp;
}

PermutationGroupPtr
MatrixIndex::expand_on_rows(const PermutationGroupPtr& pgrp) const
{
    make(newgrp, PermutationGroup, nindex(), pgrp, rowindices_);
    return newgrp;
}

usi
MatrixIndex::ncolindex() const
{
    return ncolindex_;
}

usi
MatrixIndex::nindex() const
{
    return ncolindex_ + nrowindex_;
}

usi
MatrixIndex::nrowindex() const
{
    return nrowindex_;
}

void
MatrixIndex::print(ostream& os) const
{
    os << Env::indent << "Matrix Index" << endl;

    ++Env::indent;

    os << Env::indent << "Row Indices: ";
    for (usi i=0; i < nrowindex_; ++i)
        os << " " << rowindices_[i];

    os << Env::indent << "Col Indices: ";
    for (usi i=0; i < ncolindex_; ++i)
        os << " " << colindices_[i];

    --Env::indent;
}

usi
MatrixIndex::rowindex(usi i) const
{
    return rowindices_[i];
}

const usi*
MatrixIndex::rows() const
{
    return rowindices_;
}


MatrixMap::MatrixMap(uli nrows, uli ncols)
    : nrows_(nrows), ncols_(ncols), blocks_(0)
{
    blocks_ = new CountableArray<Matrix>(nrows * ncols);
}

MatrixMap::~MatrixMap()
{
    delete blocks_;
}

Matrix*
MatrixMap::get(uli row, uli col) const
{
    uli idx = row * ncols_ + col;
    return blocks_->get(idx);
}

void
MatrixMap::insert(uli row, uli col, Matrix* m)
{
    uli idx = row * ncols_ + col;
    return blocks_->insert(idx, m);
}

void
MatrixMap::append(Matrix *m)
{
    insert(m->rowindex(), m->colindex(), m);
}

Matrix**
MatrixMap::begin() const
{
    return blocks_->begin();
}

Matrix**
MatrixMap::end() const
{
    return blocks_->end();
}
