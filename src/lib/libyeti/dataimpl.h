#ifndef data_impl_h
#define data_impl_h

#include "class.h"
#include "data.h"
#include "cache.h"
#include "sort.h"
#include "tile.h"
#include "index.h"
#include "sortimpl.h"
#include "malloc.h"
#include "tensor.h"
#include "exception.h"
#include "env.h"
#include "filler.h"
#include "runtime.h"

#include <algorithm>

#define ASSERT_DATA_NONNULL(x) if (!x) raise(SanityCheckError, "data is null");


namespace yeti {

template <typename data_t>
class DataTemplate : public Data {

    private:
        template <typename assign_t>
        void _assign(const assign_t* d);

        template <typename assign_t>
        void _assign(Data* d, const SortPtr& sort);

        template <typename assign_t>
        void _transfer(assign_t* data);

        template <typename assign_t>
        void _data(assign_t** d) const;

        template <typename acc_t>
        void _accumulate(const acc_t* d, acc_t scale);

        template <typename acc_t>
        void _accumulate(Data* d, const SortPtr& sort);

        template <typename cast_t>
        cast_t* _cast(cast_t* data);

        template <typename get_t>
        void _get(get_t** data) const;

    protected:
        data_t* data_;

    public:
        DataTemplate(size_t n);

        ~DataTemplate();

        void accumulate(const double* data, double scale = 1, const SortPtr& sort = 0);

        void accumulate(const int* data, int scale = 1, const SortPtr& sort = 0);

        void accumulate(const float* data, float scale = 1, const SortPtr& sort = 0);

        void accumulate(const quad* data, quad scale = 1, const SortPtr& sort = 0);

        void accumulate(Data* data, double scale = 1, const SortPtr& sort = 0);

        void assign(const double* data, const SortPtr& sort = 0);

        void assign(const int* data, const SortPtr& sort = 0);

        void assign(const float* data, const SortPtr& sort = 0);

        void assign(const quad* data, const SortPtr& sort = 0);

        char* buffer() const;

        double* cast(double *data);

        float* cast(float* data);

        int* cast(int* data);

        quad* cast(quad* data);

        void compute(Tile* tile, TileElementComputer* filler);

        bool equals(Tile* tile, TileElementComputer* filler);

        void debug_abort() const;

        void get(double** d) const;

        void get(int** d) const;

        void get(float** d) const;

        void get(quad** d) const;

        void assign(Data* data, const SortPtr& sort = 0);

        void transfer(double* data);

        void transfer(float* data);

        void transfer(int* data);

        void transfer(quad* data);

        void assign(data_t** d) const;

        void reference(void* data);

        void data(double **data) const;

        void data(int **data) const;

        void data(quad** data) const;

        void free();

        void malloc();

        void memset();

        float max_log() const;

        bool equals(const void* data);

        void print(std::ostream& os = std::cout) const;

        size_t n() const;

        size_t size() const;

        const void* pointer() const;

        TemplateInfo::type_t type();

        const char* type_name();

        void allocate();

        void clear();

        void sort(
            const SortPtr& sort,
            void* buffer
        );

        bool null() const;

        bool nonnull() const;

        void debug_fill(int mod, int denom);

        size_t element_size() const;

};


template <class data_t>
class MemoryBlockFactory : public DataBlockFactory {

    private:
        MemoryAllocationPtr mem_;

    public:
        MemoryBlockFactory();

        MemoryBlockFactory(uli storage);

        ~MemoryBlockFactory();

        Data::storage_t storage_type() const;

        DataBlock* get_block(Tile* tile);

        /**
            Computes the total size of all memory blocks in the tensor and
            allocated the pool of memory to be used for allocating blocks
            @param tensor
        */
        void configure(Tensor *tensor);

        void allocate(DataBlock* d);

        DataBlockFactory* copy() const;

};

template <class data_t>
class LocalDiskBlockFactory : public DataBlockFactory {

    private:
        DiskBufferPtr buffer_;

        uli offset_;

    public:
        LocalDiskBlockFactory();

        ~LocalDiskBlockFactory();

        Data::storage_t storage_type() const;

        DataBlock* get_block(Tile* tile);

        /**
            Computes the total size of all memory blocks in the tensor and
            allocated the pool of memory to be used for allocating blocks
            @param tensor
        */
        void configure(Tensor *tensor);

        void open_buffer(const std::string& name);

        void allocate(DataBlock* block);

        DataBlockFactory* copy() const;

};

template <class data_t>
class RecomputedBlockFactory :
    public DataBlockFactory
{

    private:
        ThreadedTileElementComputerPtr filler_;

    public:
        RecomputedBlockFactory(
            const ThreadedTileElementComputerPtr& filler
        );

        RecomputedBlockFactory(const TileElementComputerPtr& filler);

        RecomputedBlockFactory();

        Data::storage_t storage_type() const;

        DataBlock* get_block(Tile* tile);

        void configure(Tensor *tensor);

        DataBlockFactory* copy() const;

};

template <typename data_t>
void
DataBlock::allocate()
{
    if (is_retrieved())
        raise(SanityCheckError, "cannot allocate a block that is currently retrieved");

    if (data_)
    {
        std::cerr << "cannot reallocate data" << std::endl;
        abort();
    }
    else
    {
        data_ = new DataTemplate<data_t>(n_);
    }
}

template <typename data_t>
data_t*
Data::get() const
{
    data_t* data;
    get(&data);
    return data;
}

template <typename data_t>
DataTemplate<data_t>::DataTemplate(size_t n)
    : Data(n), data_(0)
{

}

template <typename data_t>
DataTemplate<data_t>::~DataTemplate()
{
    TypeInfo<data_t>::name();
}

template <typename data_t>
size_t
DataTemplate<data_t>::n() const
{
    return n_;
}

template <typename data_t>
size_t
DataTemplate<data_t>::size() const
{
    return n_ * sizeof(data_t);
}

template <typename data_t>
bool
DataTemplate<data_t>::null() const
{
    return !data_;
}

template <typename data_t>
bool
DataTemplate<data_t>::nonnull() const
{
    return data_;
}

template <typename data_t>
void
DataTemplate<data_t>::clear()
{
    data_ = 0;
}

template <typename data_t>
char*
DataTemplate<data_t>::buffer() const
{
    return reinterpret_cast<char*>(data_);
}

template <typename data_t>
void
DataTemplate<data_t>::debug_abort() const
{
    abort();
}

template <typename data_t>
template <typename assign_t>
void
DataTemplate<data_t>::_data(assign_t** d) const
{
    if (TypeCheck<data_t, assign_t>::mismatch)
        raise(SanityCheckError, "invalid data type assign");

    (*d) = reinterpret_cast<assign_t*>(data_);
}

template <typename data_t>
void
DataTemplate<data_t>::data(double **d) const
{
    _data<double>(d);
}

template <typename data_t>
void
DataTemplate<data_t>::data(int **d) const
{
    _data<int>(d);
}

template <typename data_t>
void
DataTemplate<data_t>::data(quad **d) const
{
    _data<quad>(d);
}

template <typename data_t>
const void*
DataTemplate<data_t>::pointer() const
{
    return data_;
}

template <typename data_t>
void
DataTemplate<data_t>::print(std::ostream& os) const
{
    if (data_ == 0)
        return;

    usi print_width = 5;
    data_t* dptr = data_;
    os << Env::indent;
    for (size_t i=0; i < n_; ++i, ++dptr)
    {
        if (i % print_width == 0 && i != 0)
            os << std::endl << Env::indent;

        os << " " << std::stream_printf(TypeInfo<data_t>::printf_str, *dptr);
    }
}

template <typename data_t>
bool
DataTemplate<data_t>::equals(const void *data)
{
    data_t* myptr = reinterpret_cast<data_t*>(data_);
    const data_t* dptr = reinterpret_cast<const data_t*>(data);
    for (size_t i=0; i < n_; ++i, ++dptr, ++myptr)
    {
        usi me = *myptr;
        usi d = *dptr;
        if ( d != me )
        {
            return false;
        }
    }
    return true;
}

template <typename data_t>
void
DataTemplate<data_t>::accumulate(const double *data, double scale, const SortPtr &sort)
{
    data_t dscale = scale;
    if (sort)
    {
        sort->accumulate(data, data_, dscale);
    }
    else
    {
        _accumulate<double>(data, dscale);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::accumulate(const quad *data, quad scale, const SortPtr &sort)
{
    data_t dscale = scale;
    if (sort)
    {
        sort->accumulate(data, data_, dscale);
    }
    else
    {
        _accumulate<quad>(data, dscale);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::accumulate(const int *data, int scale, const SortPtr& sort)
{
    if (sort)
        sort->template accumulate<int, data_t>(data, data_, scale);
    else
        _accumulate<int>(data, scale);
}

template <typename data_t>
void
DataTemplate<data_t>::accumulate(Data* data, double scale, const SortPtr &sort)
{
    data_t* dptr; data->get(&dptr);
    data_t dscale = scale;
    if (sort)
    {
        sort->accumulate(dptr, data_, dscale);
    }
    else
    {
        _accumulate<data_t>(dptr, dscale);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::accumulate(const float *data, float scale, const SortPtr &sort)
{
    data_t dscale = scale;
    if (sort)
    {
        sort->accumulate(data, data_, dscale);
    }
    else
    {
        _accumulate<float>(data, dscale);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::assign(const double *data, const SortPtr& sort)
{
    if (sort)
    {
        sort->sort(data, data_);
    }
    else
    {
        _assign<double>(data);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::assign(const quad *data, const SortPtr& sort)
{
    if (sort)
    {
        sort->sort(data, data_);
    }
    else
    {
        _assign<quad>(data);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::assign(const int *data, const SortPtr& sort)
{
    if (sort)
        sort->sort(data, data_);
    else
        _assign<int>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::assign(const float *data, const SortPtr& sort)
{
    if (sort)
        sort->sort(data, data_);
    else
        _assign<float>(data);
}

template <typename data_t>
template <typename cast_t>
cast_t*
DataTemplate<data_t>::_cast(cast_t* data)
{
    if (TypeCheck<cast_t,data_t>::match) //reinterpret cast to get rid of compiler warnings
        return reinterpret_cast<cast_t*>(data_);

    _transfer<cast_t>(data);
#if YETI_SANITY_CHECK
    if (!data)
    {
        std::cerr << "data is null" << std::endl;
        abort();
    }
#endif

    return data;
}

template <typename data_t>
float*
DataTemplate<data_t>::cast(float *data)
{
    return _cast<float>(data);
}

template <typename data_t>
double*
DataTemplate<data_t>::cast(double *data)
{
    return _cast<double>(data);
}

template <typename data_t>
quad*
DataTemplate<data_t>::cast(quad* data)
{
    return _cast<quad>(data);
}

template <typename data_t>
int*
DataTemplate<data_t>::cast(int *data)
{
    return _cast<int>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::compute(
    Tile* tile,
    TileElementComputer* filler
)
{
    filler->compute(tile, data_);
}

template <typename data_t>
bool
DataTemplate<data_t>::equals(
    Tile* tile,
    TileElementComputer* filler
)
{
    return filler->equals(tile, data_);
}

template <typename data_t>
void
DataTemplate<data_t>::free()
{
    if (data_ == 0)
    {
        std::cerr << "null data being freed" << std::endl;
        abort();
    }
    YetiRuntime::lock_malloc();
    delete[] data_;
    YetiRuntime::unlock_malloc();
    data_ = 0;
}

template <typename data_t>
void
DataTemplate<data_t>::malloc()
{
    if (data_)
    {
        std::cerr << "nonnull data being allocated" << std::endl;
        abort();
    }
    YetiRuntime::lock_malloc();
    data_ = new data_t[n_];
    YetiRuntime::unlock_malloc();
    ::memset(data_, 0, n_ * sizeof(data_t));
}

template <typename data_t>
void
DataTemplate<data_t>::debug_fill(int modulus, int denominator)
{
    double d = denominator;
    for (uli i=0; i < n_; ++i)
    {
        double n = i % modulus;
        data_[i] = n / d;
    }
}

template <typename data_t>
template <typename get_t>
void
DataTemplate<data_t>::_get(get_t** data) const
{
    if (TypeCheck<get_t,data_t>::mismatch)
        raise(SanityCheckError, "mismatched types in get");

    *data = reinterpret_cast<get_t*>(data_);
}

template <typename data_t>
void
DataTemplate<data_t>::get(int** data) const
{
    _get<int>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::get(float** data) const
{
    _get<float>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::get(double** data) const
{
    _get<double>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::get(quad** data) const
{
    _get<quad>(data);
}

template <typename data_t>
float
DataTemplate<data_t>::max_log() const
{
    //data_t* max = std::max_element(data_, data_ + n_);

    data_t* dptr = data_;
    data_t max = 0;
    for (uli i=0; i < n_; ++i, ++dptr)
    {
        if ( fabs(*dptr) > max )
            max = fabs(*dptr);
    }
    if ( max < 1e-16 ) //undefined
        return LOG_ZERO;
    else
        return log10(max);
}


template <typename data_t>
template <typename assign_t>
void
DataTemplate<data_t>::_transfer(assign_t* data)
{
    if (TypeCheck<data_t, assign_t>::match)
        raise(SanityCheckError, "equivalent types being transferred");

    ASSERT_DATA_NONNULL(data_);

    data_t* myptr = data_;
    assign_t* dptr = data;
    for (size_t i=0; i < n_; ++i, ++myptr, ++dptr)
    {
        (*dptr) = (*myptr);
    }
}

template <typename data_t>
void
DataTemplate<data_t>::transfer(double* data)
{
    _transfer<double>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::transfer(quad* data)
{
    _transfer<quad>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::transfer(float* data)
{
    _transfer<float>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::transfer(int* data)
{
    _transfer<int>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::assign(Data* data, const SortPtr& sort)
{
    switch(data->type())
    {
        case TemplateInfo::double_type:
            _assign<double>(data, sort);
            break;
        case TemplateInfo::float_type:
            _assign<float>(data, sort);
            break;
        case TemplateInfo::integer_type:
            _assign<int>(data, sort);
            break;
        case TemplateInfo::quad_type:
            _assign<quad>(data, sort);
            break;
    }
}

template <typename data_t>
void
DataTemplate<data_t>::assign(data_t** d) const
{
    (*d) = data_;
}

template <typename data_t>
void
DataTemplate<data_t>::reference(void* data)
{
    data_ = reinterpret_cast<data_t*>(data);
}

template <typename data_t>
void
DataTemplate<data_t>::memset()
{
    ::memset(data_, 0, n_ * sizeof(data_t));
}



template <typename data_t>
template <typename acc_t>
void
DataTemplate<data_t>::_accumulate(Data* d, const SortPtr& sort)
{
    acc_t* dptr; d->assign(&dptr);
    accumulate(dptr, sort);
}

template <typename data_t>
template <typename acc_t>
void
DataTemplate<data_t>::_accumulate(const acc_t* d, acc_t scale)
{
    ASSERT_DATA_NONNULL(data_);
    data_t* dataptr = data_;
    const acc_t* accptr = d;
    for (size_t i=0; i < n_; ++i, ++dataptr, ++accptr)
        (*dataptr) += (*accptr) * scale;
}

template <typename data_t>
template <typename assign_t>
void
DataTemplate<data_t>::_assign(const assign_t* d)
{
    ASSERT_DATA_NONNULL(data_);
    data_t* dataptr = data_;
    const assign_t* aptr = d;
    for (size_t i=0; i < n_; ++i, ++dataptr, ++aptr)
        (*dataptr) = (*aptr);
}

template <typename data_t>
template <typename assign_t>
void
DataTemplate<data_t>::_assign(Data* data, const SortPtr &sort)
{
    assign_t* dptr; data->get(&dptr);
    if (sort)
        sort->sort(dptr, data_);
    else
        _assign<assign_t>(dptr);
}

template <typename data_t>
void
DataTemplate<data_t>::sort(
    const SortPtr& sort,
    void* buffer
)
{
    data_t* tmp = reinterpret_cast<data_t*>(buffer);
    sort->sort<data_t>(data_, tmp);
    ::memcpy(data_, tmp, n_ * sizeof(data_t));
}

template <typename data_t>
size_t
DataTemplate<data_t>::element_size() const
{
    return sizeof(data_t);
}

template <typename data_t>
TemplateInfo::type_t
DataTemplate<data_t>::type()
{
    return TypeInfo<data_t>::type();
}

template <typename data_t>
const char*
DataTemplate<data_t>::type_name()
{
    return TypeInfo<data_t>::name();
}

template <typename data_t>
MemoryBlockFactory<data_t>::MemoryBlockFactory()
    : mem_(0)
{
}

template <typename data_t>
MemoryBlockFactory<data_t>::MemoryBlockFactory(uli mem)
    : mem_(new MemoryAllocation(mem))
{
}

template <typename data_t>
MemoryBlockFactory<data_t>::~MemoryBlockFactory()
{
}

template <typename data_t>
void
MemoryBlockFactory<data_t>::allocate(DataBlock* block)
{
    Data* d = block->data();
    void* ptr = mem_->get(d->size());
    d->reference(ptr);
    d->memset();
}

template <typename data_t>
DataBlockFactory*
MemoryBlockFactory<data_t>::copy() const
{
    return new MemoryBlockFactory<data_t>(mem_->size());
}

template <typename data_t>
DataBlock*
MemoryBlockFactory<data_t>::get_block(Tile* tile)
{
    DataBlock* block = new MemoryBlock(tile->config()->get_data_mode(), tile->ndata());
    block->allocate<data_t>();
    return block;
}

template <typename data_t>
Data::storage_t
MemoryBlockFactory<data_t>::storage_type() const
{
    return Data::in_core;
}

template <typename data_t>
void
MemoryBlockFactory<data_t>::configure(Tensor* tensor)
{
    if (mem_)
    {
        std::cerr << "allocated tensor twice" << std::endl;
        std::cerr << "make sure this is correct" << std::endl;
        return;
    }

    uli size = tensor->total_data_block_size();
    mem_ = new MemoryAllocation(size);
}

template <typename data_t>
RecomputedBlockFactory<data_t>::RecomputedBlockFactory()
    : filler_(0)
{
}

template <typename data_t>
RecomputedBlockFactory<data_t>::RecomputedBlockFactory(const ThreadedTileElementComputerPtr &filler)
    : filler_(filler)
{
}

template <typename data_t>
RecomputedBlockFactory<data_t>::RecomputedBlockFactory(const TileElementComputerPtr &filler)
    : filler_(new ThreadedTileElementComputer(filler))
{
}

template <typename data_t>
DataBlock*
RecomputedBlockFactory<data_t>::get_block(Tile* tile)
{
    DataBlock* block = new RecomputedBlock(
                            filler_.get(),
                            tile
                           );
    block->allocate<data_t>();
    return block;
}

template <typename data_t>
Data::storage_t
RecomputedBlockFactory<data_t>::storage_type() const
{
    return Data::recomputed;
}

template <typename data_t>
void
RecomputedBlockFactory<data_t>::configure(Tensor *tensor)
{
}

template <typename data_t>
DataBlockFactory*
RecomputedBlockFactory<data_t>::copy() const
{
    return new RecomputedBlockFactory<data_t>;
}

template <typename data_t>
LocalDiskBlockFactory<data_t>::LocalDiskBlockFactory(
)
    :
    buffer_(0),
    offset_(0)
{
}

template <typename data_t>
void
LocalDiskBlockFactory<data_t>::configure(Tensor *tensor)
{
    open_buffer(tensor->name());
}

template <typename data_t>
void
LocalDiskBlockFactory<data_t>::open_buffer(const std::string& name)
{
    if (buffer_)
    {
        std::cerr << "Disk block factory has already allocated buffer" << std::endl;
        abort();
    }

    std::string filename = name + ".tensor";
    buffer_ = new DiskBuffer(filename);
}

template <typename data_t>
LocalDiskBlockFactory<data_t>::~LocalDiskBlockFactory()
{
}

template <typename data_t>
void
LocalDiskBlockFactory<data_t>::allocate(DataBlock* block)
{
    block->set_buffer(offset_, buffer_);
    offset_ += block->data()->size();
}

template <typename data_t>
DataBlockFactory*
LocalDiskBlockFactory<data_t>::copy() const
{
    LocalDiskBlockFactory<data_t>* factory
        = new LocalDiskBlockFactory<data_t>;
    return factory;
}

template <typename data_t>
DataBlock*
LocalDiskBlockFactory<data_t>::get_block(Tile* tile)
{
    DataBlock* block
        = new LocalDiskBlock(
            tile->config()->get_data_mode(),
            tile->ndata(),
            tile->config()->cache
        );
    block->allocate<data_t>();
    return block;
}

template <typename data_t>
Data::storage_t
LocalDiskBlockFactory<data_t>::storage_type() const
{
    return Data::disk;
}

}

#endif
