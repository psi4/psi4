#ifndef yeti_dataabstract_h
#define yeti_dataabstract_h

#include "class.h"

#include "data.hpp"
#include "cache.hpp"
#include "sort.hpp"
#include "tile.hpp"
#include "index.hpp"
#include "permutation.hpp"
#include "tensor.hpp"
#include "thread.hpp"
#include "filler.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"

#define MALLOC_CACHE_BLOCKS 0

namespace yeti {

class DataMode {

    public:
        typedef enum { read, write, accumulate } data_mode_t;

        data_mode_t flag;

        DataMode();

};

/**
    @class Data
    Abstract class for holding a block of data of a given data type.
*/
class Data :
    public Malloc<Data>
{

    protected:
        size_t n_;

        virtual void debug_abort() const = 0;

    public:
        typedef enum { in_core = 0, recomputed = 1, disk = 2 } storage_t;

        /**
            @param n The number of values
        */
        Data(uli n);

        virtual ~Data();

        /**
            @return The number of values
        */
        uli n() const;

        virtual char* buffer() const = 0;

        /**
            Accumulate a set of values.
            @param data The data to accumulate
            @param scale A scale factor to be applied
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void accumulate(const double* data, double scale = 1, const SortPtr& sort = 0);

        /**
            Accumulate a set of values.
            @param data The data to accumulate
            @param scale A scale factor to be applied
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void accumulate(const int* data, int scale = 1, const SortPtr& sort = 0);

        /**
            Accumulate a set of values.
            @param data The data to accumulate
            @param scale A scale factor to be applied
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void accumulate(const float* data, float scale = 1, const SortPtr& sort = 0);

        /**
            Accumulate a set of values.
            @param data The data to accumulate
            @param scale A scale factor to be applied
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void accumulate(const quad* data, quad scale = 1, const SortPtr& sort = 0);

        /**
            Accumulate a set of values.
            @param data The data to accumulate
            @param scale A scale factor to be applied
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void accumulate(Data* data, double scale = 1, const SortPtr& sort = 0) = 0;

        /**
            Assign a set of values to the data
            @param data The data to accumulate
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void assign(const double* data, const SortPtr& sort = 0);

        /**
            Assign a set of values to the data
            @param data The data to accumulate
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void assign(const int* data, const SortPtr& sort = 0);

        /**
            Assign a set of values to the data
            @param data The data to accumulate
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void assign(const float* data, const SortPtr& sort = 0);

        /**
            Assign a set of values to the data
            @param data The data to accumulate
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void assign(const quad* data, const SortPtr& sort = 0);

        /**
            Assign a set of values to the data
            @param data The data to accumulate
            @param sort An optional sort object for sorting the data prior to accumulation
            @throw If not reimplimented in subclass
        */
        virtual void assign(Data* data, const SortPtr& sort = 0) = 0;

        /**
            Transfer the values from my data into the given parameter buffer. This
            is used for changing the data type from, e.g. float to double.
            @param data The data buffer to insert values into
            @throw If the values are already of the given type
        */
        virtual void transfer(double* data) = 0;

        /**
            Transfer the values from my data into the given parameter buffer. This
            is used for changing the data type from, e.g. float to double.
            @param data The data buffer to insert values into
            @throw If the values are already of the given type
        */
        virtual void transfer(float* data) = 0;

        /**
            Transfer the values from my data into the given parameter buffer. This
            is used for changing the data type from, e.g. float to double.
            @param data The data buffer to insert values into
            @throw If the values are already of the given type
        */
        virtual void transfer(int* data) = 0;

        /**
            Transfer the values from my data into the given parameter buffer. This
            is used for changing the data type from, e.g. float to double.
            @param data The data buffer to insert values into
            @throw If the values are already of the given type
        */
        virtual void transfer(quad* data) = 0;

        /**
            Cast the values to a given data type.  If a type conversion is
            needed, a buffer must be provided as a parameter. If no type
            conversion is necessary, the parameter buffer is ignored and
            the underlying data is returned.
            @param data Buffer for converted values
            @return An array containing values with the appropriate data type
        */
        virtual double* cast(double* data) = 0;

        /**
            Cast the values to a given data type.  If a type conversion is
            needed, a buffer must be provided as a parameter. If no type
            conversion is necessary, the parameter buffer is ignored and
            the underlying data is returned.
            @param data Buffer for converted values
            @return An array containing values with the appropriate data type
        */
        virtual float* cast(float* data) = 0;

        /**
            Cast the values to a given data type.  If a type conversion is
            needed, a buffer must be provided as a parameter. If no type
            conversion is necessary, the parameter buffer is ignored and
            the underlying data is returned.
            @param data Buffer for converted values
            @return An array containing values with the appropriate data type
        */
        virtual long double* cast(quad* data) = 0;

        /**
            Cast the values to a given data type.  If a type conversion is
            needed, a buffer must be provided as a parameter. If no type
            conversion is necessary, the parameter buffer is ignored and
            the underlying data is returned.
            @param data Buffer for converted values
            @return An array containing values with the appropriate data type
        */
        virtual int* cast(int* data) = 0;

        /**
            Delete the memory in the block and return it to the memory pool
            it was obtained from.
        */
        virtual void free() = 0;

        /**
            Allocate a block of memory of the appropriate size
            for the block
        */
        virtual void malloc() = 0;

        /**
            Set the value pointed to by data to the beginning of the data array
            for this data block.
            @param data Pointer to beginning of array
            @throw If data type is wrong
        */
        virtual void get(int** data) const = 0;

        /**
            Set the value pointed to by data to the beginning of the data array
            for this data block.
            @param data Pointer to beginning of array
            @throw If data type is wrong
        */
        virtual void get(double** data) const = 0;

        /**
            Set the value pointed to by data to the beginning of the data array
            for this data block.
            @param data Pointer to beginning of array
            @throw If data type is wrong
        */
        virtual void get(float** data) const = 0;

        /**
            Set the value pointed to by data to the beginning of the data array
            for this data block.
            @param data Pointer to beginning of array
            @throw If data type is wrong
        */
        virtual void get(quad** data) const = 0;

        /**
            Point the underlying data to the block of memory pointed to by
            the parameter
            @param data
        */
        virtual void reference(void* data);

        /**
            @return Runtime value giving the underlying data type
        */
        virtual TemplateInfo::type_t type() = 0;

        /**
            @return Descriptive string giving type name
        */
        virtual const char* type_name() = 0;

        /**
            @return The total size of the memory block. This is equivalent to
            n * sizeof(datatype)
        */
        virtual size_t size() const = 0;

        /**
            @return The underlying data pointer
        */
        virtual const void* pointer() const = 0;

        virtual void print(std::ostream& os = std::cout) const = 0;

        /**
            @param data Block of data to compare
            @return Whether the block of data contains equivalent values
        */
        virtual bool equals(const void* data) = 0;

        /**
            Simultaneous compute and fill for cases in which data from
            filler cannot be made persistent past the function call
            @param tile The tile to compute the filler with
            @param filler The filler to compute values
        */
        virtual void compute(
            Tile* tile,
            TileElementComputer* filler
        ) = 0;

        virtual bool equals(
            Tile* tile,
            TileElementComputer* filler
        ) = 0;

        /**
            Set all values to zero
        */
        virtual void memset() = 0;

        /**
            Set the underlying data to a null pointer
        */
        virtual void clear() = 0;

        /**
            @return Whether the underlying data is null
        */
        virtual bool null() const = 0;

        /**
            @return Whether the underlying data is non-null
        */
        virtual bool nonnull() const = 0;

        /**
            Sort the data. This assumes the sort has been pre-configured
            for a given tile
            @param sort
            @param buffer A buffer to temporarily hold the sorted values
        */
        virtual void sort(const SortPtr& sort, void* buffer) = 0;

        /**
            @return The log of the maximum absolute value
        */
        virtual float max_log() const = 0;

        virtual void debug_fill(int modulus, int denominator) = 0;

        /**
            @return The size of an individual element in memory, e.g. 4 for 32-bit integer
        */
        virtual size_t element_size() const = 0;

        /**
            @return The data array
            @throw If incompatiable data type is requested
        */
        template <class data_t>
        data_t*
        get() const;

};

class Buffer :
    public YetiRuntimeCountable
{

    public:
        virtual void allocate_region(uli offset, uli size) = 0;

        virtual void read(uli offset, uli size, char* buffer) = 0;

        virtual void write(uli offset, uli size, const char* buffer) = 0;
};

/**
    @class DataBlock
    Abstract class for fetching data for Tile and Matrix classes
*/
class DataBlock :
    public Malloc<DataBlock>,
    public YetiRuntimeCountable
{

    protected:
        Data* data_;

        const DataMode* mode_;

        uli n_;

        /**
            These should be called after a local block of memory has either been
            allocated from malloc or the cache.  Init functions will expect memory
            to be available
        */
        virtual void _retrieve(uli threadnum);

        virtual void _release(uli threadnum);

        virtual void init_read(uli threadnum);

        virtual void init_write(uli threadnum);

        virtual void init_accumulate(uli threadnum);

        virtual void finalize_read(uli threadnum);

        virtual void finalize_write(uli threadnum);

        virtual void finalize_accumulate(uli threadnum);

    public:
        /**
            @param n The number of elements in the data block
        */
        DataBlock(
            const DataMode* mode,
            uli n
        );

        virtual ~DataBlock();

        /**
            @return The number of elements in the data block
        */
        uli n() const;

        /**
            @return The underlying data array
        */
        Data* data() const;

        virtual void print(std::ostream& os = std::cout) const = 0;

        /**
            Allocate data of a given type.
            @param Runtime value for the data type
        */
        void allocate(TemplateInfo::type_t type);

        /**
            Allocate data of a given type.
        */
        template <typename data_t>
        void allocate();

        /**
            This assumes sort has been pre-configured
            @param sort
            @param buffer Tempoary buffer for holding sorted values
        */
        virtual void sort(const SortPtr& sort, void* buffer);

        virtual void set_buffer(uli offset, const BufferPtr& buffer);

};

/**
    @param MemoryBlock
    Block of data that is persistent in memory
*/
class MemoryBlock : public DataBlock {

    protected:
        void _retrieve(uli threadnum);

        void _release(uli threadnum);

    public:
        MemoryBlock(
            const DataMode* mode,
            uli n
        );

        ~MemoryBlock();

        void sort(const SortPtr& sort, void* buffer);

        void print(std::ostream &os = std::cout) const;
};

/**
    @class CachedDataBlock
    Block of data that is allocated from a central cache system.
    Caching is used to minimize fetches from disk/remote locations.
*/
class CachedDataBlock : public DataBlock {

    protected:
        friend class DataCacheEntry;

        LayeredDataCachePtr main_cache_;

        /**
            A data cache holding blocks of the appropriate size
            for this data block
        */
        DataCache* cache_;

        /**
            A cache entry containing a block of data that
            this cache entry is linked to
        */
        DataCacheEntry* cache_entry_;

        void _retrieve(uli threadnum);

        void _release(uli threadnum);

        void finalize();

    public:
        CachedDataBlock(
            const DataMode* mode,
            uli n,
            const LayeredDataCachePtr& cache
        );

        /**
        */
        virtual ~CachedDataBlock();

        void clear();

        virtual void print(std::ostream &os = std::cout) const;

};

/**
    @class RecomputedBlock
    Cached data block which fetches values by performing some sort of recompute
*/
class RecomputedBlock : public CachedDataBlock {

    private:
        ThreadedTileElementComputerPtr filler_;

        Tile* tile_;

    protected:
        void init_read(uli threadnum);

        void init_write(uli threadnum);

        void init_accumulate(uli threadnum);

    public:

        RecomputedBlock(
            const ThreadedTileElementComputerPtr& filler,
            Tile* tile
        );

        ~RecomputedBlock();

        /**
            Debug constructor
        */
        RecomputedBlock(
            uli ndata,
            const ThreadedTileElementComputerPtr& filler,
            const LayeredDataCachePtr& cache,
            Tile* tile
        );

        void print(std::ostream &os = std::cout) const;

};

/**
    @class SortedBlock
    A block which upon leaving the cache accumulates its values
    in a sorted fasion to a parent data block.  Upon reading, it
    sorts values as fetched.
    The caching is meant to minimize the number of sorts.  This is mostly
    in used in tensor contractions requiring the values to be
    symmetrized.
*/
class SortedBlock : public CachedDataBlock {

    private:
        DataBlockPtr parent_;

        SortPtr sort_;

        IndexRangeTuplePtr tuple_;

        /**
            This pointer is not owned by the list
        */
        Permutation** perms_;

        uli nperms_;

    protected:
        void init_read(uli threadnum);

        void init_accumulate(uli threadnum);

        void finalize_accumulate(uli threadnum);

    public:
        /**
            @param parent The data block whose values are to be sorted
            @param cache
            @param tuple  The index range tuple used to configure the sort
            @param perms  The list of permutations defining the accumulations
            @param nperms   The number of permutations
        */
        SortedBlock(
            const DataMode* mode,
            const DataBlockPtr& parent,
            const LayeredDataCachePtr& cache,
            const IndexRangeTuplePtr& tuple,
            Permutation** perms,
            uli nperms
        );

        ~SortedBlock();

        uli nperms() const;

        Permutation** perms() const;

        void print(std::ostream &os = std::cout) const;

};

class DiskBuffer :
    public Buffer
{

    private:
        std::string filename_;

        int fileno_;

        uli size_;

    protected:
        void _retrieve(uli threadnum);

        void _release(uli threadnum);

    public:
        DiskBuffer(const std::string& filename);

        ~DiskBuffer();

        uli size() const;

        void allocate_region(uli offset, uli size);

        void read(uli offset, uli size, char* buffer);

        void write(uli offset, uli size, const char* buffer);

        const std::string& filename() const;

};

class LocalDiskBlock :
    public CachedDataBlock
{
    private:
        BufferPtr buffer_;

        uli offset_;

    protected:

        void init_read(uli threadnum);

        void init_write(uli threadnum);

        void init_accumulate(uli threadnum);

        void finalize_read(uli threadnum);

        void finalize_write(uli threadnum);

        void finalize_accumulate(uli threadnum);

    public:
        LocalDiskBlock(
            const DataMode* mode,
            uli n,
            const LayeredDataCachePtr& cache
        );

        ~LocalDiskBlock();

        void set_buffer(uli offset, const BufferPtr& buffer);
};


/**
    @class DataBlockFactory
*/
class DataBlockFactory : public smartptr::Serializable {

    public:

        virtual ~DataBlockFactory();

        /** Allocate the space for the acutal data to be held by the data block.
            By default, no space is allocated.  Generally, only data blocks held
            permantently in memory will have space allocated by this method.
            @param d The data block to allocate memory for
        */
        virtual void allocate(DataBlock* d);

        /**
            Depending on the tensor, certain aspects of the factory may need to be configured.
            @param tensor
        */
        virtual void configure(Tensor* tensor) = 0;

        virtual DataBlock* get_block(Tile* tile) = 0;

        virtual Data::storage_t storage_type() const = 0;

        virtual DataBlockFactory* copy() const = 0;


};


}

#endif
