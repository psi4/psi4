#ifndef yeti_malloc_h
#define yeti_malloc_h

#include "class.h"

#include "malloc.hpp"
#include "thread.hpp"

namespace yeti {

static const uli indexset_size = 6 * sizeof(unsigned long);

static const uli nindexset = 5000000;

/**
    @class FastMalloc
    Class that owns a particular memory pool, and allocates blocks quickly
    from an already malloc'd block.  All memory blocks are of uniform size
    to minimize memory fragmentation.
*/
class FastMalloc {

    private:
        static std::map<std::string, int>* registration_queue_;

        ThreadLock* lock_;

        char* mallocd_;

        char* dataptr_;

        char* mallocptr_;

        uli index_;

        uli size_;

        uli n_;

        uli refcount_;

        char* data_;

        std::string name_;

        std::map<const void*, std::string> tracker_;

        void search();

    public:
        /**
            @param size The size of the individual blocks
            @param n The number of blocks to allocate
            @param name A descriptive string for the type of blocks being allocated.
                        This usually corresponds to a class name.
        */
        FastMalloc(
            uli size,
            uli n,
            const std::string& name
        );

        /**
            Create the object but do not yet allocate the blocks
            @param name
        */
        FastMalloc(const std::string& name);

        /**
            Allocate blocks for the memory pool
            @param size
            @param n
        */
        void init(uli size, uli n);

        ~FastMalloc();

        /**
            Aborts if no free blocks are found
            @return A pointer to a free block in the memory pool
        */
        void* malloc();

        /**
            Remove lock on the given block of memory, allowing
            other objects to use it from the memory pool
            @param ptr Pointer to the block to free
        */
        void free(const void* ptr);

        void* tracked_malloc(const char* file, int line);

        void register_malloc(const void* ptr, const char* file, int line);

        /**
            Register a particular queue.  This is used in YetiRuntime
            sanity checks to ensure that a given memory pool is actually allocated.
            @param name
        */
        static void queue(const std::string& name);

        /**
            Unregister a particular queue.  This is used in YetiRuntime
            sanity checks to ensure that a given memory pool is actually allocated.
            @param name
        */
        static void unqueue(const std::string& name);

        /**
            Make sure all registered memory pools have been allocated.
            Aborts on failure.
        */
        static void check_queue();

};

/**
    @class MemoryAllocation
    Class encapsulating a memory pool which hands out blocks of data
*/
class MemoryAllocation : public smartptr::Countable {

    private:
        size_t size_;

        size_t remaining_;

        void* data_;

        char* ptr_;


    public:
        MemoryAllocation(size_t size);

        ~MemoryAllocation();

        void* get(size_t size);

        uli size() const;

};






}

#endif

