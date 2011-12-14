#ifndef yeti_malloc_h
#define yeti_malloc_h

#include "class.h"

#include "malloc.hpp"
#include "thread.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

static const uli indexset_size = 6 * sizeof(unsigned long);

static const uli nindexset = 5000000;

/** @ingroup MemoryManagement
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

        uli size_;

        uli n_;

        uli refcount_;

        char* data_;

        std::string name_;

        uli offset_;

        bool static_mem_;

        uli
        search(
            uli start,
            uli stop
        );

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
            No thread lock is created.  This must only ever be called from a single thread.
        */
        FastMalloc(
            char* data,
            char* mallocd,
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

        void* malloc_no_lock();

        void* get_object(uli malloc_number);

        uli get_malloc_number(void* obj);

        uli nfree() const;

        /**
            Remove lock on the given block of memory, allowing
            other objects to use it from the memory pool
            @param ptr Pointer to the block to free
        */
        void free(const void* ptr);

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

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

