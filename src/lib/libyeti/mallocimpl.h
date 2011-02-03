#ifndef yeti_MALLOCIMPL_H
#define yeti_MALLOCIMPL_H

#include "malloc.h"

#define USE_MALLOC_CLASS 1

namespace yeti {

/**
    @class FastMallocTemplate
    Class used for providing a memory pool for a specific class. This allocates
    blocks of the appropriate size for the given template parameter.
*/
template <class T>
class FastMallocTemplate {

    private:
        static uli nblocks_;

        static uli blocksize_;

        static std::string classname_;

    public:
        static FastMalloc* malloc;

        /**
            Constructor for registering a particular memory pool.
            Template parameter T correponds to the parent class.  All
            derived classes will have a special memory allocation.  The block
            sizes, though, will be determined by the constructor template
            parameter U.  For example, DataBlock is T and RecomputedBlock is U.
            This ensures that allocated block sizes for all derived classes of
            DataBlock will be large enough since the largest derived class is
            RecomputedBlock.
            @param u Pointer used for implicit template instantiaion
            @param classname
        */
        template <class U>
        FastMallocTemplate(U* u, const std::string& classname)
        {
            blocksize_ = sizeof(U);
            classname_ = classname;
            FastMalloc::queue(classname);
        }

        ~FastMallocTemplate();

        /**
            Allocate the memory pool

        */
        static void init(uli nblocks)
        {
            FastMalloc::unqueue(classname_);
            nblocks_ = nblocks;
            malloc = new FastMalloc(blocksize_, nblocks_, classname_);
        }

        static size_t blocksize()
        {
            return blocksize_;
        }

};

/**
    @class MallocOverride
    Data type used for placement new. This actually calls malloc in lieu of memory pool
    for particular objects.  This is needed in rare exceptions.
*/
class MallocOverride {
};

/**
    @class MallocTrack
*/
class MallocTrack {

    public:
        const char* file;
        int line;

        MallocTrack(
           const char* file,
           int line
        );
};

template <class T>
class Malloc {

    private:
        static bool initialized_;

    public:

#if USE_MALLOC_CLASS
        void operator delete(void* ptr)
        {
            if (initialized_)
                FastMallocTemplate<T>::malloc->free(ptr);
            else; //leak it
        }

        void* operator new(size_t size)
        {
#if YETI_SANITY_CHECK
            if (size > FastMallocTemplate<T>::blocksize())
            {
                std::cerr << "Blocksize requested is too large" << std::endl;
                abort();
            }
#endif
            if (initialized_)
                return FastMallocTemplate<T>::malloc->malloc();
            else
                return ::malloc(size);
        }

        void* operator new(size_t size, MallocOverride& o)
        {
            return malloc(size);
        }

        void* operator new(size_t size, MallocTrack& p)
        {
            void* ptr = FastMallocTemplate<T>::malloc->tracked_malloc(p.file, p.line);
            return ptr;
        }
#endif

        void register_new_object(
            const char* file,
            int line
        )
        {
            FastMallocTemplate<T>::malloc->register_malloc(this, file, line);
        }

        static void init_malloc(uli nblocks)
        {
#if USE_MALLOC_CLASS
            FastMallocTemplate<T>::init(nblocks);
            initialized_ = true;
#endif
        }

        static void uninit_malloc()
        {
            initialized_ = false;
        }

};

template <class T>
FastMallocTemplate<T>::~FastMallocTemplate()
{
    Malloc<T>::uninit_malloc();
    delete malloc;
}


#define DECLARE_SUBMALLOC(cls, mcls) \
typedef FastMallocTemplate<cls> Malloc##cls; \
template<> bool Malloc<cls>::initialized_ = false; \
template<> uli FastMallocTemplate<cls>::blocksize_ = 0; \
template<> uli FastMallocTemplate<cls>::nblocks_ = 0; \
template<> std::string FastMallocTemplate<cls>::classname_ = ""; \
template<> FastMalloc* FastMallocTemplate<cls>::malloc = 0; \
mcls* malloc_init_##cls = 0; \
static Malloc##cls malloc_decl_##cls (malloc_init_##cls, #cls)

#define DECLARE_MALLOC(cls) DECLARE_SUBMALLOC(cls,cls)

#if YETI_DEBUG_MALLOC
#define yeti_register_new(x) x->register_new_object(__FILE__,__LINE__)
#else
#define yeti_register_new(x)
#endif

}



#endif // MALLOCIMPL_H
