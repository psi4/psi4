#ifndef _yeti_runtime_h
#define _yeti_runtime_h

#include <iostream>
#include <map>

#include "class.h"
#include "index.hpp"
#include "thread.hpp"
#include "data.hpp"


namespace yeti {

class TensorTest;
class YetiRuntime {

    private:

        static std::map<std::string, IndexRangePtr> ranges_;

        //static std::map<IndexRange*, std::map<IndexRange*, usi> > subrange_depths_;

        static size_t* max_range_sizes_;

        static size_t* max_block_sizes_;
        
        static usi max_nindex_;

        static usi max_depth_;

        static uli nthread_;

        static uli nproc_;

        static uli memory_;

        static uli remaining_mem_;

        static uli cache_mem_;

        static uli me_;

        static Thread** threads_;

        static bool threaded_runtime_;

        static ThreadGroup* threadgrp_;

        static ThreadLock* printlock_;

        static ThreadLock* malloc_lock_;

        static void init_sizes();

        static void init_malloc();

    public:

        /**
            Get an index tuple specified by a given string, e.g.
            The string "i,j,k" would return a tuple of size 3
            where i,j,k specify specific index ranges.
            @param A string with index ranges comma separated
        */
        static IndexRangeTuple* get_index_tuple(const std::string& str);

        static IndexRange* get_range(const std::string& id);

        static void init(
            uli me,
            uli nproc,
            usi nindex,
            uli nthread,
            uli memory,
            uli cache
        );

        static void finalize();

        static void init_threads(uli nthread);

        static void reconfigure(uli me, uli nproc);

        static const size_t* max_block_sizes();

        static usi max_depth();

        static usi max_nindex();

        static const size_t* max_range_sizes();

        static uli me();

        static uli memory();

        static uli remaining_mem();

        static uli cache_mem();

        static void* allocate(size_t n);

        static void free(void* ptr, size_t n);

        static float matrix_multiply_cutoff;

        static uli nproc();

        static void register_index_range(
            IndexRange *range,
            const char* id1,
            const char* id2 = 0,
            const char* id3 = 0,
            const char* id4 = 0,
            const char* id5 = 0,
            const char* id6 = 0,
            const char* id7 = 0
        );

        static void register_index_range(
            const std::string& id,
            IndexRange* range
        );

        static void register_subrange(
            IndexRange* parent,
            IndexRange* child
        );

        static uli nthread();

        static void add_thread(Thread* thr);

        static ThreadGroup* get_thread_grp();

        static void lock_print();

        static void unlock_print();

        static void lock_malloc();

        static void unlock_malloc();

        static void set_threaded_runtime(bool flag);

        static bool is_threaded_runtime();

};

usi* usi_allocate(usi i, usi j);

uli* yeti_malloc_tile_location();
uli* yeti_malloc_indexset();
usi* yeti_malloc_perm();
void* yeti_malloc_indexptr();
void* yeti_malloc_matrix_generator();
void* yeti_malloc_matrix_generator_set();

void yeti_free_indexset(const uli* ptr);
void yeti_free_perm(usi* ptr);
void yeti_free_indexptr(void* ptr);
void yeti_free_tile_location(void* ptr);

void yeti_free_matrix_generator(void* ptr);
void yeti_free_matrix_generator_set(void* ptr);



}

#endif // Header Guard
