/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef yeti_tensorbranch_h
#define yeti_tensorbranch_h

#include <libsmartptr/serialize.h>
#include "mapimpl.h"

#include "data.hpp"
#include "index.hpp"
#include "tensor.hpp"
#include "tensorbranch.hpp"
#include "tensorblock.hpp"
#include "node.hpp"
#include "permutation.hpp"
#include "tensorparser.hpp"
#include "filler.hpp"
#include "elementop.hpp"
#include "matrix.hpp"
#include "contraction.hpp"

#include "cache.h"

#include "gigmatrix.h"

#include <list>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif



namespace yeti {

class TensorBranch :
    public smartptr::Countable,
    public MempoolVirtualAddressMalloc
{

    private:
        TemplateInfo::type_t element_type_;

        size_t element_size_;

        TensorBlock* parent_block_;

        uli parent_malloc_number_;
        
        TensorElementComputer* filler_;

        TensorBranch* last_branch_location_;

        MetaDataNode* mdnode_;

        MemoryPool* mempool_;

        uli malloc_number_;

        TensorDataController* first_controller_;

        TensorDataController* last_controller_;

        void allocate_data_controller(size_t blocksize);

        TensorIndexDescr* descr_;

    public:
        TensorBranch(
            const uli* indexset,
            usi depth,
            TensorBlock* parent
        );

        ~TensorBranch();

        void set_parent(TensorBlock* parent);
        
        void allocate_data_controller(StorageBlock* block);

        void set_element_computer(TensorElementComputer* filler);

        TensorController* get_tensor_controller() const;

        TensorDataController* first_data_controller() const;

        TensorIndexDescr* get_descr() const;

        void set_descr(TensorIndexDescr* descr);
        
        TensorElementComputer* get_element_computer() const;

        TensorBranch* get_last_branch_location() const;

        TensorBlock* get_parent_block() const;

        MemoryPool* get_metadata_mempool() const;

        MetaDataNode* get_node() const;

        uli get_parent_malloc_number() const;

        void set_malloc_number(uli num);

        uli get_malloc_number() const;

        void load_metadata(long offset);

        void allocate(DataNode* node);

        void set_element_type(TemplateInfo::type_t type);

        uli ncontrollers_nonzero() const;

        size_t size_data_wasted();

        void sort_metadata_into(
            Sort* sort,
            TensorBranch* new_branch
        );

        void sort_data_into(
            Sort* sort,
            TensorBranch* new_branch
        );

        TemplateInfo::type_t element_type() const;

        /**
            The size of elements in the given block... so
            potentially double, float, etc whatever
        */
        size_t element_size() const;

        void realign_memory_pool(
            TensorBranch* old_branch,
            TensorBranch* new_branch
        );

        void set_metadata_mempool(MemoryPool* mempool);

};

class TensorDataController :
    public smartptr::Countable,
    public MempoolVirtualAddressMalloc
{

    private:
        TensorBranch* branch_;

        DataNode* first_data_node_;

        DataNode* last_data_node_;

        size_t remaining_;

        size_t total_size_;

        char* data_;


    public:
        TensorDataController(
            TensorBranch* branch,
            StorageBlock* storage_block
        );

        ~TensorDataController();

        char* get_data() const;

        size_t get_total_size() const;

        size_t get_data_size() const;

        size_t get_remaining() const;

        DataNode* first() const;

        void realign_memory_pool(
            TensorBranch* old_branch,
            TensorBranch* new_branch
        );

        bool register_node(DataNode* node, size_t blocksize);

        void set_data(char* data);

        TensorDataController* next;

};

}

#endif
