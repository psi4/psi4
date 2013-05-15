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

#ifndef DATAPOLICIES_H
#define DATAPOLICIES_H

#include "tensor.hpp"

namespace yeti {

class DoNothingDataRetrieve {
    public:
        void retrieve(
            DataNode* node,
            const MetaDataNode* parent,
            const uli* indices
        );

        bool need_recompute() const;
};

class RecomputeDataRetrieve {
    private:
        template <typename data_t>
        void compute(
            TensorElementComputer* filler,
            DataNode* node,
            const uli* indices
        );
        
    public:
        void retrieve(
            DataNode* node,
            const MetaDataNode* parent,
            const uli* indices
        );

        bool need_recompute() const;
};

class DoNothingDataRelease {
    public:
        void release(
            DataNode* node,
            const MetaDataNode* parent
        );
};

class DoNothingDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class DeleteAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class CommitAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class ClearAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class DeleteAllDataClear {
    public:
        void clear(
            TensorBlock* block
        );
};

class ClearAllData {
    public:
        void clear(
            TensorBlock* block
        );
};

}

#endif // DATAPOLICIES_H
