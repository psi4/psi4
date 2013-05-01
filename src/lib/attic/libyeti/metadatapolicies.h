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

#ifndef METADATAPOLICIES_H
#define METADATAPOLICIES_H

#include "tensor.hpp"

namespace yeti {

class DoNothingMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class RealignMemoryPoolMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class RecomputeMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class SortMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class DoNothingMetaDataRelease {
    public:
        void release(
            MetaDataNode* mdnode
        );
};

}

#endif // METADATAPOLICIES_H
