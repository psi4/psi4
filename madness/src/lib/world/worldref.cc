/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

#include <world/worldref.h>
#include <world/worldmutex.h>
#include <iostream>

namespace madness {
    namespace detail {
        RemoteCounter::pimpl_mapT RemoteCounter::pimpl_map_;

        /// Clean-up the implementation object

        /// Here we check that the pimpl has been initialized, and if so, we
        /// release the current reference. If the count drops to zero, then
        /// this is the last reference to the pimpl and it should be deleted.
        void RemoteCounter::destroy() {
            if(pimpl_.is_local()) {
                if(pimpl_->release()) {
                    // No one else is referencing this pointer.
                    // We can safely dispose of it.

#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                    print(">>> RemoteCounter::unregister_ptr_: key=", pimpl_->key(), ", value=", pimpl_);
#endif
                    unregister_ptr_(pimpl_->key());
                    delete pimpl_.get();
                }
            }

            pimpl_ = WorldPtr<implT>();
        }

        RemoteCounter& RemoteCounter::operator=(const RemoteCounter& other) {
            WorldPtr<implT> temp = other.pimpl_;

            if(pimpl_ != temp) {
                if(temp)
                    temp->add_ref();
                destroy();
                pimpl_ = temp;
            }

            return *this;
        }

        std::ostream& operator<<(std::ostream& out, const RemoteCounter& counter) {
            out << "RemoteCounter( owner=" << counter.owner() << " worldid=" <<
                    counter.get_worldid() << " use_count=" << counter.use_count() << ")";
            return out;
        }
    } // namespace detail
} // namespace madness

