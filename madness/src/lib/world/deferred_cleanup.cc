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


  $Id: $
*/

#include <world/deferred_cleanup.h>
#include <world/world.h>

namespace madness {
    namespace detail {

        void DeferredCleanup::destroy(bool mode) {
            mutex_.lock();
            destroy_ = true;
            mutex_.unlock();
        }

        bool DeferredCleanup::destroy() const {
            mutex_.lock();
            bool mode = destroy_;
            mutex_.unlock();
            return mode;
        }

        void DeferredCleanup::add(const void_ptr& obj) {
            ScopedMutex<RecursiveMutex> lock(mutex_);
            // if we do not store obj, it will be destroyed when the calling
            // function (DeferredDeleter::operator()) exits.
            if(! destroy_)
                deferred_.push_back(obj);
        }

        void DeferredCleanup::do_cleanup() {
            std::list<void_ptr> cleaner;

            do {
                cleaner.clear(); // Cleanup memory
                mutex_.lock();
                cleaner.swap(deferred_);
                mutex_.unlock();
            } while(! cleaner.empty());
        }


        std::shared_ptr<DeferredCleanup> DeferredCleanup::get_deferred_cleanup(const World& w) {
            return w.gop.deferred;
        }

    }  // namespace detail
}  // namespace madness
