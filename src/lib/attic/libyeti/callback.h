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

#ifndef yeti_CALLBACK_H
#define yeti_CALLBACK_H

namespace yeti {

class Callback {

    public:
        virtual void callback() = 0;
};


template <class T>
class tmpl_Callback :
    public Callback
{
    private:
        typedef void (T::*void_fxn_ptr)(void);

        T* obj_;

        void_fxn_ptr fxn_;

    public:
        tmpl_Callback(
            T* obj,
            void_fxn_ptr fxn
        ) : obj_(obj),
            fxn_(fxn)
        {
        }

        void callback()
        {
            ((*obj_).*fxn_)();
        }

        void* operator new(size_t size, char* ptr)
        {
            return ptr;
        }

        void operator delete(void* ptr, size_t size)
        {
            //do nothing
        }
};

}

#endif // CALLBACK_H
