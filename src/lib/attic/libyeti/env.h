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

#ifndef _tiled_tensor_Env_h
#define _tiled_tensor_Env_h

#include <iostream>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class IndentTracker {

    private:
        int indentsize_;

        int n_;

    public:
        IndentTracker();

        IndentTracker
        operator++();

        IndentTracker
        operator--();

        void
        add_indent(std::ostream& os) const;

};



/**
 * The Env class determines how I/O should be performed
 * based on the current node number and the desired
 * locality of the print statement
 */
class Env {

 protected:
   /// Whether the Env object has been initialized
   static int initialized_;
   /// The current node number
   static int me_;
   /// The stream object for output from node 0 only
   static std::ostream *out0_;
   /// The stream object for output from all nodes
   static std::ostream *outn_;

 public:
   static IndentTracker& indent;

   static void init(int me, const char *outFileName = "");
   static void init(int me, std::ostream *stream);
   static void free();
   /// Return nonzero if Env has been initialized.
   static int initialized() { return initialized_; }
   /// Return an ostream that writes from all nodes.
   static std::ostream &outn() { return *outn_; }
   /// Return an ostream for error messages that writes from all nodes.
   static std::ostream &errn() { return *outn_; }
   /// Return an ostream that writes from node 0.
   static std::ostream &out0() { return *out0_; }
   /// Return an ostream for error messages that writes from node 0.
   static std::ostream &err0() { return *out0_; }
   
};

std::ostream&
operator<<(std::ostream& os, const IndentTracker& indent);

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // Header Guard
