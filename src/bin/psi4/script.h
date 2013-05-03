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

#ifndef _psi_src_lib_libscript_script_h_
#define _psi_src_lib_libscript_script_h_

#include <boost/shared_ptr.hpp>

namespace psi {
    
    /* 
     * Hopefully I can keep this abstract enough to handle multiple
     * scripting languages.
     */
    class Script {
    public:
        Script();
        virtual ~Script();
        
        /** Starts up the scripting language environment and adds
         * PSI4 specific functionality to the environment.
         */
        virtual void initialize() = 0;
        
        /** Shuts down the scripting language and frees any memory
         * that may have been allocated by the scripting language.
         */
        virtual void finalize() = 0;

        /** Run the input file script */
        virtual void run(FILE *input) = 0;
        
        static boost::shared_ptr<Script> language;
    };

    class Python : public Script {
    public:
        Python();
        virtual ~Python();
        
        virtual void initialize();
        virtual void finalize();
        virtual void run(FILE *input);
    };
}

#endif /* end of include guard: _psi_src_lib_libscript_script_h_ */
