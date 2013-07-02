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

#include "env.h"
#include <string.h>
#include <iostream>
#include <fstream>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define ENV_INDENT_SIZE 4

static IndentTracker _env_indent_tracker;

int Env::initialized_ = 0;
int Env::me_= -1;
ostream *Env::out0_;
ostream *Env::outn_;
IndentTracker& Env::indent = _env_indent_tracker;

// TODO figure out a better method that will free memory automagically
std::ostream&
yeti::operator<<(std::ostream& stream, const IndentTracker& indent)
{
    indent.add_indent(stream);
    return stream;
}

IndentTracker::IndentTracker()
    : n_(0), indentsize_(ENV_INDENT_SIZE)
{
}


void
IndentTracker::add_indent(std::ostream& os) const
{
    if (!n_)
        return;

    for (int i=0; i < n_ * indentsize_; ++i)
        os << " ";
}


/**
 * Set up the output stream using a filename, which should be
 * omitted if std::cout is to be used.
 */
void
Env::init(int me, const char *outFileName){
    if(initialized_) return;
    me_ = me;
    if(strcmp("", outFileName)){
        outn_ = new ofstream(outFileName, ios_base::app | ios_base::out );
    }else{
        outn_ = &cout;
    }
    if(me){
        out0_ = new ofstream("/dev/null");
    }else{
        out0_ = outn_;
    }
    initialized_ = 1;
}

/**
 * Set up the output stream using an exising stream.
 */
void
Env::init(int me, ostream *stream){
    if(initialized_) return;
    me_ = me;
    outn_ = stream;
    if(me){
        out0_ = new ofstream("/dev/null");
    }else{
        out0_ = outn_;
    }
    initialized_ = 1;
}

IndentTracker
IndentTracker::operator++()
{
    ++n_;
    return (*this);
}

IndentTracker
IndentTracker::operator--()
{
    --n_;
    return (*this);
}

/**
 * Free the memory associated with Env
 */
void
Env::free()
{
    if(initialized_){
        initialized_ = false;
    }
}

