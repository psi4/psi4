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


  $Id: archive_type_names.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#define MAD_ARCHIVE_TYPE_NAMES_CC
#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE
#include <world/world.h>
#include <cstring>

namespace madness {
    namespace archive {

        template <typename T> const unsigned char archive_typeinfo<T>::cookie;


        // Forces initialization of type names at startup
        // (breaks on shared libs ?)
        static BaseArchive fred_and_mary_sitting_under_a_tree;

        void archive_initialize_type_names() {
            static  bool initialized = false;
            if (initialized) return;
            initialized = true;


            for (int i=0; i<255; ++i) archive_type_names[i] = "invalid";
            archive_type_names[255] = "unknown/user-defined";

            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned char);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned short);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned int);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned long);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned long long);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed char);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed short);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed int);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed long);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed long long);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(bool);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(float);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(double);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(long double);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::complex<float>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::complex<double>);

            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<char>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned char>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<short>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned short>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<int>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned int>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<long>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned long>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<bool>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<float>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<double>);

            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::string);

            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<int>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<long>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<float>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<double>);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor< std::complex<float> >);
            ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor< std::complex<double> >);

        }
    }
}
