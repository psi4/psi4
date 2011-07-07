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


  $Id: nodefaults.h 1602 2009-12-27 19:53:06Z rjharrison $
*/


#ifndef MADNESS_WORLD_NODEFAULTS_H__INCLUDED
#define MADNESS_WORLD_NODEFAULTS_H__INCLUDED

/// \file nodefaults.h
/// \brief Implements NO_DEFAULTS


/// Disables default copy constructor and assignment operators

/// From http://home.pb.net/~tglenn/CPPNOTES.html. Inherit from this
/// class in order to inhibit the automatic generation of the default
/// copy constructor and the default assignment operator of the
/// derived class.
class NO_DEFAULTS {
public:
    NO_DEFAULTS() {};
private:
    // hide these - DO NOT IMPLEMENT!
    NO_DEFAULTS(const NO_DEFAULTS&);
    NO_DEFAULTS& operator=(const NO_DEFAULTS&);
};

#endif // MADNESS_WORLD_NODEFAULTS_H__INCLUDED

