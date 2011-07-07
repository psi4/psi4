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
#ifndef OUTPUTWRITER_H_
#define OUTPUTWRITER_H_

#include <fstream>

class OutputWriter
{
public:
  //***************************************************************************
  static OutputWriter* instance()
  {
    if (!_instance)
      _instance = new OutputWriter();
    return _instance;
  }
  //***************************************************************************

  //***************************************************************************
  void init_debug(char* filename)
  {
    dostr = new std::ofstream(filename);
  }
  //***************************************************************************

  //***************************************************************************
  void init_log(char* filename)
  {
    lostr = new std::ofstream(filename);
  }
  //***************************************************************************

  //***************************************************************************
  void init_eigv(char* filename)
  {
    evostr = new std::ofstream(filename);
  }
  //***************************************************************************

  //***************************************************************************
  virtual ~OutputWriter()
  {
    if (dostr) delete dostr;
    if (lostr) delete lostr;
    if (evostr) delete evostr;
  }
  //***************************************************************************

  //***************************************************************************
  std::ofstream* debug_stream()
  {
    return dostr;
  }
  //***************************************************************************

  //***************************************************************************
  std::ofstream* log_stream()
  {
    return lostr;
  }
  //***************************************************************************

  //***************************************************************************
  std::ofstream* eigv_stream()
  {
    return evostr;
  }
  //***************************************************************************

private:
  //***************************************************************************
  OutputWriter()
  {
    dostr = 0;
    lostr = 0;
    evostr = 0;
  }
  //***************************************************************************

  //***************************************************************************
  static OutputWriter* _instance;
  //***************************************************************************

  //***************************************************************************
  std::ofstream* dostr;
  //***************************************************************************

  //***************************************************************************
  std::ofstream* lostr;
  //***************************************************************************

  //***************************************************************************
  std::ofstream* evostr;
  //***************************************************************************
};

OutputWriter* OutputWriter::_instance = 0;

#endif /* OUTPUTWRITER_H_ */
