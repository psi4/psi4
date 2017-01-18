/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#ifndef PARALLELSCANNER_H_
#define PARALLELSCANNER_H_

#include "PsiInStream.h"
#include "PsiFileImpl.h"

namespace psi{

/** \brief The parallel (pun intended) of OutFile, except it reads files in.
 *
 *   As of right now the way to read data out of this file is via something
 *   like:
 *
 *   \code
 *   InFile MyFile("MyAwesomeFile.txt");
 *   while((bool)MyFile){
 *       std::string token;
 *       MyFile>>token;
 *       //Do stuff with token
 *   }
 *   \endcode
 *
 *   Note: Stuff comes out of MyFile as whatever type token is, so if, for
 *   example, you were reading integers declare token as an int.  Of course
 *   this means you need to know the format of your file, but heck, how else
 *   did you plan on reading it?
 *
 *   You want something else?  Code it (or buy Ryan a six-pack of Dortmunder
 *   Gold)
 *
 */
class InFile: public PsiInStream,private PsiFileImpl<std::ifstream>{
   private:
      ///Convient typedef of the FileBase
      typedef PsiFileImpl<std::ifstream> FileBase;

   public:
      /** \brief Constructor that opens a file with name "filename"
       *
       *  Upon construction passes NULL to base, so that the pointer is
       *  not set to cin.  This avoids free-ing cin in the Open call.
       *  By default no file is opened, but supplying a filename will open
       *  that file at the beginning mode.  Supplying a mode as well as a filename
       *  opens the file in that mode instead.  The mechanism of the
       *  constructor is to call Open so behavior is identical between
       *  initializing an object with a file and calling Open after calling
       *  the default constructor.
       *
       *  \param[in] filename Name of file to open (path relative to directory
       *                      Psi was called from)
       *
       *  \param[in] mode     The mode the file will be opened in
       *
       *
       */
      InFile(const std::string& filename="",const FileMode& mode=NOFILEMODE);

      ///Calls Close to release memory
      ~InFile(){Close();}

      /** \brief Opens "filename" in write mode "mode"
       *
       *  File contents are read from filename, and broadcast to
       *  all MPI processes and made available in Buffer_.  Stream_
       *  is then closed.
       */
      void Open(const std::string& filename,const FileMode& mode);

      ///Closes file and frees Stream_ pointer
      void Close(){FileBase::Close(Stream_);}
};

}



#endif /* PARALLELSCANNER_H_ */