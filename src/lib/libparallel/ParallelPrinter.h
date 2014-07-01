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
#ifndef PARALLELPRINTER_H_
#define PARALLELPRINTER_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
namespace psi{


enum FileModes{NOMODE,READ,WRITE};

enum FileOptions{NOOPTIONS,END,APPEND,TRUNCATE,BINARY};

typedef std::fstream::openmode stdFMode;
typedef std::ostream& (*StreamManips)(std::ostream&);

/** \brief The base class for our new ParallelStream filesystem
 *
 * The baseclass is never intended to be created, hence the reason it has no
 * public members.  You should interface to this class through the appropriate
 * derived class:
 *
 * OutFile for writing to files
 * InFile  for reading from files
 * OutputFile for writing to the equivalent of output.dat
 *
 * TODO: Write PsiIO in terms of this unified file manager.  Likely you would
 * derive from OutFile and InFile: BinOutFile and BinInFile respectively.  The
 * interface to ParallelStream may need alterned slightly.
 */
class ParallelStream{
   private:
      ///Where data for all processes exists. Master process bcasts to here.
      std::stringstream Buffer;

      ///The file buffer (only master process reads/writes to here)
      std::fstream File;

      ///What Mode are we in
      FileModes Mode;

      inline stdFMode GetFileOptions(const FileOptions& Options);
      inline stdFMode ReadorWrite(const FileModes& Mode);
   protected:
      ///True if file exists and is open
      bool IsGood(){return File.good();}
      ///True if file is open (may be a new file)
      bool IsOpen(){return File.is_open();}
      ///True if the file has been read out
      bool IsEoF(){return Buffer.eof();}
      bool FileStatus(){return File;}
      void OpenImpl(const std::string& name,const FileModes InMode=NOMODE,
                const FileOptions=NOOPTIONS);

      ///Close Implementation.  Free memory. Set Output to NULL. Clears buffer
      void CloseImpl();

      ///Very basic set-up, defaults to writing files
      ParallelStream(FileModes InMode=NOMODE):Mode(InMode){}

      ///Closes files if open, frees memory
      virtual ~ParallelStream(){if(IsOpen())this->CloseImpl();}

      ///Reads data into the buffer (hence is for a write operation)
      template<class T>
      void Read(T& Input){if(Mode==WRITE)Buffer<<Input;}

      ///So you may pass std::endl or whatever other stream modifiers you want
      void Read(StreamManips fp){if(Mode==WRITE)Buffer<<fp;}

      ///Opposite of Read
      template<class T>
      void Write(T& Output){if(Mode==READ)Buffer>>Output;}

      ///Returns a copy of the buffer as a stringstream
      std::stringstream CopyBuffer(){std::stringstream ACopy(Buffer);return ACopy;}
};

class OutFile:private ParallelStream{
   public:

      ///Defaults to writing in append mode
      OutFile(const std::string& name="NONE",
            const FileOptions& options=APPEND);

      ///Defaults to writing in append mode
      void Open(const std::string&name,const FileOptions& options=APPEND);

      ///Calling this actually writes the file
      void Close(){CloseImpl();}

      ///Reads a stream into the buffer
      template<typename T>
      void operator<<(T& buffer){ParallelStream::Read(buffer);}

      ///So you may pass std::endl or whatever other stream modifiers you want
      void operator<<(StreamManips fp){ParallelStream::Read(fp);}

      ///For y'all who think streams suck...
      void Printf(const char* format,...);
};

class InFile:private ParallelStream{
   private:
      bool RealEoF;
   public:
      ///Sets this class up to read from a file
      InFile(const std::string& name="NONE",
            const FileOptions& options=NOOPTIONS);

      ///Reads the data into a buffer
      void Open(const std::string&name,
            const FileOptions& options=NOOPTIONS);

      ///Calling this closes the file
      void Close(){CloseImpl();RealEoF=false;}

      ///Dumps the file out word by word,returns true if there is more to read
      template<typename T>
      bool operator>>(T& buffer){
         if(IsEoF())RealEoF=true;
         ParallelStream::Write(buffer);
         return !(IsEoF()&&RealEoF);}

      ///True if the file exists
      bool Exists(){return FileStatus();}

      std::stringstream DumpBuffer(){return CopyBuffer();}

};

}//End namespace psi
#endif /* PARALLELPRINTER_H_ */
