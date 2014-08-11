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
#ifndef STREAMBASE_H_
#define STREAMBASE_H_
#include<iostream>
#include<boost/shared_ptr.hpp>

//This is the signature of std::endl and other sorts of iomanip things
typedef std::ostream& (*StreamManips)(std::ostream&);

///A shared output stream
typedef boost::shared_ptr<std::ostream> SharedStream;

namespace psi{


/** \brief The base class for Psi4's new parallel safe printing system
 *
 *  At this level we take care of the "who gets to print" sort of details.
 *  The call PsiStreamBase::ImSpecial is the function that determines this; if
 *  down the road people want to change this, this is where to do it at.  Other
 *  than that, this class takes care of copying the Buffer_ and synching it
 *  across MPI processes.
 */
class PsiStreamBase{
   private:
      ///Makes this a copy of other, copies data in Buffer_, not &Buffer
      void Clone(const PsiStreamBase& other);

   protected:
      ///This is where each MPI task ultimately writes from or to
      std::stringstream Buffer_;

      ///Returns true if this is the lucky MPI process that gets to read/write
      bool ImSpecial()const;

   public:

      ///Calls Clone for actual copy
      PsiStreamBase(const PsiStreamBase& other){this->Clone(other);}

      ///Calls Clone for assignment iff this!=&other, returns *this
      const PsiStreamBase& operator=(const PsiStreamBase& other);

      ///Default constructor of stringstream is called
      PsiStreamBase(){}

      ///Memory not worried about until we get down the class tree to files
      virtual ~PsiStreamBase(){}
};


class PsiOutStream:public PsiStreamBase{
   private:
      ///Dumps Buffer to Stream_
      inline void DumpBuffer();

      /** \brief This is the interface for Printf and << operator to
       *           Buffer_ for most writable objects
       *
       *  \param[in] Input An object that can be passed to an ostream
       */
      template<typename T>
      void Write2Buffer(const T& Input);

      /** \brief This is the interface for Printf and << operator to
       *     Buffer_ for special iostream functions.
       *
       *     \param[in] fp iostream functions like std::endl
       */
      void Write2Buffer(StreamManips fp);

      /** \brief Copies Stream_ of other over to this
       *
       *   Note that at this level Stream_ can only be cout or cerr
       *   thus no memory leak can occur and we allow the copy
       */
      void Clone(const PsiOutStream& other);
   protected:

      ///This is where the privileged MPI task gets to write Buffer_ to
      SharedStream Stream_;

   public:

      ///Makes an OutStreamBase that defaults to std::cout
      PsiOutStream(SharedStream Stream=SharedStream());

      ///Creates this by copying other, via base copy constructor
      PsiOutStream(const PsiOutStream& other);

      ///Allows assignment of other to this, calls base assignment
      const PsiOutStream& operator=(const PsiOutStream& other);

      ///No memory to free up (we don't delete cout or cerr)
      virtual ~PsiOutStream(){}

      ///Flushes the Stream
      inline void Flush();

      /** \brief Allows printf (or fprintf) like syntax to this object
       *
       *   This function takes your format stream and converts it to a 1,000
       *   character MAXIMUM char* array.  That array is then passed to Buffer_.
       *   If the current MPI process is rank 0 on MPI_COMM_WORLD, Buffer_ is
       *   then written to Stream_.  No locks should be used around your Printf
       *   statements, nor is flushing needed.  For reference the accepted
       *   format flags are (stolen from cplusplus.com):
       *
       *   %[flags][width][.precision][length]specifier
       *
       *   specifier    Output
       *     d or i  Signed decimal integer
       *      u      Unsigned decimal integer
       *      o   Unsigned octal integer
       *      x   Unsigned hexadecimal integer
       *      X   Unsigned hexadecimal integer (uppercase)
       *      f   Decimal floating point, lowercase
       *      F   Decimal floating point, uppercase
       *      e   Scientific notation (mantissa/exponent), lowercase
       *      E   Scientific notation (mantissa/exponent), uppercase
       *      g   Use the shortest representation: %e or %f
       *      G   Use the shortest representation: %E or %F
       *      a   Hexadecimal floating point, lowercase
       *      A   Hexadecimal floating point, uppercase
       *      c   Character
       *      s   String of characters
       *      p   Pointer address
       *      n   Nothing printed
       *
       *
       *   \param[in] format The message plus format codes you want to pass
       *   \param[in] ...    The arguments for each format code given in "format"
       *
       */
      void Printf(const char* format,...);

      ///Allows c++ like interface for most objects
      template<typename T>
      void operator<<(const T& Input){Write2Buffer(Input);}

      ///Allows c++ like interface for things like std::endl
      void operator<<(StreamManips fp){Write2Buffer(fp);}

      /** \brief Makes a banner
       *
       *   \param[in] message   The message that will go in the banner
       *   \param[in] delimiter The character that will comprise the banner
       *   \param[in] width     The number of characters wide the banner will be
       *
       *
       *   Let's say you want to print "Hello World!" in a banner comprised of
       *   '*' characters, and for simplicity (your banner should be 80 chars
       *   long) you only wanted a 16 character long banner.  If you passed
       *   message="Hello World!",delimiter='*', and width=16, you would get
       *   (hopefully you are viewing this in a monospaced font):
       *
       *   0123456789ABCDEF
       *
       *   ****************
       *
       *   * Hello World! *
       *
       *   ****************
       *
       *   where if you didn't guess the columns are labeled in hexadecimal, so
       *   that they are all single characters.  For the time being your message
       *   must be width-6 characters or shorter because I have mandated that at
       *   least one delimiter is printed on each side plus two spaces between that
       *   delimiter and the actual message. The remaining space will be filled by
       *   delimiters. Consequentially, this means your message must fit on one
       *   line.  This could of course be remedied by splitting the message
       *   across multiple lines, but I'm lazy.
       */
      void MakeBanner(const std::string& Message, const char delimiter='*', const
            int width=80);
};

template<typename T>
void PsiOutStream::Write2Buffer(const T& Input){
   Buffer_<<Input;
   this->DumpBuffer();
}


}//End psi namespace
#endif /* STREAMBASE_H_ */
