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
#ifndef PSIOUTSTREAM_H_
#define PSIOUTSTREAM_H_

#include "StreamBase.h"

namespace psi{

///Specialization of PsiStream to output streams
class PsiOutStream:public PsiStreamBase<std::ostream>{
   private:
      typedef PsiStreamBase<std::ostream> BaseType;

      ///Dumps Buffer to Stream_
      void Buffer2Stream();

      /** \brief This is the interface for Printf and << operator to
       *           Buffer_ for most writable objects
       *
       *  \param[in] Input An object that can be passed to an ostream
       */
      template<typename T>
      std::ostream& Write2Buffer(const T& Input);

      /** \brief This is the interface for Printf and << operator to
       *     Buffer_ for special iostream functions.
       *
       *     \param[in] fp iostream functions like std::endl
       */
      std::ostream& Write2Buffer(StreamManips fp);

   public:

      ///Makes an OutStreamBase that defaults to std::cout
      PsiOutStream(SharedOutStream Stream=SharedOutStream());

      ///Creates this by copying other, via base copy constructor
      PsiOutStream(const PsiOutStream& other):BaseType(other){}

      ///Allows assignment of other to this, calls base assignment
      const PsiOutStream& operator=(const PsiOutStream& other){
         BaseType::operator=(other);
         return *this;
      }

      ///No memory to free up (we don't delete cout or cerr)
      virtual ~PsiOutStream(){}

      ///Flushes the Stream
      void Flush(){
         if(this->ImSpecial()){
            Stream_->flush();
         }
      }

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
      std::ostream& operator<<(const T& Input){return Write2Buffer(Input);}

      ///Allows c++ like interface for things like std::endl
      std::ostream& operator<<(StreamManips fp){return Write2Buffer(fp);}

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
      void MakeBanner(const std::string& message, const char delimiter='*', const
            int width=80);
};


template<typename T>
std::ostream& PsiOutStream::Write2Buffer(const T& Input){
   Buffer_<<Input;
   this->Buffer2Stream();
   return Buffer_;
}

}//End namespace



#endif /* PSIOUTSTREAM_H_ */