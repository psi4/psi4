/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#ifndef STRINGARRAY_H_
#define STRINGARRAY_H_
#include <boost/shared_ptr.hpp>
#include <map>
#include <sstream>
#include "exception.h"
#include "psi4-dec.h"
#include "../libparallel2/Communicator.h"
#include "../libparallel2/ParallelEnvironment.h"
namespace psi{

///Handy fxn for accessing the members of a map, as if they were const
template<typename T1,typename T2>
T2 MapDeRef(const std::map<T1,T2>& Map,const T1& Key){
   T2 returnvalue;
   if(Map.count(Key)==1)
      returnvalue=(const_cast<std::map<T1,T2>* >(&Map))->operator[](Key);
   else{
      std::stringstream Error;
      Error<<"Your key: "<<Key<<" doesn't exist"<<std::endl;
      throw PSIEXCEPTION(Error.str().c_str());
   }
   return returnvalue;
}

class PsiIOStringArray{
   private:
      void Copy(const PsiIOStringArray& other){
         this->Lengths_=other.Lengths_;
      }

   protected:

      typedef boost::shared_ptr<int[]> SharedInt;

      ///The lengths of each array
      std::map<std::string,SharedInt> Lengths_;

   public:
      ///Use pointers, to allow coupling
      void SetLength(const std::string& variable,const SharedInt& Length){
         Lengths_[variable]=Length;
      }

      ///Returns the length, as an integer
      int GetLength(const std::string& value)const{
         SharedInt temp=MapDeRef(Lengths_,value);
         return temp[0];
      }

      PsiIOStringArray(const PsiIOStringArray& other){
         this->Copy(other);
      }

      PsiIOStringArray(){}
      virtual ~PsiIOStringArray(){}
      const PsiIOStringArray& operator=(const PsiIOStringArray& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }

      virtual void Read(boost::shared_ptr<PSIO> psio,
            const std::string& Name,const int FileNum)=0;

      virtual  void Write(const boost::shared_ptr<PSIO> psio,
            const std::string& Name,const int FileNum)const=0;

      virtual void Broadcast(const std::string& Comm, const int Host,
            const std::string& Name)const=0;

      virtual void Receive(const std::string& Comm, const int Host,
            const std::string& Name)=0;

};


template<typename T>
class PsiIOSAImpl:public PsiIOStringArray{
   private:
      typedef boost::shared_ptr<T[]> SharedType;
      typedef PsiIOSAImpl<T> ThisType;

      ///Actual copy
      void Copy(const ThisType& other){
         this->Members_=other.Members_;
      }

   protected:
      ///The actual arrays we are responsible for
      std::map<std::string,SharedType> Members_;

   public:

      ///Avoids deep copy
      void SetValue(const std::string& Name,const SharedType& values){
         if(!Lengths_[Name])throw PSIEXCEPTION("Please set length first");
         Members_[Name]=values;
      }

      ///Deep copy on data
      void SetValue(const std::string& Name,const T* values){
         if(!Lengths_[Name])throw PSIEXCEPTION("Please set length first");
         int Length=GetLength(Name);
         //Avoid reallocting things of length 1 as they may couple
         //Need to give them memory if they don't have it though...
         if(Length!=1||!Members_[Name])
            Members_[Name]=SharedType(new T [Length]);
         memcpy(Members_[Name].get(),values,Length*sizeof(T));
      }

      SharedType GetValue(const std::string& value)const{
         return MapDeRef(Members_,value);
      }

      void Write(const boost::shared_ptr<PSIO> psio,
            const std::string& Name,
            const int FileNum)const{
         int Length=GetLength(Name);
         SharedType temp=GetValue(Name);
         psio->write_entry(FileNum,Name.c_str(),(char *)temp.get(),
                           sizeof(T)*Length);
      }

      void Read(boost::shared_ptr<PSIO> psio,
            const std::string& Name,
            const int FileNum){
         int Length=GetLength(Name);
         if(Length!=1||!Members_[Name])
            Members_[Name]=SharedType (new T[Length]);
         SharedType temp=GetValue(Name);
         psio->read_entry(FileNum,Name.c_str(),
                  (char *)(temp.get()),sizeof(T)*Length);
      }

      void Broadcast(const std::string& CommIn, const int Host,
            const std::string& Name)const{
         int Length=GetLength(Name);
         boost::shared_ptr<const LibParallel::Communicator> Comm=WorldComm->GetComm();
         Comm->Bcast(&Length,1,Host);
         T* temp=(GetValue(Name).get());
         Comm->Bcast(temp,Length,Host);
      }

      void Receive(const std::string& CommIn, const int Host,
            const std::string& Name){
         int Length=0;
         boost::shared_ptr<const LibParallel::Communicator> Comm=WorldComm->GetComm();
         Comm->Bcast(&Length,1,Host);
         SharedInt temp(new int [1]);
         temp[0]=Length;
         SetLength(Name,temp);
         SharedType temp2(new T[Length]);
         Comm->Bcast(temp2.get(),Length,Host);
         SetValue(Name,temp2);
      }

      const ThisType& operator=(const ThisType& other){
         PsiIOStringArray::operator=(other);
         if(this!=&other)this->Copy(other);
         return *this;
      }

      PsiIOSAImpl<T>(){}

      PsiIOSAImpl<T>(const ThisType& other):PsiIOStringArray(other){
         this->Copy(other);
      }
};


}




#endif /* STRINGARRAY_H_ */