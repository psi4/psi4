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
#ifndef SRC_LIB_LIBPSIUTIL_TABLE_H_
#define SRC_LIB_LIBPSIUTIL_TABLE_H_

namespace psi {

enum TableSide{TOP=0,RIGHT=1,BOTTOM=2,LEFT=3};

/** \brief Little class to hold a column's data
 *
 *  The class is designed to make a table out of existing data. It takes
 *  a pointer to the data (your data will only be accessed not modified)
 *  and an offset.  The offset allows your data to be stored in a matrix,
 *  such that only some columns are used. For example:
 *  \code
 *  double MyData[]={
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0,
 *     7.0, 8.0, 9.0
 *  }
 *  TableColumn<double> Col2("Column2",&MyData[1],3);
 *  \endcode
 *  would print (assuming we told the table to print 3 rows):
 *  \verbatim
 *  Column2
 *  2.0
 *  5.0
 *  8.0
 *  \endverbatim
 *
 *  The number of rows is an argument to the actual Table class and is
 *  handled there.  As you can see the pointer passed to the TableColumn
 *  class should point to the first element you want printed, if an offset
 *  is provided, say of \f$X\f$, then every \f$X\f$-th element after the
 *  original value will be printed (\f$X\f$ will be the number of columns
 *  in your matrix, unless you are doing something super fancy).
 */
template <typename T>
class TableColumn {
   private:
      ///A pointer to the first element of data we are printing
      const T* Data_;
      ///A descriptive name of the data we are printing
      std::string ColName_;
      ///The number of elements to pass over before printing one
      unsigned Offset_;
      ///The character that should appear below every element
      char VertDelim_;
      ///The character that should appear to the right of the element
      char HorDelim_;
      ///The character that will appear under the title
      char TitleDelim_;
   public:
      /** \brief Constructor for making a TableColumn
       *
       *   This constructor takes:
       *   \param[in] ColName The title of the column, to be printed in
       *              row 1
       *   \param[in] Data A A pointer to the first element you want printed
       *   \param[in] Offset Prints every \f$X\f$-th value
       *   \param[in] VertDelim The character printed between this row
       *                        and the row under it (ignored for last row)
       *   \param[in] HorDelim The character printed to the right of this
       *                       column (ignored for last column)
       *   \param[in] TitleDelim The character printed between the title
       *                         row and the first data row
       *
       *   For all delimeters passing a value of '\0' (the null character)
       *   causes the delimeter to not be printed.
       *
       *
       *
       */
      TableColumn<T>(const std::string& ColName, const T* Data,
            const unsigned Offset=1, const char VertDelim='\0',
            const char HorDelim=' ', const char TitleDelim='-');
      ///Returns the i-th element of the data set
      T operator[](const int i) const {return Data_[i*Offset_];}
      ///Returns the title delimiter
      char TitleDelim() const {return TitleDelim_;}
      ///Returns the horizontal delimiter (between cols)
      char HorDelim() const {return HorDelim_;}
      ///Returns the vertical delimiter (between rows)
      char VertDelim() const {return VertDelim_;}
      ///Returns the title
      std::string Name() const {return ColName_;}
};

///Forward declaration of our utility class
template<typename T> class TableData;

/** \brief The big deal. The actual table
 *
 *   This table class builds a table with a truly arbitrary number of
 *   columns.  Each column needs to be of type TableColumn<TYPE>, but
 *   TYPE need not be the same for each column.  The only restriction
 *   on TYPE is that it must be passable via the stream operator<< to a
 *   stringstream.
 *
 *   To use this class, say you have made four columns of data, the i-th
 *   one begin of type, Ti.  The type of this Table class is then:
 *   \code
 *   Table<TableColumn<T1>,TableColumn<T2>,TableColumn<T3>,TableColumn<T4>>
 *   \endcode
 *   Each column must have the same number of rows (pad them if they don't).
 *   Letting the number of rows be R, the constructor call is then:
 *   \code
 *   MyTable(R,C1,C2,C3,C4);
 *   \endcode
 *   where C1,C2,C3, and C4 are our instances of the four columns.  To
 *   print the table call GetTable(), which will return your table as
 *   a string.
 */
template <typename T, typename... Args>
class Table:protected Table<Args...> {
   private:
      ///The data for the current column
      TableData<T> MyData_;
   protected:
      ///Constructor for the recursion, only used internally
      Table<T, Args...>(const int NRows, const int NCols,
            const int MyCol,const T& Data, Args&... args);
      ///Fxn that actually prints out the this column
      std::string PrintOut(const int row) const;
   public:
      /** \brief Main constructor.  Builds an NRows table
       *
       *  The arguments to this class should be the number of rows,
       *  and a series of TableColumn<T> Objects (possibly with different
       *  T values).  The variadic nature allows for an arbitrary number of
       *  columns.
       *
       */
      Table<T, Args...>(const int NRows,const T& Data,Args&... args);
      ///Allows you to change the top, bottom, left and right border char
      void SetBorder(TableSide Side,const char s);
      ///Returns your table in a ready-to-print format
      std::string GetTable() const;
};

/************ Implementations **********/
template<typename T>
TableColumn<T>::
   TableColumn<T>(const std::string& ColName, const T* Data,
            const unsigned Offset, const char VertDelim,
            const char HorDelim, const char TitleDelim) :
            ColName_(ColName), Data_(Data), Offset_(Offset),
                  VertDelim_(VertDelim), HorDelim_(HorDelim),
                  TitleDelim_(TitleDelim) {}

/** \brief Class to hold the stuff common to each Table class
 *   to minimize copy/paste
 *
 *   This class is used internally within the Table class and should not
 *   be constructed anywhere else.  Ignore its existence please!!!!!
 */
template <typename T>
class TableData {
   private:
      ///The actual data
      const T& Data_;
      ///The column I'm responsible for
      int MyCol_;
      ///The number of rows in the table
      int NRows_;
      ///The number of columns in the table
      int NCols_;
      ///The character used for the top,right,bottom,left borders
      char BorderDelims_[4];
      ///The number of characters I get to print
      int MySize_;
      void CalcMySize();
   public:
      TableData<T>(const T& Data, const int NRows,
                   const int NCols, const int MyCol,
                   const char TopBorder='\0',const char BottomBorder='\0',
                   const char LeftBorder='\0',const char RightBorder='\0');
      std::string PrintOut(const int Row) const;
      int NRows()const{return NRows_;}
      void SetBorder(TableSide Side,const char s){
         BorderDelims_[Side]=s;
         CalcMySize();
      }
      char GetBorder(TableSide Side)const{return BorderDelims_[Side];}
};

template<typename T,typename...Args>
Table<T,Args...>::Table<T, Args...>(const int NRows, const int NCols,
      const int MyCol,const T& Data, Args&... args) :
   Table<Args...>(NRows, NCols, MyCol+1, args...),
    MyData_(Data, NRows, NCols, MyCol) {}


template<typename T,typename...Args>
std::string Table<T,Args...>::PrintOut(const int row) const {
   std::stringstream Result;
   Result<<MyData_.PrintOut(row)<<Table<Args...>::PrintOut(row);
   return Result.str();
}


template<typename T, typename...Args>
Table<T,Args...>::
Table<T, Args...>(const int NRows,const T& Data,Args&... args) :
            Table<Args...>(2*NRows+3, (int)sizeof...(Args)+1,1, args...),
                  MyData_(Data, 2*NRows+3, (int)sizeof...(Args)+1, 0) {}


template <typename T, typename... Args>
std::string Table<T,Args...>::GetTable()const{
   std::stringstream Result;
   for (int i=0; i<MyData_.NRows(); i++){
      std::string Temp=PrintOut(i);
      if(Temp!="")Result<<Temp<<std::endl;
   }
   return Result.str();
}

template<typename T, typename... Args>
void Table<T,Args...>::SetBorder(TableSide Side,const char s){
         Table<Args...>::SetBorder(Side,s);
         MyData_.SetBorder(Side,s);
}


///The base case of our Table (1 column of data)
template <typename T>
class Table<T> {
   private:
      TableData<T> MyData_;
   protected:
      Table<T>(const int NRows, const int NCols, const int MyCol,
            const T& Data) :
            MyData_(Data, NRows, NCols, MyCol) {
      }
      std::string PrintOut(const int row) const {
         return MyData_.PrintOut(row);
      }
      void SetBorder(TableSide Side,const char s){
         MyData_.SetBorder(Side,s);
      }
};





template <typename T>
TableData<T>::TableData<T>(const T&Data,
      const int NRows, const int NCols,
      const int MyCol, const char TopBorder,
      const char BottomBorder, const char LeftBorder,
      const char RightBorder) :
      Data_(Data), NRows_(NRows), NCols_(NCols), MyCol_(MyCol){
   BorderDelims_[TOP]=TopBorder;BorderDelims_[BOTTOM]=BottomBorder;
   BorderDelims_[LEFT]=LeftBorder;BorderDelims_[RIGHT]=RightBorder;
}

template<typename T>
void TableData<T>::CalcMySize(){
   int BorderChars=(GetBorder(LEFT)=='\0' ? 0 : 1);
          BorderChars+=(GetBorder(RIGHT)=='\0' ? 0 : 1);
      //There are NCols-1
      BorderChars+=NCols_-1;
      const int UsableChars=80-BorderChars;
      const int Remainder=(UsableChars)%NCols_;
      MySize_=(UsableChars-Remainder)/NCols_+(MyCol_==0?Remainder:0);
}

namespace TablePrint{
   std::string Repeat(const int N,const char c){
      std::stringstream Result;
      for(int i=0;i<N;i++)Result<<c;
      return Result.str();
   }

}

//This should be done with setw, but it's being a pain in my ass
template <typename T>
std::string TableData<T>::PrintOut(const int Row) const {
   std::stringstream Result;
   //The
   int Width=MySize_+(MyCol_==NCols_-1?0:1);
   if(MyCol_==0&&GetBorder(LEFT)!='\0')Result<<GetBorder(LEFT);
   //Top Border
   const bool OddRow=(Row%2==1),EvenRow=!OddRow;
   const bool FirstRow=(Row==0),ThirdRow=(Row==2),LastRow=(Row==NRows_-1);
   const bool HasVert=(Data_.VertDelim()!='\0'),
              HasTitle=(Data_.TitleDelim()!='\0');
   if (FirstRow && GetBorder(TOP)!='\0')
      Result<<TablePrint::Repeat(Width,GetBorder(TOP));
   else if (Row==1){//Col Titles
      Result<<std::setw(MySize_)<<(Data_.Name()==""?" ":Data_.Name());
   }
   else if (ThirdRow && HasTitle){ //Title bottom border
      Result<<TablePrint::Repeat(Width,Data_.TitleDelim());
   }
   else if (OddRow){ //Normal Data Row
      Result<<std::setw(MySize_)<<Data_[(Row-3)/2];
   }
   else if (EvenRow && !LastRow && HasVert){//Normal bottom delim
      Result<<std::setfill(Data_.VertDelim())<<std::setw(MySize_+1);
   }
   //The bottom of the table
   else if(LastRow &&GetBorder(BOTTOM)!='\0')
      Result<<TablePrint::Repeat(Width,GetBorder(BOTTOM));
   if(!LastRow && GetBorder(RIGHT)!='\0')Result<<GetBorder(RIGHT);
   else if(!FirstRow && !LastRow && MyCol_!=NCols_-1&&
           Data_.HorDelim()!='\0'){
      if((EvenRow && HasVert)|| OddRow)Result<<Data_.HorDelim();
   }
   return Result.str();
}

} //End namespace psi

#endif /* SRC_LIB_LIBPSIUTIL_TABLE_H_ */
