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
#ifndef TABLESPECS_H_
#define TABLESPECS_H_
#include<string>
#include<sstream>
#include<iomanip>

namespace psi {
///As the name implies this a NULL class
class Nothing {
};

///Enumerations of the possible alignments
enum Align {
   CENTER, LEFT, RIGHT
};

/** \brief TableSpecs will hold all the data needed to print a table out.
 *
 *  Because I'm being a good boy and not using C++11, this class suffers
 *  from two problems one it can't contain an arbitrary number of
 *  parameters and two every time a specialization is needed code has to
 *  be copied and pasted.  Because of the recursive nature
 *
 *  Say you want to have this class make a M column table, where each column
 *  has a data type of a different T and is labeled Tc where c is the column
 *  number, i.e. T1 is the type of column1, T2, is the type of column 2, etc.
 *  Furthermore assume that your table will have N rows. I mandate that the
 *  number of rows for each column must be the same, so use place holders if
 *  that isn't the case.  Your object definition is (a lot longer here
 *  because I have to show you all the options etc.):
 *
 *  \code
 *  ///This is the character that will be used to delimit the titles from
 *  ///the data (default='-')
 *  char TSep='*';
 *
 *  ///This is the separator between columns (default= a space)
 *  char ColSep='|';
 *
 *  ///This is the separator between rows of data (default= nothing)
 *  char RowSep="-";
 *
 *  ///This is the width of your table
 *  int Width=26;
 *
 *  ///This is your table specification
 *  TableSpecs<T1,T2,T3,...,TM> MyData(N,TSep,ColSep,RowSep,Width);
 *
 *  ///These are pointers to the first elements of each Data array
 *  T1* Data1;
 *  T2* Data2;
 *  ...
 *  TM* Data3;
 *
 *  MyData.Init(Data1,Data2,...,DataM);
 *
 *  // This is a vector of the column titles
 *  //
 *  // Each title may span more than one line.  If they do Titles should
 *  // be an array, such that if the longest title spans L lines it is of
 *  // L*M length, and Titles[i*M+j] is the i-th line of the j-th
 *  // title
 *  std::vector<std::string> Titles;
 *  Titles.push_back("Title1");
 *  Titles.push_back("Title2");
 *  ...
 *  Titles.push_back("TitleM");
 *
 *  MyData.SetTitles(Titles);//If not specified no titles printed
 *
 *  ///This is how you can control the alignments of each column
 *  ///Choices are CENTER,LEFT,RIGHT defaulting to CENTER
 *  std::vector<Align> Alignments;
 *  Alignments.push_back(Align1);
 *  Alignments.push_back(Align2);
 *  ...
 *  Alignments.push_back(Align3);
 *
 *  MyData.SetAlignments(Alignments);
 *
 *  ///This is how you get your Table, then you are free to do with it what
 *  ///you will
 *  std::string MyTable=MyData.Table();
 *  \endcode
 *
 *  The above code should generate something along the lines of (blank lines
 *  between rows so hopefully doxygen doesn't destroy it):
 *
 *  Column #:    0123456789ABCDEFGHIJKLMNOP (assumed 26 columns)
 *
 *                Title1 | Title2 | Title3
 *
 *               **************************
 *
 *               Data1[0]|Data2[0]|Data3[0]
 *
 *               --------------------------
 *
 *               Data1[1]|Data2[1]|Data3[1]
 *
 *               --------------------------
 *
 *               ...
 *
 *               --------------------------
 *
 *               Data1[N]|Data2[N]|Data3[N]
 *
 *
 *
 *
 *  In order to minimize the amount of duplicate code this class is defined
 *  recursively, which may make it hard to understand at first.
 *  Consequentially, I have included a pretty thorough description below
 *  if you are hell-bent on understanding it or modifying it read on.
 *
 *
 *  The line this comment is actually attached to is a forward declaration
 *  of a TableSpecs class with up to 6 Columns.  All of these classes
 *  together are make up the TableSpecs class and how it works is defined
 *  here.  If in the future more columns are needed simply add more template
 *  parameters and update the following classes so that they are consistent
 *  with the number of template parameters you changed it to.  You also will
 *  need to modify StreamBase.h's MakeTable function.
 *
 *  The following discussion assumes that
 *  the class is still built with 6 columns; if it is not the generalization
 *  to the current number of columns should be straightforward.
 *
 *  Only one implementation of a TableSpecs class is actually defined.  That
 *  is the case we call the general case and it corresponds to the case when
 *  6 arguments are specified.  Every case between 6 and 0 Columns is then
 *  filled in via recursion.  At each level of recursion we pull of the
 *  first template parameter and define it as a data set and pass the
 *  remaining ones to the base class.  This ends when we hit the class
 *  TableSpecs<Nothing,Nothing,Nothing,Nothing,Nothing,Nothing>, which we
 *  term the base case.
 *
 *  So how does this work more concretely, say we want three columns of
 *  integers, we then define an object TableSpecs<int,int,int>.  This expands
 *  to the general case of TableSpecs<int,int,int,Nothing,Nothing,Nothing>,
 *  which has a member std::vector<int> Data and a base class:
 *  TableSpces<int,int,Nothing,Nothing,Nothing,Nothing>.  This base class has
 *  a member std::vector<int> Data and a base class TableSpecs<int,Nothing...>
 *  where the ellipses continue for 4 more Nothings.  This base class also
 *  has a member std::vector<int> Data, but now it's base class is the base
 *  recursion case.  Ultimately we use typedefs to avoid the messy types that
 *  arise.
 *
 *  Let us assume in our example above are TableSpecs objects are of types:
 *  General, Generalm1 (General "minus" 1), Generalm2, and Base for the class
 *  hierarchy from TableSpecs<int,int,int> down to the base case respectively.
 *  Column 1's data is then General::Data, Column 2's data is then
 *  Generalm1::Data, and Column 3's data is then Generalm2::Data.
 *
 *  Also note that for the for the Table() function to work, each class
 *  needs to have access to the data.  The easiest way to do this is to
 *  use setters and getters to go to and from the base respectively.
 */
template <typename T1=Nothing, typename T2=Nothing, typename T3=Nothing,
      typename T4=Nothing, typename T5=Nothing, typename T6=Nothing>
class TableSpecs;

///Base case
template <>
class TableSpecs<Nothing, Nothing, Nothing, Nothing, Nothing, Nothing> {
   private:
      std::string TokenSep_;
   protected:

      ///Returns true because this is the base
      virtual bool IsBase() const {
         return true;
      }

      ///The number of rows
      int NRows_;
      virtual int NRows() const {
         return NRows_;
      }

      ///This is the title of the column, just stored in general class
      std::vector<std::string> Title_;
      virtual std::vector<std::string> Title() const {
         return Title_;
      }

      ///This is the alignment of each column
      std::vector<Align> Alignments_;
      virtual std::vector<Align> Alignments() const {
         return Alignments_;
      }

      ///This is the Title separator character, defaults to '-'
      char TitleSep_;
      virtual char TitleSep() const {
         return TitleSep_;
      }

      ///This is the column separator character, defaults to " "
      char ColSep_;
      virtual char ColSep() const {
         return ColSep_;
      }

      ///This is the separator between rows, defaults to null string
      char RowSep_;
      virtual char RowSep() const {
         return RowSep_;
      }

      ///This is the width (in columns) of the table, defaults to 80
      int Width_;
      virtual int Width() const {
         return Width_;
      }

      ///This is the user requested widths, value of 0 means calculate it
      std::vector<int> ReqWidth_;
      virtual std::vector<int> ReqWidth()const{
         return ReqWidth_;
      }
   public:
      ///Sets the titles
      void SetTitles(const std::vector<std::string>& Titles) {
         Title_=Titles;
      }

      ///Sets the alignment of the columns
      void SetAlignments(const std::vector<Align>& Alignments) {
         Alignments_=Alignments;
      }

      ///Sets the width of each column
      void SetWidths(const std::vector<int>& Widths){
         ReqWidth_=Widths;
      }
      ///Returns a NULL string
      virtual std::string GetData(int row) const {
         return std::string(" "+TokenSep_+" ");
      }

      ///Returns the number of columns
      virtual int NCols() const {
         return 0;
      }

      ///This fxn should never be called, but is needed so that it compiles
      void Init(void* Data1=NULL, void* Data2=NULL, void* Data3=NULL,
            void* Data4=NULL, void* Data5=NULL, void* Data6=NULL) const {
      }

      TableSpecs<Nothing, Nothing, Nothing, Nothing, Nothing, Nothing>(
            const int NRows, const char TitleSep='-', const char ColSep=' ',
            const char RowSep='\0', int Width=80) :
            NRows_(NRows), TitleSep_(TitleSep), ColSep_(ColSep),
                  RowSep_(RowSep), Width_(Width), TokenSep_("RMRStringSepRMR") {
      }

      ///Get rid of warning about non-virtual destructor...
      virtual ~TableSpecs<Nothing, Nothing, Nothing, Nothing, Nothing, Nothing>() {
      }
};





///General case, class that is instantiated any time all T's are not Nothing
template <typename T1, typename T2, typename T3, typename T4, typename T5,
      typename T6>
class TableSpecs:public TableSpecs<T2, T3, T4, T5, T6> {
   private:
      ///Typedef of the base for compactness
      typedef TableSpecs<T2, T3, T4, T5, T6> BaseType;
      std::string TokenSep_;
   protected:
      ///This class does not take ownership of your data
      T1* Data_;

      ///Tells us we are no longer the base recursion class
      virtual bool IsBase() const {
         return false;
      }

      ///Returns the number of columns
      virtual int NCols() const {
         int cols=BaseType::NCols();
         return ++cols;
      }

      ///Returns the title vector
      virtual std::vector<std::string> Title() const {
         return BaseType::Title();
      }

      ///Returns the Alignment vector
      virtual std::vector<Align> Alignments() const {
         return BaseType::Alignments();
      }

      ///Returns the width
      virtual int Width() const {
         return BaseType::Width();
      }

      ///Returns the RowSeparator
      virtual char RowSep() const {
         return BaseType::RowSep();
      }

      ///Returns the Column separator
      virtual char ColSep() const {
         return BaseType::ColSep();
      }

      ///Returns the Title separator
      virtual char TitleSep() const {
         return BaseType::TitleSep();
      }

      ///Returns the number of rows
      virtual int NRows() const {
         return BaseType::NRows();
      }

      ///Returns the requested width of the column
      virtual std::vector<int> ReqWidth()const{
         return BaseType::ReqWidth();
      }
   public:

      ///Sets the titles
      virtual void SetTitles(const std::vector<std::string>& Titles) {
         BaseType::SetTitles(Titles);
      }

      ///Sets the width of each column
      virtual void SetWidths(const std::vector<int>& Widths){
         BaseType::SetWidths(Widths);
      }

      ///Sets the alignment of the columns
      virtual void SetAlignments(const std::vector<Align>& Alignments) {
         BaseType::SetAlignments(Alignments);
      }

      virtual std::string Table()const;

      ///Default constructor, sets Data to NULL
      TableSpecs<T1, T2, T3, T4, T5, T6>(const int NRows, const char
         TitleSep='-', const char ColSep=' ', const char RowSep='\0',
         int Width=80) :
             Data_(NULL), TokenSep_("RMRStringSepRMR"),
             TableSpecs<T2, T3, T4, T5, T6>(NRows, TitleSep, ColSep, RowSep,
                   Width) {}

      virtual ~TableSpecs<T1, T2, T3, T4, T5, T6>() {
      }

      ///The call that actually initializes this class
      void Init(T1* Data1, T2* Data2=NULL, T3* Data3=NULL, T4* Data4=NULL,
            T5* Data5=NULL, T6* Data6=NULL) {
         Data_=Data1;
         if (Data2!=NULL) BaseType::Init(Data2, Data3, Data4, Data5, Data6);
      }

      ///Overloaded function to catch doubles
      template<typename NotDouble>
      std::string PrintData(const NotDouble Data)const{
         std::stringstream RawData;
         RawData<<Data<<" "<<TokenSep_<<" ";
         return RawData.str();
      }


      std::string PrintData(const double Data)const{
         std::stringstream RawData;
         RawData<<std::fixed<<std::setprecision(16)<<Data<<" "<<TokenSep_<<" ";
         return RawData.str();
      }

      ///Returns the Data for row "row" as a string (non-EOL terminated)
      virtual std::string GetData(const int row) const {
         std::stringstream RowData;
         RowData<<PrintData(Data_[row]);
         if (!BaseType::IsBase()) RowData<<BaseType::GetData(row);
         return RowData.str();
      }


};

static inline std::string Alignment(const int Width, const std::string Text,
      const Align& AlignType=CENTER) {
   std::string ReturnString;
   int TitleSize=Text.size();
   if (TitleSize>Width)      //Truncation of title, sorry
   ReturnString=Text.substr(0, Width);
   else {
      int nspaces=Width-TitleSize;
      int lspaces=0;
      int rspaces=0;
      if (AlignType==CENTER) {
         lspaces=(nspaces-nspaces%2)/2+(nspaces%2);
         rspaces=nspaces-lspaces;
      }
      else {
         lspaces=(AlignType==LEFT ? 0 : nspaces);
         rspaces=(AlignType==RIGHT ? 0 : nspaces);
      }
      std::string Lspaces(lspaces, ' ');
      std::string Rspaces(rspaces, ' ');
      ReturnString=Lspaces;
      ReturnString+=Text;
      ReturnString+=Rspaces;
   }
   return ReturnString;
}

static inline std::string ParseSStream(std::stringstream& tokenizer,
      const std::string& TokenSep) {
   std::string line;
   bool done=false;
   int itr=0;
   while (!done&&!tokenizer.eof()) {
      std::string token;
      tokenizer>>token;
      if (token!=TokenSep){
         if(itr++>=1)line+=" ";
         line+=token;
      }
      else done=true;
   }
   return line;

}

static inline std::string MakeRow(const int ncols, std::stringstream& tokenizer,
      const std::vector<Align>& Aligns, const std::vector<int>& ColspCell,
      const char colsep, const std::string& TokenSep) {
   std::stringstream DaTable;
   for (int j=0; j<ncols-1; j++) {
      std::string line=ParseSStream(tokenizer,TokenSep);
      Align align=(Aligns.size()>0 ? Aligns[j] : CENTER);
      DaTable<<Alignment(ColspCell[j], line, align)<<colsep;
   }
   std::string line=ParseSStream(tokenizer,TokenSep);
   Align align=(Aligns.size()>0 ? Aligns.back() : CENTER);
   DaTable<<Alignment(ColspCell.back(), line, align)<<std::endl;
   return DaTable.str();
}


template <typename T1, typename T2, typename T3, typename T4, typename T5,
      typename T6>
std::string TableSpecs<T1, T2, T3, T4, T5, T6>::Table()const{///Returns the Table as a string
                std::stringstream DaTable;//This will be the table

                ///Start with the number of rows
                int ncols=NCols();
                ///Then the widths the user set up
                std::vector<int> ColspCell=ReqWidth();
                std::vector<bool> OrigSet(ncols,false);
                int tempsize=ColspCell.size();
                if(tempsize<ncols){
                   for(int i=0;i<(ncols-tempsize);i++)
                      ColspCell.push_back(0);
                }
                int nColsByUser=0;
                int accountedforcolumns=0;
                for(int i=0;i<ncols;i++){
                   if(ColspCell[i]!=0){
                      OrigSet[i]=true;
                      nColsByUser++;
                      accountedforcolumns+=ColspCell[i];
                   }
                }
                int width=Width();
                int nseparators=ncols-1;
                int EffWidth=width-nseparators-accountedforcolumns;
                ncols-=nColsByUser;
                if(ncols==0&&EffWidth!=0)throw
                 PSIEXCEPTION("You micro-managed the columns and didn't count"
                       "them out right, noob...");
                int LeftOverCols=EffWidth%ncols;
                int othercolumns=(EffWidth-LeftOverCols)/ncols;
                ncols=NCols();
                for(int i=0;i<ncols;i++){
                   if(!OrigSet[i])ColspCell[i]=othercolumns;
                }
                for (int i=0; i<LeftOverCols; i++)
                   if(!OrigSet[i])ColspCell[i]++;

                //Add the titles
                std::vector<std::string>Titles=Title();
                std::vector<Align>Aligns=Alignments();
                int NTitleRows=Titles.size()/ncols;
                char colsep=ColSep();
                char rowsep='\0';
                std::stringstream Buffer;
                for (int i=0; i<Titles.size(); i++) {
                   Buffer<<Titles[i];
                   Buffer<<" "<<TokenSep_<<" ";
                }

                //Titles
                for (int i=0; i<NTitleRows; i++)
                   DaTable<<MakeRow(ncols, Buffer, Aligns, ColspCell, colsep, TokenSep_);

                //Title seperator
                std::string FirstSep(width, TitleSep());
                DaTable<<FirstSep<<std::endl;

                //Data
                int nrows=NRows();
                rowsep=RowSep();
                for (int i=0; i<nrows; i++) {
                   std::stringstream tokenizer(GetData(i));
                   DaTable<<MakeRow(ncols, tokenizer, Aligns, ColspCell, colsep, TokenSep_);
                   if (rowsep!='\0') {
                      std::string RowSep(width, rowsep);
                      DaTable<<RowSep<<std::endl;
                   }
                }
                return DaTable.str();
             }


}      //End namespace psi
#endif /* TABLESPECS_H_ */