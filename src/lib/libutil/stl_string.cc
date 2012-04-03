#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <sstream>

#include "libutil.h"

namespace psi {

using namespace std;

/*********************
 String manipulation
 ********************/

std::vector<std::string> split(const std::string& str){
  // Split a string
  typedef std::string::const_iterator iter;
  std::vector<std::string> splitted_string;
  iter i = str.begin();
  while(i != str.end()){
    // Ignore leading blanks
    i = find_if(i,str.end(), not_space);
    // Find the end of next word
    iter j = find_if(i,str.end(),space);
    // Copy the characters in [i,j)
    if(i!=str.end())
      splitted_string.push_back(std::string(i,j));
    i = j;
  }
  return(splitted_string);
}

bool opening_square_bracket(char c);
bool closing_square_bracket(char c);

std::vector<std::string> split_indices(const std::string& str){
  // Split a string
  typedef string::const_iterator iter;
  strvec splitted_string;
  iter i = str.begin();
  while(i != str.end()){
    // Ignore leading blanks
    i = find_if(i,str.end(), opening_square_bracket);
    // Find the end of next word
    iter j = find_if(i,str.end(),closing_square_bracket);
    // Copy the characters in [i,j]
    if(i!=str.end())
      splitted_string.push_back(std::string(i,j+1));
    i = j;
  }
  return(splitted_string);
}

bool opening_square_bracket(char c)
{
  return (c=='[');
}

bool closing_square_bracket(char c)
{
  return (c==']');
}

bool space(char c)
{
  return isspace(c);
}

bool not_space(char c)
{
  return !isspace(c);
}


std::string find_and_replace(std::string & source, const std::string & target, const std::string & replace )
{
  string str = source;
  string::size_type pos = 0;   // where we are now
  string::size_type found;     // where the found data is

  if (target.size () > 0)   // searching for nothing will cause a loop
    {
    while ((found = str.find (target, pos)) != string::npos)
      {
      str.replace (found, target.size (), replace);
      pos = found + replace.size ();
      }
    }
  return str;
}

void trim_spaces(std::string& str)
{
  // Trim Both leading and trailing spaces
  size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
  size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

  // if all spaces or empty return an empty string
  if(( string::npos == startpos ) || ( string::npos == endpos))
  {
    str = "";
  }
  else
    str = str.substr( startpos, endpos-startpos+1 );
}


/*****************
 String conversion
 *****************/

void to_lower(std::string& str)
{
  std::transform( str.begin(), str.end(), str.begin(),::tolower);
}

void to_upper(std::string& str)
{
  std::transform( str.begin(), str.end(), str.begin(),::toupper);
}

double ToDouble(const std::string str)
{
  return std::strtod( str.c_str(), NULL);
}

double to_double(const std::string str)
{
  return std::atof(str.c_str());
}

std::string to_string(const int val)
{
    std::stringstream strm;
    strm <<  val;
    return strm.str();
}

std::string to_string(const double val)
{
    std::stringstream strm;
    strm << setprecision(25) << setw(35)  << val;
    return strm.str();
}

int string_to_integer(const std::string inString)
{
  int i = 0;
  char* end;
  i = static_cast<int>(std::strtod( inString.c_str(), &end ));
  return i;
}

std::string add_reference(std::string& str, int reference)
{
  return(str + "{" + to_string(reference) + "}");
}

void append_reference(std::string& str, int reference)
{
  str += "{" + to_string(reference) + "}";
}

/*********************************************************
 Memory Allocation
 *********************************************************/

/**
 * Convert the size of a doubles array in Mb using the definition 1Mb = 1048576 bytes
 * @param n size of the array
 * @return
 */
double to_MB(size_t n)
{
  return(double(n*sizeof(double))/1048576.0);
  // Using this definition 1 Mb has ca. 5% more than 1000000 bytes
}

unsigned long int init_smatrix(short**& matrix,int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = static_cast<unsigned long int>(size1);
  uli_size2 = static_cast<unsigned long int>(size2);
  size=uli_size1*uli_size2;
  if(!uli_size1 || !uli_size2){
    matrix=NULL;
  }else{
    matrix = new short*[uli_size1];
    short* vector = new short[size];
    for(unsigned long int i=0;i<size;i++) vector[i]=0;
    for(unsigned long int i=0;i<uli_size1;i++)
      matrix[i]=&(vector[i*uli_size2]);
  }
  return(size*sizeof(short));
}

unsigned long int free_smatrix(short**& matrix, int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = static_cast<unsigned long int>(size1);
  uli_size2 = static_cast<unsigned long int>(size2);
  size=uli_size1*uli_size2;
  if(matrix == NULL) return(0);
  delete[] matrix[0];
  delete[] matrix;
  return(size*sizeof(short));
}

unsigned long int init_smatrix(short***& matrix,int size1, int size2, int size3)
{
  unsigned long int size = static_cast<unsigned long int>(size1*size2*size3);
  matrix = new short**[size1];
  for(int i=0;i<size1;i++){
      matrix[i] = new short*[size2];
  }
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          matrix[i][j] = new short[size3];
      }
  }
  return(size*sizeof(short));
}

unsigned long int free_smatrix(short*** matrix, int size1, int size2, int size3)
{
  unsigned long int size = static_cast<unsigned long int>(size1*size2*size3);
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          delete[] matrix[i][j];
      }
  }
  for(int i=0;i<size1;i++){
      delete[] matrix[i];
  }
  delete[] matrix;
  return(size*sizeof(short));
}

}
