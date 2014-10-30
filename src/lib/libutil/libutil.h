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

#ifndef _psi_src_lib_libutil_libutil_h_
#define _psi_src_lib_libutil_libutil_h_

#include <string>
#include <vector>
#include <sys/time.h>

namespace psi {

typedef std::vector<std::string> strvec;

std::string file_to_string(std::string const& name);


bool space(char c);
bool not_space(char c);
std::vector<std::string> split(const std::string& str);
std::vector<std::string> split_indices(const std::string& str);
void to_lower(std::string& str);
void to_upper(std::string& str);
std::string to_string(const int val);
std::string to_string(const double val);

double to_double(const std::string str);

double to_double(const std::string str);
double ToDouble(const std::string inString);
int string_to_integer(const std::string inString);
void append_reference(std::string& str, int reference);
std::string add_reference(std::string& str, int reference);
void append_reference(std::string& str, int reference);

std::string find_and_replace(std::string & source, const std::string & target, const std::string & replace);
void trim_spaces(std::string& str);

/**
 * @brief Compute the Levenshtein distance between two strings
 * @param s1
 * @param s2
 * @return the distance as an unsigned integer
 */
unsigned int edit_distance(const std::string& s1, const std::string& s2);


class Timer
{
public:
    Timer() : ___start(), ___end(), ___dummy(),
              delta_time_seconds(0), delta_time_hours(0), delta_time_days(0)
        {gettimeofday(&___start,&___dummy);}
    double get() {gettimeofday(&___end,&___dummy);
      delta_time_seconds=(___end.tv_sec - ___start.tv_sec) + (___end.tv_usec - ___start.tv_usec)/1000000.0;
      delta_time_hours=delta_time_seconds/3600.0;
      delta_time_days=delta_time_hours/24.0;
      return(delta_time_seconds);}
private:
  struct timeval ___start,___end;
  struct timezone ___dummy;
  double delta_time_seconds;
  double delta_time_hours;
  double delta_time_days;
};

void generate_combinations(int n, int k, std::vector<std::vector<int> >& combinations);

}

#endif // _psi_src_lib_libutil_libutil_h_
