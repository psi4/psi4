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
