#include "class.h"

using namespace yeti;
using namespace std;

const char* TypeInfo<double>::printf_str = "%12.8f";
const char* TypeInfo<float>::printf_str = "%12.8f";
const char* TypeInfo<quad>::printf_str = "%20.14Lf";
const char* TypeInfo<int>::printf_str = "%6d";

float TestEquals<float>::cutoff = 1e-12;
double TestEquals<double>::cutoff = 1e-12;
quad TestEquals<quad>::cutoff = 1e-12;

