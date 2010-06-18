#include <stdio.h>
#include <stdarg.h>
#include "parallel.h"

using namespace psi;

void p_fprintf(FILE * __restrict __stream, const char * __restrict __format, ...)
{
    va_list args;
    va_start(args, __format);
    int status = 0;
    
    status = vfprintf(__stream, __format, args);

    va_end(args);
}

