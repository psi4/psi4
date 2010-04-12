#include <stdio.h>
#include <stdarg.h>
#include "parallel.h"

using namespace psi;

extern "C" {

int fprintf(FILE * __restrict __stream, const char * __restrict __format, ...)
{
    va_list args;
    va_start(args, __format);
    int status = 0;

    if (Communicator::world->me() == 0) {
        status = vfprintf(__stream, __format, args);
    }

    va_end(args);
    return status;
}

}
