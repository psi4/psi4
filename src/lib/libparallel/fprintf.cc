#include <stdio.h>
#include <stdarg.h>
#include "parallel.h"

using namespace psi;

extern "C" {

    void p_fprintf(FILE * __restrict __stream, const char * __restrict __format, ...)
    {
        va_list args;
        va_start(args, __format);
        int status = 0;

        status = vfprintf(__stream, __format, args);

        va_end(args);
    }

    int fprintf(FILE * __restrict __stream, const char * __restrict __format, ...)
    {
        va_list args;
        va_start(args, __format);
        int status = 0;

        if (Communicator::world.get() != NULL && Communicator::world->me() == 0) {
            status = vfprintf(__stream, __format, args);
        }
        else if (Communicator::world.get() == NULL) {
            p_fprintf(stderr, "Communicator object does not exist.\n");
            status = vfprintf(__stream, __format, args);
        }

        va_end(args);
        return status;
    }

}
