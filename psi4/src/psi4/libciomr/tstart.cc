/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
** \file
** \brief Controls starting and stopping of timers
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>

#ifdef _MSC_VER
#include <Winsock2.h>
#include <winsock.h>
#else
#include <unistd.h>
#endif

#include "psi4/psi4-dec.h"
#include "psi4/times.h"

#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

std::time_t time_start, time_end;
std::time_t time_start_overall;
int running = 0;
double user_start, sys_start;
double user_start_overall, sys_start_overall;
double user_stop, sys_stop;

/*!
** tstart(): Starts a timer
**
** \ingroup CIOMR
*/
void PSI_API tstart() {
    int error;
    struct tms total_tmstime;
    const long clk_tck = sysconf(_SC_CLK_TCK);
    times(&total_tmstime);

    // host name has up to HOST_NAME_MAX(==64) + 1 bytes, and must end in the null byte
    std::vector<char> name(65);
    error = gethostname(name.data(), 65);
    if (error != 0) strncpy(name.data(), "nohostname", 11);
    if (name.back() != '\0') name.push_back('\0');

    /// start a global timer
    if (!running) {
        time_start_overall = std::time(nullptr);
        user_start_overall = ((double)total_tmstime.tms_utime) / clk_tck;
        sys_start_overall = ((double)total_tmstime.tms_stime) / clk_tck;
        running = 1;
    }

    /// start module timers
    time_start = std::time(nullptr);
    user_start = ((double)total_tmstime.tms_utime) / clk_tck;
    sys_start = ((double)total_tmstime.tms_stime) / clk_tck;

    outfile->Printf("\n*** tstart() called on %s\n", name.data());
    outfile->Printf("*** at %s\n", ctime(&time_start));
}

/*!
** tstop(): Stop timer
**
** \ingroup CIOMR
*/
void PSI_API tstop() {
    int error;
    std::time_t total_time;
    std::time_t total_time_overall;
    struct tms total_tmstime;
    double user_s, sys_s;

    // host name has up to HOST_NAME_MAX(==64) + 1 bytes, and must end in the null byte
    std::vector<char> name(65);
    error = gethostname(name.data(), 65);
    if (error != 0) strncpy(name.data(), "nohostname", 11);
    if (name.back() != '\0') name.push_back('\0');

    time_end = std::time(nullptr);
    total_time = time_end - time_start;
    total_time_overall = time_end - time_start_overall;

    times(&total_tmstime);
    const long clk_tck = sysconf(_SC_CLK_TCK);
    user_stop = ((double)total_tmstime.tms_utime) / clk_tck;
    sys_stop = ((double)total_tmstime.tms_stime) / clk_tck;

    user_s = user_stop - user_start;
    sys_s = sys_stop - sys_start;

    outfile->Printf("\n*** tstop() called on %s at %s", name.data(), ctime(&time_end));

    /// print all module timings
    outfile->Printf("Module time:\n");
    outfile->Printf("\tuser time   = %10.2f seconds = %10.2f minutes\n", user_s, user_s / 60.0);
    outfile->Printf("\tsystem time = %10.2f seconds = %10.2f minutes\n", sys_s, sys_s / 60.0);
    outfile->Printf("\ttotal time  = %10d seconds = %10.2f minutes\n", (int)total_time, ((double)total_time) / 60.0);

    user_s = user_stop - user_start_overall;
    sys_s = sys_stop - sys_start_overall;

    /// print all overall timings
    outfile->Printf("Total time:\n");
    outfile->Printf("\tuser time   = %10.2f seconds = %10.2f minutes\n", user_s, user_s / 60.0);
    outfile->Printf("\tsystem time = %10.2f seconds = %10.2f minutes\n", sys_s, sys_s / 60.0);
    outfile->Printf("\ttotal time  = %10d seconds = %10.2f minutes\n", (int)total_time_overall,
                    ((double)total_time_overall) / 60.0);
}
}  // namespace psi
