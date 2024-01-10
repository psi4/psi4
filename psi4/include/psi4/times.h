/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef psi_include_times_h
#define psi_include_times_h

#ifdef _MSC_VER
// Fake Windows implementation of the system/user timer
struct tms {
    double tms_stime;
    double tms_utime;
};
static void times(struct tms *time) {
    time->tms_stime = 0;
    time->tms_utime = 0;
}
#define _SC_CLK_TCK 0
static long sysconf(int name) { return (long)name; }
#else
#include <sys/times.h>
#include <unistd.h>
#endif

#endif