/* SPDX-License-Identifier: BSD-3-Clause */
/* Psi4's libpsio, reimplemented on libspill.
 *
 * The signatures below are Psi4's, unchanged. That is the whole point: §7 of
 * DESIGN.md says the API is wrong if a port needs more than a mechanical
 * translation of call sites, so this port changes none of them. Psi4 keeps
 * psio.h and psio.hpp as they are; what goes away is the ~3.9k lines behind
 * them.
 *
 * Three facts from the Psi4 tree at a0e6ba5c4 make this possible, and each was
 * checked rather than assumed:
 *
 *   1. `psio_address` is linear. psio_get_address(start, shift) is exactly
 *      start + shift in (page, offset) coordinates, so page*PSIO_PAGELEN +
 *      offset is a faithful byte offset. No consumer does arithmetic on .page:
 *      zero occurrences outside libpsio.
 *   2. `psio_tocscan` is an existence test. Its psio_tocentry* result is never
 *      dereferenced outside libpsio -- ->sadd and ->eadd appear zero times --
 *      so the entry the shim hands back need only be valid and correctly
 *      filled, not a node in a live linked list.
 *   3. The on-disk table of contents is nobody's business. rd_toclen, tocread
 *      and toclen have no consumers outside libpsio; tocwrite has one. So
 *      libspill owning the table of contents costs Psi4 nothing.
 */
#ifndef PSIO_LIBSPILL_H
#define PSIO_LIBSPILL_H

#include <string>

/* Two ways to build this shim.
 *
 * By default it uses psio_types.h, a copy of Psi4's public types, so the port
 * can be compiled and tested without a Psi4 tree.
 *
 * With -DPSIO_USE_PSI4_HEADERS it includes Psi4's own psio.h instead and
 * declares nothing that Psi4 already declares. That is the conformance build:
 * every entry point below is then checked against Psi4's real declaration
 * rather than against our copy of it, and a signature that has drifted shows up
 * as a link error rather than at integration time. See
 * test/psio_conformance_test.cc. */
#ifdef PSIO_USE_PSI4_HEADERS
#include "psi4/libpsio/psio.h"
#else
#include "psio_types.h"
#endif

namespace psi {

#ifndef PSIO_USE_PSI4_HEADERS
psio_address psio_get_address(psio_address start, size_t shift);
psio_address psio_get_global_address(psio_address entry_start, psio_address rel_address);

void psio_error(size_t unit, size_t errval, std::string prev_msg = "");

int psio_init();

int  psio_open(size_t unit, int status);
int  psio_close(size_t unit, int keep);
int  psio_open_check(size_t unit);

int psio_write(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd, psio_address *eadd);
int psio_read(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd, psio_address *eadd);
int psio_write_entry(size_t unit, const char *key, char *buffer, size_t size);
int psio_read_entry(size_t unit, const char *key, char *buffer, size_t size);

psio_tocentry *psio_tocscan(size_t unit, const char *key);
bool psio_tocentry_exists(size_t unit, const char *key);
void psio_tocprint(size_t unit);
int  psio_tocwrite(size_t unit);
size_t psio_rd_toclen(size_t unit);
#endif

/* Not declared by Psi4's psio.h: psio_done has no counterpart there, and tocdel
 * and zero_disk are PSIO class methods rather than free functions. */
int psio_done();
bool psio_tocdel(size_t unit, const char *key);
void psio_zero_disk(size_t unit, const char *key, size_t rows, size_t cols);

/* Not part of Psi4's API: lets the test point the shim at a scratch directory,
 * a job PSIOManager does in Psi4 itself. */
void psio_set_scratch_dir(const char *dir);

}  // namespace psi
#endif
