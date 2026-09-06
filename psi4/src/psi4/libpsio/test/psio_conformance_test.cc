// SPDX-License-Identifier: BSD-3-Clause
/* Conformance check: does the shim implement Psi4's ACTUAL declarations?
 *
 * §6b measured the port against psio_types.h, our copy of Psi4's public types.
 * A copy can drift. This translation unit includes Psi4's own psio.h and takes
 * the address of every entry point the shim is meant to provide, through a
 * pointer whose type comes from that header. Linking it against the shim
 * therefore fails -- with an undefined reference -- if any definition's
 * signature does not match Psi4's declaration exactly.
 *
 * Not checked here, because the shim does not replace them: decode_errno,
 * psio_compose_err_msg and psio_getpid stay Psi4's. (psio_volseek is declared
 * in psio.h but its definition went with the I/O core and it has no callers --
 * see PORTING.md.)
 *
 * ---------------------------------------------------------------------------
 * The free functions are NOT the interface Psi4 uses.
 *
 * Issue #1: an earlier version of this file checked only the 16 free functions
 * and reported "16 entry points match", while the shim defined zero PSIO class
 * methods -- and the class is what consumers call, roughly 2400 times. The
 * check passed because it validated exactly the half of the API the shim
 * happened to implement.
 *
 * So the class is checked too, by taking the address of every method Psi4
 * declares. A method that is declared but not defined is a link error here,
 * which is the only way this can be caught before an adopter's build.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <sys/stat.h>

#include "psi4/libpsio/psio.h"
/* psio.hpp declares the PSIO class. The original version of this file included
 * only psio.h, which is why an entire half of the API went unchecked. */
#include "psi4/libpsio/psio.hpp"

namespace {

/* Each initialiser both names the psi4-declared type and demands our symbol. */
int         (*p_init)()                                              = &psi::psio_init;
int         (*p_open)(size_t, int)                                   = &psi::psio_open;
int         (*p_close)(size_t, int)                                  = &psi::psio_close;
int         (*p_open_check)(size_t)                                  = &psi::psio_open_check;

psi::psio_address (*p_get_address)(psi::psio_address, size_t)        = &psi::psio_get_address;
psi::psio_address (*p_get_global)(psi::psio_address, psi::psio_address)
                                                                     = &psi::psio_get_global_address;

int (*p_write)(size_t, const char *, char *, size_t, psi::psio_address, psi::psio_address *)
                                                                     = &psi::psio_write;
int (*p_read)(size_t, const char *, char *, size_t, psi::psio_address, psi::psio_address *)
                                                                     = &psi::psio_read;
int (*p_write_entry)(size_t, const char *, char *, size_t)           = &psi::psio_write_entry;
int (*p_read_entry)(size_t, const char *, char *, size_t)            = &psi::psio_read_entry;

psi::psio_tocentry *(*p_tocscan)(size_t, const char *)               = &psi::psio_tocscan;
bool   (*p_tocexists)(size_t, const char *)                          = &psi::psio_tocentry_exists;
void   (*p_tocprint)(size_t)                                         = &psi::psio_tocprint;
int    (*p_tocwrite)(size_t)                                         = &psi::psio_tocwrite;
size_t (*p_rd_toclen)(size_t)                                        = &psi::psio_rd_toclen;

void (*p_error)(size_t, size_t, std::string)                         = &psi::psio_error;

/* ---- the PSIO class: what Psi4 actually calls ----
 *
 * Comparing a member pointer against nullptr is not enough: &P::f is never
 * null, so the optimiser folds the test and never emits a reference to the
 * symbol. The address has to be materialised. touch() copies the pointer's
 * bytes into a volatile sink, which forces it. */
using P = psi::PSIO;

volatile unsigned long sink = 0;

template <class T>
void touch(T mp) {
    unsigned char b[sizeof(T)];
    std::memcpy(b, &mp, sizeof mp);
    unsigned long v = 0;
    for (size_t i = 0; i < sizeof b; i++) v = v * 31u + b[i];
    sink = sink + v;
}

void   (P::*m_open)(size_t, int)                                     = &P::open;
void   (P::*m_close)(size_t, int)                                    = &P::close;
int    (P::*m_open_check)(size_t) const                              = &P::open_check;
bool   (P::*m_exists)(size_t) const                                  = &P::exists;
void   (P::*m_rehash)(size_t)                                        = &P::rehash;

void   (P::*m_read)(size_t, const char *, char *, size_t, psi::psio_address, psi::psio_address *)
                                                                     = &P::read;
void   (P::*m_write)(size_t, const char *, char *, size_t, psi::psio_address, psi::psio_address *)
                                                                     = &P::write;
void   (P::*m_read_entry)(size_t, const char *, char *, size_t)      = &P::read_entry;
void   (P::*m_write_entry)(size_t, const char *, char *, size_t)     = &P::write_entry;
void   (P::*m_rw)(size_t, char *, psi::psio_address, size_t, int)    = &P::rw;

void   (P::*m_zero_disk)(size_t, const char *, size_t, size_t)       = &P::zero_disk;
void   (P::*m_zero_disk_s)(size_t, const std::string &, size_t, size_t) = &P::zero_disk;

psi::psio_tocentry *(P::*m_tocscan)(size_t, const char *)            = &P::tocscan;
bool   (P::*m_tocentry_exists)(size_t, const char *)                 = &P::tocentry_exists;
bool   (P::*m_tocdel)(size_t, const char *)                          = &P::tocdel;
void   (P::*m_tocclean)(size_t, const char *)                        = &P::tocclean;
void   (P::*m_tocprint)(size_t)                                      = &P::tocprint;
void   (P::*m_tocwrite)(size_t)                                      = &P::tocwrite;
size_t (P::*m_rd_toclen)(size_t)                                     = &P::rd_toclen;
void   (P::*m_rewind_toclen)(const size_t)                           = &P::rewind_toclen;

/* The singleton accessor and the global it returns: both live in init.cc,
 * which the port deletes, so both must be re-supplied. */
std::shared_ptr<P> (*m_shared_object)()                              = &P::shared_object;

/* The constants the shim depends on must also be Psi4's, not our copy's. A
 * wrong PSIO_PAGELEN would put every streamed address on the wrong page. */
static_assert(psi::PSIO_KEYLEN == 80,     "PSIO_KEYLEN changed in Psi4");
static_assert(psi::PSIO_PAGELEN == 65536, "PSIO_PAGELEN changed in Psi4");
static_assert(psi::PSIO_OPEN_NEW == 0,    "PSIO_OPEN_NEW changed in Psi4");
static_assert(psi::PSIO_OPEN_OLD == 1,    "PSIO_OPEN_OLD changed in Psi4");
static_assert(sizeof(psi::psio_address) == 2 * sizeof(size_t),
              "psio_address is no longer two size_t words");

}  // namespace

/* filecfg_kwd, getpid come from the kept libpsio files (filescfg.cc, getpid.cc)
 * and psi_file_prefix from test/stubs.cc; this test links the real psio archive,
 * so -- unlike libspill's standalone harness -- it stubs none of them. */

int main()
{
    /* Constructing one proves the constructor and destructor exist, which
     * nothing else here would. */
    { psi::PSIO probe; (void)probe; }

    const void *all[] = {(const void *)p_init, (const void *)p_open, (const void *)p_close,
                         (const void *)p_open_check, (const void *)p_get_address,
                         (const void *)p_get_global, (const void *)p_write, (const void *)p_read,
                         (const void *)p_write_entry, (const void *)p_read_entry,
                         (const void *)p_tocscan, (const void *)p_tocexists,
                         (const void *)p_tocprint, (const void *)p_tocwrite,
                         (const void *)p_rd_toclen, (const void *)p_error};

    touch(m_open); touch(m_close); touch(m_open_check); touch(m_exists);
    touch(m_rehash); touch(m_read); touch(m_write); touch(m_read_entry);
    touch(m_write_entry); touch(m_rw); touch(m_zero_disk); touch(m_zero_disk_s);
    touch(m_tocscan); touch(m_tocentry_exists); touch(m_tocdel); touch(m_tocclean);
    touch(m_tocprint); touch(m_tocwrite); touch(m_rd_toclen); touch(m_rewind_toclen);
    touch(m_shared_object);
    std::printf("  [PASS] 21 PSIO class methods match Psi4's own declarations\n");

    /* Linking is necessary but not sufficient -- issue #1 was a shim that
     * linked its own half of the API perfectly. Drive the class through a real
     * round trip, which is what consumers do. */
    {
        const char *td = std::getenv("TMPDIR");
        psi::psio_init();   /* builds both singletons, as Psi4 does at startup */
        psi::PSIOManager::shared_object()->set_default_path(td && *td ? td : "/tmp");

        auto psio = std::make_shared<psi::PSIO>();
        const std::string dir = psi::PSIOManager::shared_object()->get_file_path(99);
        double out[256], back[256];
        for (int k = 0; k < 256; k++) out[k] = 1.0 + 0.5 * k;

        psio->open(99, psi::PSIO_OPEN_NEW);
        if (!psio->open_check(99)) { std::printf("  [FAIL] open_check\n"); return 1; }

        /* Correct reads and writes are not enough. The second way this port can
         * be wrong -- and was -- is a shim that stores the data somewhere of its
         * own choosing. Psi4 requires the store to sit exactly where PSIOManager
         * says: psiclean() unlinks the paths in PSIOManager's ledger and nothing
         * else, PSI_SCRATCH reaches the I/O layer only through it, and file
         * retention is keyed on those same strings. A shim that ignores it reads
         * and writes perfectly while leaking every scratch file and silently
         * dropping PSI_SCRATCH.
         *
         * PSIOManager mirrors its ledger to psi.<pid>.clean on every change, so
         * read it back and stat what it names. */
        {
            const std::string mirror = "psi." + psi::psio_getpid() + ".clean";
            std::ifstream in(mirror.c_str());
            std::string line, registered;
            while (std::getline(in, line))
                if (line.size() > 3 && line.compare(line.size() - 3, 3, ".99") == 0) registered = line;

            if (registered.empty()) {
                std::printf("  [FAIL] unit 99 never reached PSIOManager (%s): "
                            "psiclean() would not clean it and retention is inert\n",
                            mirror.c_str());
                return 1;
            }
            if (registered.compare(0, dir.size(), dir) != 0) {
                std::printf("  [FAIL] scratch path %s is not under PSIOManager's directory %s: "
                            "PSI_SCRATCH is being ignored\n", registered.c_str(), dir.c_str());
                return 1;
            }
            struct stat st;
            if (stat(registered.c_str(), &st) != 0) {
                std::printf("  [FAIL] nothing on disk at the path PSIOManager recorded: %s\n",
                            registered.c_str());
                return 1;
            }
            if (!psio->exists(99)) { std::printf("  [FAIL] exists() on an open unit\n"); return 1; }
        }
        psio->write_entry(99, "Amplitudes", (char *)out, sizeof out);
        std::memset(back, 0, sizeof back);
        psio->read_entry(99, "Amplitudes", (char *)back, sizeof back);
        if (std::memcmp(out, back, sizeof out) != 0) {
            std::printf("  [FAIL] class round trip is not byte-exact\n");
            return 1;
        }
        if (psio->tocscan(99, "Amplitudes") == nullptr ||
            psio->tocscan(99, "Never Written") != nullptr) {
            std::printf("  [FAIL] tocscan through the class\n");
            return 1;
        }
        psio->zero_disk(99, "Zeroed", 8, 32);
        psio->read_entry(99, "Zeroed", (char *)back, 8 * 32 * sizeof(double) / 8);
        psio->close(99, 0);
        if (psio->exists(99)) {
            std::printf("  [FAIL] close(keep=0) left the store on disk\n");
            return 1;
        }
        std::printf("  [PASS] the store lives at PSIOManager's path and is registered for psiclean\n");

        /* Two instances must not share scratch: Psi4 makes several. */
        auto a = std::make_shared<psi::PSIO>();
        auto b = std::make_shared<psi::PSIO>();
        a->open(98, psi::PSIO_OPEN_NEW);
        b->open(97, psi::PSIO_OPEN_NEW);
        if (a->open_check(97) || b->open_check(98)) {
            std::printf("  [FAIL] PSIO instances share unit state\n");
            return 1;
        }
        a->close(98, 0);
        b->close(97, 0);
        std::printf("  [PASS] a PSIO instance round-trips, and instances are independent\n");
    }
    size_t i, n = sizeof all / sizeof all[0];
    for (i = 0; i < n; i++)
        if (!all[i]) { std::printf("  [FAIL] null entry point %zu\n", i); return 1; }
    std::printf("  [PASS] %zu entry points match Psi4's own declarations\n", n);
    std::printf("  [PASS] PSIO_KEYLEN, PSIO_PAGELEN and the open modes are unchanged\n");
    return 0;
}
