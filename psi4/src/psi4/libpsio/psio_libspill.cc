// SPDX-License-Identifier: BSD-3-Clause
#include "psio_libspill.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

extern "C" {
#include "libspill.h"
}

#ifdef PSIO_USE_PSI4_HEADERS
#include "psi4/libpsio/psio.hpp"
extern "C" {
extern PSI_API char *psi_file_prefix;
}
#endif

namespace psi {

#ifndef PSIO_USE_PSI4_HEADERS
psio_address PSIO_ZERO = {0, 0};
#endif

namespace {

std::string g_dir;                       // PSIOManager's job in Psi4
std::mutex  g_lk;

// One libspill store per (PSIO instance, unit).
struct Unit {
    ls_store *store = nullptr;
    // Entries handed back by tocscan. Psi4 never dereferences these, but a
    // function returning a pointer should return a valid one; per-unit storage
    // gives them the lifetime the old linked-list TOC gave.
    std::map<std::string, psio_tocentry> toc;
    std::vector<std::string> order;       // insertion order, for tocclean
};

// State lives per PSIO INSTANCE, not per process. Psi4 creates several --
// std::make_shared<PSIO>() appears in dfmp2, mcscf, DiskJK and fnocc among
// others -- so a file-global unit table would have them share scratch. It is
// keyed on `this` rather than held in a member because psio.hpp is unchanged
// by this port, which is the point of it.
// Keyed on void* rather than PSIO* so the same table serves the standalone
// build, where there is no PSIO class at all and the free functions use the
// nullptr key.
using Units = std::map<size_t, Unit>;
std::map<const void *, Units> g_state;

Units &units_of(const void *self) { return g_state[self]; }

Unit *find_unit(Units &u, size_t unit) {
    auto it = u.find(unit);
    return it == u.end() ? nullptr : &it->second;
}

inline uint64_t to_bytes(psio_address a) {
    return static_cast<uint64_t>(a.page) * static_cast<uint64_t>(PSIO_PAGELEN) + a.offset;
}

inline psio_address to_address(uint64_t b) {
    psio_address a;
    a.page = static_cast<size_t>(b / PSIO_PAGELEN);
    a.offset = static_cast<size_t>(b % PSIO_PAGELEN);
    return a;
}

int psio_err_of(int rc, int on_read) {
    switch (rc) {
        case LS_ERR_NOKEY: return PSIO_ERROR_NOTOCENT;
        case LS_ERR_RANGE: return PSIO_ERROR_BLKEND;
        case LS_ERR_INVAL: return PSIO_ERROR_KEYLEN;
        default:           return on_read ? PSIO_ERROR_READ : PSIO_ERROR_WRITE;
    }
}

void log_cb(int err, const char *key, uint64_t off, size_t nbytes, const char *msg, void *) {
    char buf[LS_ERRBUF_MIN];
    std::fprintf(stderr, "PSIO/libspill: %s (key %s, offset %llu, %zu bytes): %s\n",
                 msg ? msg : "operation", key ? key : "-",
                 static_cast<unsigned long long>(off), nbytes,
                 ls_strerror(err, buf, sizeof buf));
}

std::string name_of(size_t unit) { return "psi." + std::to_string(unit); }

#ifdef PSIO_USE_PSI4_HEADERS
// libspill owns the file layout, but PSIO::exists has to answer "is there
// something on disk for this unit" without opening it -- opening would create
// it. Reproducing the one path rule libspill applies is the least-bad way.
std::string path_of(size_t unit) {
    std::string dir = g_dir;
    if (dir.empty()) {
        const char *t = std::getenv("TMPDIR");
        dir = (t && *t) ? t : "/tmp";
    }
    return dir + "/" + name_of(unit) + ".libspill";
}
#endif

// ---------------------------------------------------------------------------
// The logic, independent of who is calling it. The PSIO class methods and the
// free functions are both thin wrappers around these, so the two APIs cannot
// drift apart -- which is how issue #1 happened in the first place.
// Every impl_* expects g_lk to be held, except impl_exists which takes it.

void impl_open(Units &u, size_t unit, int status) {
    if (unit >= static_cast<size_t>(PSIO_MAXUNIT)) psio_error(unit, PSIO_ERROR_MAXUNIT);
    if (find_unit(u, unit)) return;                      // already open

    ls_opts o;
    ls_opts_default(&o);
    if (!g_dir.empty()) o.dir = g_dir.c_str();
    o.log = log_cb;

    // PSIO_OPEN_NEW must not inherit a previous run's contents; libspill
    // reopens a kept store when one is there, so remove it first.
    if (status == PSIO_OPEN_NEW) {
        int err = 0;
        ls_store *old = ls_open(name_of(unit).c_str(), &o, &err);
        if (old) ls_close(old, 0);                       // keep=0 unlinks
    }
    int err = 0;
    ls_store *s = ls_open(name_of(unit).c_str(), &o, &err);
    if (!s) psio_error(unit, PSIO_ERROR_OPEN);
    u[unit].store = s;
}

void impl_close(Units &u, size_t unit, int keep) {
    Unit *cu = find_unit(u, unit);
    if (!cu) psio_error(unit, PSIO_ERROR_UNOPENED);
    int rc = ls_close(cu->store, keep);
    u.erase(unit);
    if (rc != LS_OK) psio_error(unit, PSIO_ERROR_CLOSE);
}

void impl_rw(Units &u, size_t unit, const char *key, char *buffer, size_t size,
             psio_address sadd, psio_address *eadd, bool is_write) {
    Unit *cu = find_unit(u, unit);
    if (!cu) psio_error(unit, PSIO_ERROR_UNOPENED);
    if (!key || std::strlen(key) >= static_cast<size_t>(PSIO_KEYLEN))
        psio_error(unit, PSIO_ERROR_KEYLEN);

    const uint64_t off = to_bytes(sadd);
    const int rc = is_write ? ls_write(cu->store, key, off, size, buffer)
                            : ls_read (cu->store, key, off, size, buffer);
    if (rc != LS_OK) psio_error(unit, psio_err_of(rc, !is_write));
    if (is_write && cu->toc.find(key) == cu->toc.end()) cu->order.push_back(key);
    if (eadd) *eadd = to_address(off + size);
}

void impl_zero_disk(Units &u, size_t unit, const char *key, size_t rows, size_t cols) {
    Unit *cu = find_unit(u, unit);
    if (!cu) psio_error(unit, PSIO_ERROR_UNOPENED);
    const int rc = ls_reserve(cu->store, key,
                              static_cast<uint64_t>(rows) * cols * sizeof(double));
    if (rc != LS_OK) psio_error(unit, psio_err_of(rc, 0));
    if (cu->toc.find(key) == cu->toc.end()) cu->order.push_back(key);
}

psio_tocentry *impl_tocscan(Units &u, size_t unit, const char *key) {
    Unit *cu = find_unit(u, unit);
    if (!cu || !key) return nullptr;
    uint64_t nbytes = 0;
    if (ls_size(cu->store, key, &nbytes) != LS_OK) return nullptr;

    psio_tocentry &e = cu->toc[key];
    std::snprintf(e.key, PSIO_KEYLEN, "%s", key);
    e.sadd = PSIO_ZERO;
    e.eadd = to_address(nbytes);
    e.next = nullptr;
    e.last = nullptr;
    return &e;
}

bool impl_tocentry_exists(Units &u, size_t unit, const char *key) {
    Unit *cu = find_unit(u, unit);
    if (!cu || !key) return false;
    int found = 0;
    return ls_exists(cu->store, key, &found) == LS_OK && found;
}

bool impl_tocdel(Units &u, size_t unit, const char *key) {
    Unit *cu = find_unit(u, unit);
    if (!cu || !key) return false;
    cu->toc.erase(key);
    for (auto it = cu->order.begin(); it != cu->order.end(); ++it)
        if (*it == key) { cu->order.erase(it); break; }
    return ls_erase(cu->store, key) == LS_OK;
}

void impl_tocprint(Units &u, size_t unit) {
    Unit *cu = find_unit(u, unit);
    if (!cu) return;
    char **keys = nullptr;
    size_t n = 0;
    if (ls_keys(cu->store, &keys, &n) != LS_OK) return;
    for (size_t i = 0; i < n; i++) {
        uint64_t sz = 0;
        ls_size(cu->store, keys[i], &sz);
        std::printf("%-80s %llu\n", keys[i], static_cast<unsigned long long>(sz));
    }
    ls_keys_free(keys, n);
}

size_t impl_toclen(Units &u, size_t unit) {
    Unit *cu = find_unit(u, unit);
    if (!cu) return 0;
    char **keys = nullptr;
    size_t n = 0;
    if (ls_keys(cu->store, &keys, &n) != LS_OK) return 0;
    ls_keys_free(keys, n);
    return n;
}

}  // namespace

void psio_set_scratch_dir(const char *dir) { g_dir = dir ? dir : ""; }

// ---------------------------------------------------------------- addresses

psio_address psio_get_address(psio_address start, size_t shift) {
    return to_address(to_bytes(start) + shift);
}

psio_address psio_get_global_address(psio_address entry_start, psio_address rel_address) {
    return to_address(to_bytes(entry_start) + to_bytes(rel_address));
}

void psio_error(size_t unit, size_t errval, std::string prev_msg) {
    std::fprintf(stderr, "PSIO_ERROR: unit %zu, errval %zu %s\n", unit, errval,
                 prev_msg.c_str());
    throw std::runtime_error("PSIO error " + std::to_string(errval) + " on unit " +
                             std::to_string(unit) + " " + prev_msg);
}

// ------------------------------------------------------------- the PSIO class
//
// This is the interface Psi4 uses: roughly 2400 call sites go through a
// shared_ptr<PSIO>, against a few hundred through the free functions below.
// Issue #1 was that an earlier version of this shim implemented only the free
// layer, so the port did not link.

#ifdef PSIO_USE_PSI4_HEADERS

// All three live in init.cc, which this port deletes.
psio_address PSIO_ZERO = {0, 0};
std::string PSIO::default_namespace_;
std::shared_ptr<PSIO> _default_psio_lib_;

std::shared_ptr<PSIO> PSIO::shared_object() { return _default_psio_lib_; }

// Reproduces init.cc's constructor exactly, because the files this port keeps
// depend on what it sets: get_filename.cc and filemanager.cc read pid_, and
// filescfg.cc reads files_keywords_, which filecfg_kwd populates.
PSIO::PSIO() {
    psio_unit = (psio_ud *)malloc(sizeof(psio_ud) * PSIO_MAXUNIT);
    if (psio_unit == nullptr) throw std::runtime_error("Error in PSIO_INIT()!\n");
    state_ = 1;
    for (int i = 0; i < PSIO_MAXUNIT; i++) {
        psio_unit[i].vol.path = nullptr;
        psio_unit[i].vol.stream = -1;
        psio_unit[i].toclen = 0;
        psio_unit[i].toc = nullptr;
    }
    filecfg_kwd("DEFAULT", "NAME", -1, psi_file_prefix);
    pid_ = getpid();
    std::lock_guard<std::mutex> g(g_lk);
    g_state[this];                       // give this instance its own table
}

PSIO::~PSIO() {
    {
        std::lock_guard<std::mutex> g(g_lk);
        auto it = g_state.find(this);
        if (it != g_state.end()) {
            for (auto &kv : it->second)
                if (kv.second.store) ls_close(kv.second.store, 0);
            g_state.erase(it);
        }
    }
    if (psio_unit) free(psio_unit);
    psio_unit = nullptr;
    state_ = 0;
}

void PSIO::open(size_t unit, int status) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_open(units_of(this), unit, status);
}

void PSIO::close(size_t unit, int keep) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_close(units_of(this), unit, keep);
}

int PSIO::open_check(size_t unit) const {
    std::lock_guard<std::mutex> g(g_lk);
    return find_unit(units_of(this), unit) != nullptr ? 1 : 0;
}

bool PSIO::exists(size_t unit) const {
    {
        std::lock_guard<std::mutex> g(g_lk);
        if (find_unit(units_of(this), unit)) return true;
    }
    struct stat st;
    return stat(path_of(unit).c_str(), &st) == 0 && st.st_size > 0;
}

// libspill's table of contents is live and in memory; there is nothing to
// re-read. Kept so call sites need not change.
void PSIO::rehash(size_t) {}

void PSIO::read(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
                psio_address *eadd) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_rw(units_of(this), unit, key, buffer, size, sadd, eadd, false);
}

void PSIO::write(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
                 psio_address *eadd) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_rw(units_of(this), unit, key, buffer, size, sadd, eadd, true);
}

void PSIO::read_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    read(unit, key, buffer, size, PSIO_ZERO, &end);
}

void PSIO::write_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    write(unit, key, buffer, size, PSIO_ZERO, &end);
}

// Psi4 writes rows*cols zeroed doubles one row at a time -- `rows` separate
// calls through the whole stack. ls_reserve is one call, and on Linux one
// fallocate. 59 call sites.
void PSIO::zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_zero_disk(units_of(this), unit, key, rows, cols);
}

void PSIO::zero_disk(size_t unit, const std::string &key, size_t rows, size_t cols) {
    zero_disk(unit, key.c_str(), rows, cols);
}

// Raw positional access, bypassing the table of contents: Psi4 uses it to read
// and write its own on-disk TOC. libspill owns that now, so this addresses a
// reserved key rather than raw file offsets.
void PSIO::rw(size_t unit, char *buffer, psio_address address, size_t size, int wrt) {
    if (wrt)
        write(unit, "__psio_raw__", buffer, size, address, nullptr);
    else
        read(unit, "__psio_raw__", buffer, size, address, nullptr);
}

psio_tocentry *PSIO::tocscan(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocscan(units_of(this), unit, key);
}

bool PSIO::tocentry_exists(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocentry_exists(units_of(this), unit, key);
}

bool PSIO::tocdel(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocdel(units_of(this), unit, key);
}

// Psi4 drops `key` and everything written after it. libspill's table has no
// order of its own, so the shim tracks insertion order to reproduce this.
void PSIO::tocclean(size_t unit, const char *key) {
    std::vector<std::string> doomed;
    {
        std::lock_guard<std::mutex> g(g_lk);
        Unit *cu = find_unit(units_of(this), unit);
        if (!cu) return;
        bool hit = (key == nullptr);
        for (const auto &k : cu->order) {
            if (!hit && k == key) hit = true;
            if (hit) doomed.push_back(k);
        }
    }
    for (const auto &k : doomed) tocdel(unit, k.c_str());
}

void PSIO::tocprint(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_tocprint(units_of(this), unit);
}

// The on-disk table of contents is libspill's now. tocwrite has one consumer
// outside libpsio and rd_toclen has none, which is why that ownership transfer
// costs Psi4 nothing.
void PSIO::tocwrite(size_t) {}
void PSIO::rewind_toclen(const size_t) {}

size_t PSIO::rd_toclen(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_toclen(units_of(this), unit);
}

// Private helpers, declared in psio.hpp and used only from within libpsio.
psio_tocentry *PSIO::toclast(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    Unit *cu = find_unit(units_of(this), unit);
    if (!cu || cu->order.empty()) return nullptr;
    auto it = cu->toc.find(cu->order.back());
    return it == cu->toc.end() ? nullptr : &it->second;
}

// #3515 made toclen() const; it cannot delegate to the non-const rd_toclen, so
// it calls the shared impl directly (the same read-only key count).
size_t PSIO::toclen(size_t unit) const {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_toclen(units_of(this), unit);
}
void PSIO::wt_toclen(size_t, size_t) {}
void PSIO::tocread(size_t) {}

// ------------------------------------------------------- the free functions
//
// Psi4's own structure: these delegate to the default instance. They are a
// secondary API -- a few hundred call sites against the class's ~2400 -- but
// they are part of psio.h and consumers do use them.

int psio_init() {
    if (_default_psio_lib_.get() == 0) _default_psio_lib_ = std::make_shared<PSIO>();
    return 1;
}

int psio_done() {
    _default_psio_lib_.reset();
    return 1;
}

int psio_open(size_t unit, int status) {
    psio_init();
    _default_psio_lib_->open(unit, status);
    return 1;
}
int psio_close(size_t unit, int keep) { _default_psio_lib_->close(unit, keep); return 1; }
int psio_open_check(size_t unit) { return _default_psio_lib_ ? _default_psio_lib_->open_check(unit) : 0; }

int psio_write(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
               psio_address *eadd) {
    _default_psio_lib_->write(unit, key, buffer, size, sadd, eadd);
    return 1;
}
int psio_read(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
              psio_address *eadd) {
    _default_psio_lib_->read(unit, key, buffer, size, sadd, eadd);
    return 1;
}
int psio_write_entry(size_t unit, const char *key, char *buffer, size_t size) {
    _default_psio_lib_->write_entry(unit, key, buffer, size);
    return 1;
}
int psio_read_entry(size_t unit, const char *key, char *buffer, size_t size) {
    _default_psio_lib_->read_entry(unit, key, buffer, size);
    return 1;
}
psio_tocentry *psio_tocscan(size_t unit, const char *key) {
    return _default_psio_lib_ ? _default_psio_lib_->tocscan(unit, key) : nullptr;
}
bool psio_tocentry_exists(size_t unit, const char *key) {
    return _default_psio_lib_ ? _default_psio_lib_->tocentry_exists(unit, key) : false;
}
bool psio_tocdel(size_t unit, const char *key) {
    return _default_psio_lib_ ? _default_psio_lib_->tocdel(unit, key) : false;
}
void psio_tocprint(size_t unit) { if (_default_psio_lib_) _default_psio_lib_->tocprint(unit); }
int  psio_tocwrite(size_t) { return 1; }
size_t psio_rd_toclen(size_t unit) { return _default_psio_lib_ ? _default_psio_lib_->rd_toclen(unit) : 0; }
void psio_zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    _default_psio_lib_->zero_disk(unit, key, rows, cols);
}

#else   // standalone: no PSIO class to delegate to, so drive the impl layer
        // directly through one process-wide table.

namespace { Units &default_units() { return units_of(nullptr); } }

int psio_init() { return 1; }
int psio_done() {
    std::lock_guard<std::mutex> g(g_lk);
    for (auto &kv : default_units()) if (kv.second.store) ls_close(kv.second.store, 0);
    g_state.erase(static_cast<const void *>(nullptr));
    return 1;
}
int psio_open(size_t unit, int status) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_open(default_units(), unit, status);
    return 1;
}
int psio_close(size_t unit, int keep) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_close(default_units(), unit, keep);
    return 1;
}
int psio_open_check(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    return find_unit(default_units(), unit) != nullptr ? 1 : 0;
}
int psio_write(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
               psio_address *eadd) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_rw(default_units(), unit, key, buffer, size, sadd, eadd, true);
    return 1;
}
int psio_read(size_t unit, const char *key, char *buffer, size_t size, psio_address sadd,
              psio_address *eadd) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_rw(default_units(), unit, key, buffer, size, sadd, eadd, false);
    return 1;
}
int psio_write_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
}
int psio_read_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    return psio_read(unit, key, buffer, size, PSIO_ZERO, &end);
}
psio_tocentry *psio_tocscan(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocscan(default_units(), unit, key);
}
bool psio_tocentry_exists(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocentry_exists(default_units(), unit, key);
}
bool psio_tocdel(size_t unit, const char *key) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_tocdel(default_units(), unit, key);
}
void psio_tocprint(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_tocprint(default_units(), unit);
}
int psio_tocwrite(size_t) { return 1; }
size_t psio_rd_toclen(size_t unit) {
    std::lock_guard<std::mutex> g(g_lk);
    return impl_toclen(default_units(), unit);
}
void psio_zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    std::lock_guard<std::mutex> g(g_lk);
    impl_zero_disk(default_units(), unit, key, rows, cols);
}

#endif  // PSIO_USE_PSI4_HEADERS

}  // namespace psi
