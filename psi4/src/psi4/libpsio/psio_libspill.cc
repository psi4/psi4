// SPDX-License-Identifier: BSD-3-Clause
#include "psio_libspill.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "libspill.h"
}

#ifdef PSIO_USE_PSI4_HEADERS
#include "psi4/libpsio/psio.hpp"
// psi_file_prefix is psi::psi_file_prefix with C++ linkage (core.cc defines it,
// the test's stubs.cc stands in for it); init.cc reached it through this header.
#include "psi4/psi4-dec.h"
#endif

namespace psi {

#ifndef PSIO_USE_PSI4_HEADERS
psio_address PSIO_ZERO = {0, 0};
#endif

namespace {

#ifndef PSIO_USE_PSI4_HEADERS
std::string g_dir;                       // in Psi4 this is PSIOManager's job
#endif
std::mutex  g_lk;

// One libspill store per (PSIO instance, unit).
struct Unit {
    ls_store *store = nullptr;
    // The path PSIOManager was told about, kept so close() can report the same
    // string back to it. The old open.cc kept it in psio_ud::vol.path.
    std::string path;
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

// Where a unit's store goes.
//
// In Psi4 this is not libspill's decision and not the shim's: PSIOManager owns
// it. It computes one path per unit -- PSI_SCRATCH or a per-file override, plus
// the pid/namespace decoration -- hands that path to whoever opens the file,
// and psiclean() later unlinks exactly the paths it was told about. So the shim
// must open the path PSIOManager names, byte for byte, which is what
// ls_opts.exact_name is for: it suppresses libspill's own
// <dir>/<name>.libspill layout. Without it the scratch files land in TMPDIR
// under names PSIOManager never sees -- PSI_SCRATCH ignored, psiclean()
// cleaning nothing, retention inert.
struct StorePath {
    std::string dir;   // no trailing separator
    std::string name;  // basename; taken verbatim when exact is set
    bool exact = false;
    std::string full() const { return dir.empty() ? name : dir + "/" + name; }
};

ls_opts opts_for(const StorePath &p) {
    ls_opts o;
    ls_opts_default(&o);
    if (!p.dir.empty()) o.dir = p.dir.c_str();   // borrows p; p outlives o
    o.exact_name = p.exact ? 1 : 0;
    o.log = log_cb;
    return o;
}

#ifdef PSIO_USE_PSI4_HEADERS
// Split PSIO::get_unit_filename()'s "<path><name>.<unit>" into the directory
// and basename libspill wants. Psi4 composes it (get_filename.cc, and the
// helper #3515 factored out of open.cc); the shim only takes it apart.
StorePath split_path(const std::string &full) {
    StorePath p;
    const std::string::size_type slash = full.find_last_of('/');
    if (slash == std::string::npos) {
        p.name = full;
    } else {
        p.dir = full.substr(0, slash);
        p.name = full.substr(slash + 1);
    }
    p.exact = true;
    return p;
}
#else
// Standalone: no PSIOManager to defer to, so let libspill name the file.
StorePath standalone_path(size_t unit) {
    StorePath p;
    p.dir = g_dir;
    p.name = "psi." + std::to_string(unit);
    return p;
}
#endif

// ---------------------------------------------------------------------------
// The logic, independent of who is calling it. The PSIO class methods and the
// free functions are both thin wrappers around these, so the two APIs cannot
// drift apart -- which is how issue #1 happened in the first place.
// Every impl_* expects g_lk to be held, except impl_exists which takes it.

void impl_open(Units &u, size_t unit, int status, const StorePath &p) {
    if (unit >= static_cast<size_t>(PSIO_MAXUNIT)) psio_error(unit, PSIO_ERROR_MAXUNIT);
    if (find_unit(u, unit)) return;                      // already open

    ls_opts o = opts_for(p);

    // PSIO_OPEN_NEW must not inherit a previous run's contents; libspill
    // reopens a kept store when one is there, so remove it first.
    if (status == PSIO_OPEN_NEW) {
        int err = 0;
        ls_store *old = ls_open(p.name.c_str(), &o, &err);
        if (old) ls_close(old, 0);                       // keep=0 unlinks
    }
    int err = 0;
    ls_store *s = ls_open(p.name.c_str(), &o, &err);
    if (!s) psio_error(unit, PSIO_ERROR_OPEN);
    Unit &nu = u[unit];
    nu.store = s;
    nu.path = p.full();
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

// Not part of Psi4's API. In Psi4 the scratch directory is PSIOManager's, so
// this sets it there rather than keeping a second, divergent copy; g_dir is
// only the standalone build's fallback.
void psio_set_scratch_dir(const char *dir) {
#ifdef PSIO_USE_PSI4_HEADERS
    psio_init();
    PSIOManager::shared_object()->set_default_path(dir ? dir : "");
#else
    g_dir = dir ? dir : "";
#endif
}

// ---------------------------------------------------------------- addresses

psio_address psio_get_address(psio_address start, size_t shift) {
    return to_address(to_bytes(start) + shift);
}

psio_address psio_get_global_address(psio_address entry_start, psio_address rel_address) {
    return to_address(to_bytes(entry_start) + to_bytes(rel_address));
}

#ifndef PSIO_USE_PSI4_HEADERS
// Psi4 keeps its own error.cc, whose messages are far better than this; the
// standalone build has no such thing, so it gets a minimal stand-in.
void psio_error(size_t unit, size_t errval, std::string prev_msg) {
    std::fprintf(stderr, "PSIO_ERROR: unit %zu, errval %zu %s\n", unit, errval,
                 prev_msg.c_str());
    throw std::runtime_error("PSIO error " + std::to_string(errval) + " on unit " +
                             std::to_string(unit) + " " + prev_msg);
}
#endif

// ------------------------------------------------------------- the PSIO class
//
// This is the interface Psi4 uses: roughly 2400 call sites go through a
// shared_ptr<PSIO>, against a few hundred through the free functions below.
// Issue #1 was that an earlier version of this shim implemented only the free
// layer, so the port did not link.

#ifdef PSIO_USE_PSI4_HEADERS

// All of these live in init.cc, which this port deletes. _default_psio_manager_
// is not optional decoration: psio.hpp declares it extern, filemanager.cc is
// the only other file that names it, and every scratch path in the shim now
// comes from it.
psio_address PSIO_ZERO = {0, 0};
std::string PSIO::default_namespace_;
std::shared_ptr<PSIO> _default_psio_lib_;
std::shared_ptr<PSIOManager> _default_psio_manager_;

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
    files_keywords_.clear();          // done.cc did this
}

void PSIO::open(size_t unit, int status) {
    // get_unit_filename() is Psi4's composition of PSIOManager's directory and
    // the pid/namespace-decorated name; it touches no shim state, so it is
    // computed outside the lock.
    const StorePath p = split_path(get_unit_filename(unit));
    {
        std::lock_guard<std::mutex> g(g_lk);
        impl_open(units_of(this), unit, status, p);
    }
    // Registering the file is what makes psiclean() and per-file retention
    // work at all -- PSIOManager only ever unlinks paths in this ledger. The
    // old open.cc registered here too.
    PSIOManager::shared_object()->open_file(p.full(), unit);
}

void PSIO::close(size_t unit, int keep) {
    std::string path;
    {
        std::lock_guard<std::mutex> g(g_lk);
        Units &u = units_of(this);
        if (Unit *cu = find_unit(u, unit)) path = cu->path;
        impl_close(u, unit, keep);       // ls_close(keep=0) unlinks the file
    }
    if (!path.empty()) PSIOManager::shared_object()->close_file(path, unit, keep != 0);
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
    // The old exists() opened the path O_RDWR and closed it again. Asking
    // libspill is the same question without the shim having to reconstruct a
    // layout it no longer owns -- and without creating the store, which is why
    // this is ls_store_exists rather than ls_open.
    const StorePath p = split_path(get_unit_filename(unit));
    ls_opts o = opts_for(p);
    int found = 0;
    return ls_store_exists(p.name.c_str(), &o, &found) == LS_OK && found;
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
    if (_default_psio_manager_.get() == 0) _default_psio_manager_ = std::make_shared<PSIOManager>();
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
    impl_open(default_units(), unit, status, standalone_path(unit));
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
