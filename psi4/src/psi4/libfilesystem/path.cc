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

//
// Created by Justin Turney on 9/1/16.
//

#include "path.h"

#include <stdexcept>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <system_error>

#ifdef _MSC_VER
#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <stdlib.h>
#define S_IRUSR 0
#define S_IWUSR 0
#define S_IXUSR 0
#define S_ISREG(m) (((m)&S_IFMT) == S_IFREG)
#define S_ISDIR(m) (((m)&S_IFMT) == S_IFDIR)
#define SYSTEM_GETCWD ::_getcwd
#define SYSTEM_MKDIR(D, P) ::_mkdir((D))
#define PATH_MAX _MAX_PATH
#define SYSTEM_REALPATH(N, R) ::_fullpath((R), (N), _MAX_PATH)
static int SYSTEM_TRUNCATE(const char *path, off_t length) {
    int descriptor = ::_open(path, _O_BINARY | _O_WRONLY);
    ::_chsize(descriptor, length);
    return ::_close(descriptor);
}
#define PATH_SEPARATOR "\\"
#else
#include <unistd.h>
#define SYSTEM_GETCWD ::getcwd
#define SYSTEM_MKDIR ::mkdir
#define SYSTEM_REALPATH ::realpath
#define SYSTEM_TRUNCATE ::truncate
#define PATH_SEPARATOR "/"
#endif
#include <sys/stat.h>
#include <climits>

namespace psi {
namespace filesystem {

bool create_directory(const path &p) { return SYSTEM_MKDIR(p.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR) == 0; }

path path::make_absolute() const {
    // TODO: Handle ~ and environment variables (e.g. $HOME)
    int path_max;
#ifdef PATH_MAX
    path_max = PATH_MAX;
#else
    path_max = pathconf(path, _PC_PATH_MAX);
    if (path_max <= 0) path_max = 4096;
#endif

    auto *temp = new char[path_max];
    if (SYSTEM_REALPATH(str().c_str(), temp) == nullptr) {
        // Ignore errors relating to a file or directory component not existing
        if (errno != (int)std::errc::no_such_file_or_directory && errno != (int)std::errc::not_a_directory) {
            throw std::runtime_error("path::make_absolute: " + std::string(strerror(errno)));
        }
    }
    path ptemp(temp);
    delete[] temp;
    return ptemp;
}

bool path::exists() const {
    struct stat sb;
    return stat(str().c_str(), &sb) == 0;
}

void path::set(const std::string &str) {
    path_ = tokenize(str, PATH_SEPARATOR);
#ifdef _MSC_VER
    absolute_ = !str.empty() && str[1] == ':';
#else
    absolute_ = !str.empty() && str[0] == PATH_SEPARATOR[0];
#endif
}

std::string path::str() const {
    std::ostringstream ss;

#ifndef _MSC_VER
    if (absolute_) ss << PATH_SEPARATOR;
#endif

    for (size_t i = 0; i < path_.size(); ++i) {
        ss << path_[i];
        if (i + 1 < path_.size()) ss << PATH_SEPARATOR;
    }

    return ss.str();
}

bool path::is_directory() const {
    struct stat sb;
    if (stat(str().c_str(), &sb)) return false;
    return S_ISDIR(sb.st_mode);
}

bool path::is_file() const {
    struct stat sb;
    if (stat(str().c_str(), &sb)) return false;
    return S_ISREG(sb.st_mode);
}

std::string path::stem() const {
    std::string path = filename();
    return path.substr(0, path.find_last_of("."));
}

std::string path::filename() const {
    if (empty()) return "";
    const std::string &last = path_[path_.size() - 1];
    return last;
}

std::string path::extension() const {
    const std::string &name = filename();
    size_t pos = name.find_last_of(".");
    if (pos == std::string::npos) return "";
    return name.substr(pos + 1);
}

path path::parent_path() const {
    path result;
    result.absolute_ = absolute_;

    if (path_.empty()) {
        if (!absolute_) result.path_.push_back("..");
    } else {
        size_t until = path_.size() - 1;
        for (size_t i = 0; i < until; i++) result.path_.push_back(path_[i]);
    }
    return result;
}

path path::operator/(const path &other) const {
    if (other.absolute_) throw std::runtime_error("path::operator/(): expected a relative path");

    path result(*this);

    for (size_t i = 0; i < other.path_.size(); i++) result.path_.push_back(other.path_[i]);

    return result;
}

path &path::operator=(const path &path) {
    path_ = path.path_;
    absolute_ = path.absolute_;
    return *this;
}

path &path::operator=(path &&path) {
    if (this != &path) {
        path_ = std::move(path.path_);
        absolute_ = path.absolute_;
    }
    return *this;
}

bool path::remove_file() { return std::remove(str().c_str()) == 0; }

bool path::resize_file(size_t target_length) { return SYSTEM_TRUNCATE(str().c_str(), (off_t)target_length) == 0; }

path path::getcwd() {
    char temp[PATH_MAX];
    if (SYSTEM_GETCWD(temp, PATH_MAX) == nullptr)
        throw std::runtime_error("path::getcwd(): " + std::string(strerror(errno)));
    return path(temp);
}

}  // filesystem
}  // psi
