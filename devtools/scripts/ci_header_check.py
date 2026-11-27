#!/usr/bin/env python3
"""Check that psi4 headers (and optionally sources) compile on their own.

Unity builds hide missing includes: a header (or .cc) that forgets an include still compiles
when something earlier in the same translation unit provided it. This parses each header
standalone (one generated `#include "the/header.h"` TU each) with -fsyntax-only, using the
compiler and flags CMake recorded in compile_commands.json for a source in the same directory.
Unity build directories are fine: CMake's Unity/unity_*.cxx files are expanded to their members.

Headers are skipped when their directory has no compiled sources in this build, or when the
same-named .cc exists but isn't compiled (optional dependency not enabled, e.g. siminteri.h).

  ci_header_check.py BUILD_DIR SRC_ROOT [-j N] [--sources] [--exclude REGEX] [--times FILE]

  BUILD_DIR  dir holding compile_commands.json (for psi4's superbuild: <build>/psi4-core-prefix/src/psi4-core-build)
  SRC_ROOT   headers under here are checked (e.g. psi4/src/psi4)
  --sources  also parse every .cc on its own (several times costlier; about the parsing
             half of a full non-unity build)
"""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


_PCH_FLAGS = ("-Winvalid-pch", "-fpch-instantiate-templates", "-fno-pch-timestamp", "-fpch-validate-input-files-content")


def strip_pch(flags):
    """Drop precompiled-header flags (CMake's `-include .../cmake_pch.hxx`; clang's `-Xclang -include-pch -Xclang x.pch`),
    which would otherwise supply headers to every TU and hide missing includes. Returns (flags, stripped?)."""
    units, i = [], 0
    while i < len(flags):  # keep `-Xclang X` pairs together
        if flags[i] == "-Xclang" and i + 1 < len(flags):
            units.append((flags[i:i + 2], flags[i + 1]))
            i += 2
        else:
            units.append(([flags[i]], flags[i]))
            i += 1
    out, stripped, i = [], False, 0
    while i < len(units):
        toks, val = units[i]
        if val in ("-include", "-include-pch") and i + 1 < len(units) and re.search(r"cmake_pch|\.pch$|\.gch$", units[i + 1][1]):
            stripped, i = True, i + 2
            continue
        if re.match(r"-include.*cmake_pch", val) or val in _PCH_FLAGS:
            stripped, i = True, i + 1
            continue
        out += toks
        i += 1
    return out, stripped


def load_commands(build_dir):
    """[(source, compiler+flags without -c/-o/source/PCH, cwd)], unity TUs expanded to their member sources."""
    with open(os.path.join(build_dir, "compile_commands.json")) as fp:
        entries = json.load(fp)
    jobs, pch = [], False
    for e in entries:
        args = e["arguments"] if "arguments" in e else shlex.split(e["command"])
        src = e["file"]
        flags, skip = [], False
        for a in args:
            if skip:
                skip = False
                continue
            if a in ("-o", "-MF", "-MT", "-MQ"):
                skip = True
                continue
            if a in ("-c", "-MD", "-MMD") or a == src or (not a.startswith("-") and os.path.basename(a) == os.path.basename(src)):
                continue
            flags.append(a)
        if not src.endswith((".cc", ".cxx", ".cpp")):
            continue
        flags, stripped = strip_pch(flags)
        pch |= stripped
        if "/Unity/" in src:
            with open(src) as fp:
                members = re.findall(r'#include\s+"([^"]+)"', fp.read())
            jobs += [(os.path.normpath(os.path.join(os.path.dirname(src), m)), flags, e["directory"]) for m in members]
        else:
            jobs.append((os.path.normpath(src), flags, e["directory"]))
    if pch:
        print("note: build uses precompiled headers; PCH flags stripped so they can't supply missing includes")
    return jobs


def header_jobs(jobs, src_root, tmpdir, exclude):
    # A header with a same-named .cc uses that source's target flags (e.g. libmints/matrix.h gets mints's,
    #   not l2intf's); any other header uses the target that compiles the most sources in its directory.
    by_stem, dir_counts = {}, {}
    for src, flags, cwd in jobs:
        by_stem[os.path.splitext(src)[0]] = (flags, cwd)
        counts = dir_counts.setdefault(os.path.dirname(src), {})
        key = (tuple(flags), cwd)
        counts[key] = counts.get(key, 0) + 1
    by_dir = {d: (list(max(c, key=c.get)[0]), max(c, key=c.get)[1]) for d, c in dir_counts.items()}
    compiled = set(by_stem)
    out, skipped = [], []
    for root, _, files in os.walk(src_root):
        for f in sorted(files):
            if not f.endswith((".h", ".hpp")):
                continue
            hdr = os.path.join(root, f)
            stem = os.path.splitext(hdr)[0]
            sibling_uncompiled = any(os.path.exists(stem + ext) for ext in (".cc", ".cpp")) and stem not in compiled
            if root not in by_dir or sibling_uncompiled or (exclude and re.search(exclude, hdr)):
                skipped.append(hdr)
                continue
            flags, cwd = by_stem.get(stem, by_dir[root])
            tu = os.path.join(tmpdir, hdr.strip("/").replace("/", "__") + ".cc")
            with open(tu, "w") as fp:
                fp.write(f'#include "{hdr}"\n')
            out.append((tu, flags, cwd, hdr))
    return out, skipped


def run_one(src, flags, cwd, label):
    t0 = time.monotonic()
    p = subprocess.run(flags + ["-fsyntax-only", src], cwd=cwd, capture_output=True, text=True)
    return label, p.returncode, time.monotonic() - t0, p.stderr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build_dir")
    ap.add_argument("src_root")
    ap.add_argument("-j", type=int, default=2)
    ap.add_argument("--sources", action="store_true", help="also parse each .cc on its own")
    ap.add_argument("--exclude", help="regex of header/source paths to skip")
    ap.add_argument("--times", help="write per-TU seconds to this file")
    args = ap.parse_args()

    jobs = load_commands(args.build_dir)
    src_root = os.path.abspath(args.src_root)
    tmp = tempfile.TemporaryDirectory()
    work, skipped = header_jobs(jobs, src_root, tmp.name, args.exclude)
    if args.sources:
        work += [(s, f, c, s) for s, f, c in jobs
                 if s.startswith(src_root) and not (args.exclude and re.search(args.exclude, s))]
    for s in skipped:
        print(f"skip {os.path.relpath(s, src_root)}")

    t0 = time.monotonic()
    fails, times = [], []
    with ThreadPoolExecutor(args.j) as ex:
        for fut in as_completed([ex.submit(run_one, *w) for w in work]):
            name, rc, dt, err = fut.result()
            times.append((dt, name))
            if rc:
                fails.append(name)
                print(f"FAIL {os.path.relpath(name, src_root)}\n{err[:3000]}", flush=True)

    if args.times:
        with open(args.times, "w") as fp:
            fp.writelines(f"{dt:8.2f} {name}\n" for dt, name in sorted(times, reverse=True))
    print(f"\nchecked {len(work)} TUs ({len(skipped)} headers skipped), {len(fails)} failed; "
          f"wall {time.monotonic() - t0:.0f} s at -j{args.j}")
    for name in sorted(fails):
        print(f"  failed: {os.path.relpath(name, src_root)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
