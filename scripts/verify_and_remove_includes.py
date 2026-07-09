#!/usr/bin/env python3
"""Compile-verified removal of unused #include lines flagged by clang-include-cleaner.

For each target file, candidate unused includes (as reported in unused_includes.txt,
produced by scripts/find_unused_includes.py) are removed one at a time, in descending
line-number order, and the translation unit is recompiled after each removal using its
exact command from compile_commands.json (redirected to an isolated temp object file,
so the real build directory is never touched). A removal is kept only if the file still
compiles; otherwise the line is restored and logged as a false positive / genuinely
needed include.

Successful files (>=1 include actually removed) are committed individually to git so
the change is bisectable. Everything is logged to logs/<prefix>.json and
logs/<prefix>.md for reproducibility.

Usage:
    python3 scripts/verify_and_remove_includes.py --top 20 --commit
    python3 scripts/verify_and_remove_includes.py --files path/to/a.cc path/to/b.cc
"""
import argparse
import json
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path("/home/tamago/psi4")
BUILD_DIR = REPO_ROOT / "objdir_psidev/psi4-core-prefix/src/psi4-core-build"
COMPILE_COMMANDS = BUILD_DIR / "compile_commands.json"
DEFAULT_REPORT = REPO_ROOT / "unused_includes.txt"
TMP_OBJ_DIR = Path("/tmp/verify_and_remove_obj")

INCLUDE_LINE_RE = re.compile(r'^\s*#\s*include\s*[<"]')
CANDIDATE_RE = re.compile(r'^\s*-\s*(.+?)\s*@Line:(\d+)\s*$')


def load_compile_db():
    entries = json.loads(COMPILE_COMMANDS.read_text())
    return {e["file"]: e for e in entries}


def load_candidates(report_path):
    """Returns dict: absolute file path -> list of (line_no, header_text), plus ordered file list."""
    text = Path(report_path).read_text()
    blocks = text.split("\n\n")
    candidates = {}
    order = []
    for b in blocks:
        lines = b.splitlines()
        if not lines or not lines[0].startswith("/"):
            continue
        fname = lines[0].strip()
        cands = []
        for line in lines[1:]:
            m = CANDIDATE_RE.match(line)
            if m:
                header, lineno = m.group(1), int(m.group(2))
                cands.append((lineno, header))
        if cands:
            candidates[fname] = cands
            order.append((len(cands), fname))
    order.sort(reverse=True)
    return candidates, [f for _, f in order]


def build_argv(entry, tmp_out):
    argv = shlex.split(entry["command"])
    try:
        oidx = argv.index("-o")
        argv[oidx + 1] = str(tmp_out)
    except ValueError:
        argv += ["-o", str(tmp_out)]
    return argv, entry["directory"]


def compile_check(entry, tmp_out):
    argv, cwd = build_argv(entry, tmp_out)
    start = time.time()
    proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True)
    elapsed = time.time() - start
    ok = proc.returncode == 0
    stderr_tail = "\n".join(proc.stderr.strip().splitlines()[-25:]) if proc.stderr else ""
    return ok, stderr_tail, elapsed


def process_file(path, entry, cands, tmp_out, do_commit):
    p = Path(path)
    original_text = p.read_text()
    lines = original_text.split("\n")

    record = {"file": path, "candidates": [], "removed": [], "kept": [], "errors": []}

    # descending line order: removing a line never shifts not-yet-processed candidates above it
    for lineno, header in sorted(cands, key=lambda x: -x[0]):
        idx = lineno - 1
        if idx < 0 or idx >= len(lines):
            record["errors"].append({"line": lineno, "header": header, "reason": "line index out of range"})
            continue
        current_line = lines[idx]
        if not INCLUDE_LINE_RE.match(current_line):
            record["errors"].append({
                "line": lineno, "header": header,
                "reason": f"line content mismatch, expected #include, got: {current_line!r}",
            })
            continue

        removed_line = lines[idx]
        del lines[idx]
        p.write_text("\n".join(lines))

        ok, stderr_tail, elapsed = compile_check(entry, tmp_out)

        if ok:
            record["removed"].append({"line": lineno, "header": header, "text": removed_line.strip(), "compile_seconds": round(elapsed, 2)})
        else:
            lines.insert(idx, removed_line)
            p.write_text("\n".join(lines))
            record["kept"].append({
                "line": lineno, "header": header, "text": removed_line.strip(),
                "reason": "compile failed without it", "compiler_error_tail": stderr_tail,
                "compile_seconds": round(elapsed, 2),
            })

    # final sanity recompile of resulting file
    final_ok, final_err, final_elapsed = compile_check(entry, tmp_out)
    record["final_compile_ok"] = final_ok
    record["final_compile_seconds"] = round(final_elapsed, 2)
    if not final_ok:
        # Should not happen given incremental verification, but guard anyway: fully restore.
        p.write_text(original_text)
        record["final_compile_error_tail"] = final_err
        record["reverted_to_original"] = True
        record["removed"] = []
        return record

    record["reverted_to_original"] = False

    if record["removed"] and do_commit:
        rel = str(p.relative_to(REPO_ROOT))
        subprocess.run(["git", "add", rel], cwd=REPO_ROOT, check=True)
        n = len(record["removed"])
        msg = f"Remove {n} unused include{'s' if n != 1 else ''} from {rel} (clang-include-cleaner + compile-verified)"
        commit = subprocess.run(["git", "commit", "-m", msg], cwd=REPO_ROOT, capture_output=True, text=True)
        record["committed"] = commit.returncode == 0
        if commit.returncode == 0:
            sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
            record["commit_sha"] = sha
    else:
        record["committed"] = False

    return record


def write_logs(records, log_prefix):
    logs_dir = REPO_ROOT / "logs"
    logs_dir.mkdir(exist_ok=True)
    json_path = logs_dir / f"{log_prefix}.json"
    md_path = logs_dir / f"{log_prefix}.md"

    json_path.write_text(json.dumps(records, indent=2))

    total_removed = sum(len(r["removed"]) for r in records)
    total_kept = sum(len(r["kept"]) for r in records)
    files_changed = sum(1 for r in records if r["removed"])

    lines = []
    lines.append(f"# Unused includes verification pass — {log_prefix}")
    lines.append("")
    lines.append(f"- Files processed: {len(records)}")
    lines.append(f"- Files with >=1 include actually removed: {files_changed}")
    lines.append(f"- Total includes removed (compile-verified): {total_removed}")
    lines.append(f"- Total includes kept (false positive / needed): {total_kept}")
    lines.append("")
    lines.append("| File | Removed | Kept (needed) | Commit |")
    lines.append("|---|---|---|---|")
    for r in records:
        sha = r.get("commit_sha", "-")
        lines.append(f"| `{r['file']}` | {len(r['removed'])} | {len(r['kept'])} | {sha} |")
    lines.append("")

    for r in records:
        if not r["removed"] and not r["kept"]:
            continue
        lines.append(f"## `{r['file']}`")
        if r["removed"]:
            lines.append("")
            lines.append("Removed:")
            for item in r["removed"]:
                lines.append(f"- `{item['header']}` (was line {item['line']})")
        if r["kept"]:
            lines.append("")
            lines.append("Kept (compile failed without it):")
            for item in r["kept"]:
                lines.append(f"- `{item['header']}` (line {item['line']})")
        lines.append("")

    md_path.write_text("\n".join(lines))
    return json_path, md_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default=str(DEFAULT_REPORT))
    ap.add_argument("--top", type=int, default=None, help="process the top N worst-offender files")
    ap.add_argument("--files", nargs="*", default=None, help="explicit list of absolute file paths")
    ap.add_argument("--commit", action="store_true", help="git commit each successfully cleaned file")
    ap.add_argument("--log-prefix", default=None)
    args = ap.parse_args()

    TMP_OBJ_DIR.mkdir(parents=True, exist_ok=True)
    compile_db = load_compile_db()
    candidates, ranked_order = load_candidates(args.report)

    if args.files:
        target_files = args.files
    elif args.top:
        target_files = ranked_order[: args.top]
    else:
        target_files = ranked_order

    log_prefix = args.log_prefix or f"unused_includes_pass1_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    records = []
    for i, f in enumerate(target_files, 1):
        if f not in compile_db:
            print(f"[{i}/{len(target_files)}] {f} -> SKIP (not in compile_commands.json)", file=sys.stderr)
            continue
        entry = compile_db[f]
        tmp_out = TMP_OBJ_DIR / (Path(f).name + ".o")
        rec = process_file(f, entry, candidates[f], tmp_out, args.commit)
        print(
            f"[{i}/{len(target_files)}] {f} -> removed {len(rec['removed'])}, kept {len(rec['kept'])}"
            + (f", COMMIT {rec.get('commit_sha')}" if rec.get("committed") else ""),
            file=sys.stderr,
        )
        records.append(rec)

    json_path, md_path = write_logs(records, log_prefix)
    print(f"\nLogs written to {json_path} and {md_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
