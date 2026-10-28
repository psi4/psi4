#!/usr/bin/env python3

import argparse
import ast
import ctypes
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple


XC_HEADER_RE = re.compile(r"^\s*#define\s+XC_([A-Z0-9_]+_XC_[A-Z0-9_]+)\s+\d+", re.MULTILINE)
XC_FLAGS_NEEDS_LAPLACIAN = 1 << 15
XC_FLAGS_NEEDS_TAU = 1 << 16
XC_FLAGS_DEVELOPMENT = 1 << 14


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[4]
    default_dft_dir = root / "psi4" / "driver" / "procrouting" / "dft"
    default_libxc_file = default_dft_dir / "libxc_functionals.py"

    parser = argparse.ArgumentParser(
        description=(
            "Sync LibXC full XC functionals from xc_funcs.h into "
            "psi4/driver/procrouting/dft/libxc_functionals.py"
        )
    )
    parser.add_argument(
        "--header",
        required=True,
        type=Path,
        help="Path to libxc xc_funcs.h",
    )
    parser.add_argument(
        "--dft-dir",
        default=default_dft_dir,
        type=Path,
        help="Path to psi4/driver/procrouting/dft directory",
    )
    parser.add_argument(
        "--libxc-file",
        default=default_libxc_file,
        type=Path,
        help="Path to libxc_functionals.py to update",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply edits to libxc_functionals.py (default is dry-run)",
    )
    parser.add_argument(
        "--no-safety-check",
        action="store_true",
        help="Skip probing XC_build for candidate functionals",
    )
    parser.add_argument(
        "--no-runtime-check",
        action="store_true",
        help="Skip runtime SCF probe for proposed additions",
    )
    parser.add_argument(
        "--no-runtime-check-uhf",
        action="store_true",
        help="Skip unrestricted runtime SCF probe (UHF) for proposed additions",
    )
    parser.add_argument(
        "--allow-development",
        action="store_true",
        help="Allow LibXC functionals marked XC_FLAGS_DEVELOPMENT",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional path to write a JSON report",
    )
    return parser.parse_args(argv)


def _parse_ast(path: Path) -> ast.AST:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _dict_get_constant_str(dct: ast.Dict, key: str) -> Optional[str]:
    for k, v in zip(dct.keys, dct.values):
        if isinstance(k, ast.Constant) and isinstance(k.value, str) and k.value == key:
            if isinstance(v, ast.Constant) and isinstance(v.value, str):
                return v.value
    return None


def _dict_get_aliases(dct: ast.Dict) -> List[str]:
    aliases: List[str] = []
    for k, v in zip(dct.keys, dct.values):
        if isinstance(k, ast.Constant) and isinstance(k.value, str) and k.value == "alias":
            if isinstance(v, (ast.List, ast.Tuple)):
                for elt in v.elts:
                    if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                        aliases.append(elt.value)
    return aliases


def _dict_get_xc_keys(dct: ast.Dict) -> List[str]:
    keys: List[str] = []
    for k, v in zip(dct.keys, dct.values):
        if isinstance(k, ast.Constant) and isinstance(k.value, str) and k.value == "xc_functionals":
            if isinstance(v, ast.Dict):
                for xk in v.keys:
                    if isinstance(xk, ast.Constant) and isinstance(xk.value, str):
                        keys.append(xk.value)
    return keys


def collect_registered_data(dft_dir: Path) -> Tuple[Set[str], Set[str]]:
    registered_xc_keys: Set[str] = set()
    registered_names: Set[str] = set()

    for pyf in sorted(dft_dir.glob("*_functionals.py")):
        try:
            tree = _parse_ast(pyf)
        except SyntaxError:
            continue

        for node in ast.walk(tree):
            if isinstance(node, ast.Dict):
                for xc_key in _dict_get_xc_keys(node):
                    if "_XC_" in xc_key:
                        registered_xc_keys.add(xc_key.upper())

                name = _dict_get_constant_str(node, "name")
                if name:
                    registered_names.add(name.lower())
                for alias in _dict_get_aliases(node):
                    registered_names.add(alias.lower())

    return registered_xc_keys, registered_names


def parse_header_xc_tokens(header_path: Path) -> List[str]:
    text = header_path.read_text(encoding="utf-8", errors="ignore")
    tokens = [m.group(1).upper() for m in XC_HEADER_RE.finditer(text)]
    # Keep order while dropping duplicates.
    return list(dict.fromkeys(tokens))


def proposed_method_name(token: str) -> str:
    return token.split("_XC_", 1)[1].replace("_", "-")


def probe_supported(tokens: Sequence[str]) -> Dict[str, Optional[str]]:
    try:
        import psi4
    except Exception as exc:
        raise RuntimeError(
            "Safety check requested but psi4 import failed. "
            "Run from a psi4 environment or pass --no-safety-check. "
            f"Import error: {type(exc).__name__}: {exc}"
        ) from exc

    reasons: Dict[str, Optional[str]] = {}
    for token in tokens:
        full = "XC_" + token
        try:
            psi4.core.SuperFunctional.XC_build(full, True)
            reasons[token] = None
        except Exception as exc:
            msg = str(exc).strip()
            if msg:
                reasons[token] = f"{type(exc).__name__}: {msg.splitlines()[0]}"
            else:
                reasons[token] = f"{type(exc).__name__}: {repr(exc)}"
    return reasons


def probe_libxc_flags(tokens: Sequence[str]) -> Dict[str, Dict[str, object]]:
    class _xc_func_type(ctypes.Structure):
        pass

    class _xc_func_info_type(ctypes.Structure):
        pass

    lib = ctypes.CDLL("libxc.so")
    lib.xc_func_alloc.restype = ctypes.POINTER(_xc_func_type)
    lib.xc_func_free.argtypes = [ctypes.POINTER(_xc_func_type)]
    lib.xc_func_init.argtypes = [ctypes.POINTER(_xc_func_type), ctypes.c_int, ctypes.c_int]
    lib.xc_func_init.restype = ctypes.c_int
    lib.xc_func_end.argtypes = [ctypes.POINTER(_xc_func_type)]
    lib.xc_functional_get_number.argtypes = [ctypes.c_char_p]
    lib.xc_functional_get_number.restype = ctypes.c_int
    lib.xc_func_get_info.argtypes = [ctypes.POINTER(_xc_func_type)]
    lib.xc_func_get_info.restype = ctypes.POINTER(_xc_func_info_type)
    lib.xc_func_info_get_flags.argtypes = [ctypes.POINTER(_xc_func_info_type)]
    lib.xc_func_info_get_flags.restype = ctypes.c_int
    lib.xc_func_info_get_n_ext_params.argtypes = [ctypes.POINTER(_xc_func_info_type)]
    lib.xc_func_info_get_n_ext_params.restype = ctypes.c_int
    lib.xc_func_info_get_ext_params_name.argtypes = [ctypes.POINTER(_xc_func_info_type), ctypes.c_int]
    lib.xc_func_info_get_ext_params_name.restype = ctypes.c_char_p
    lib.xc_func_info_get_ext_params_description.argtypes = [ctypes.POINTER(_xc_func_info_type), ctypes.c_int]
    lib.xc_func_info_get_ext_params_description.restype = ctypes.c_char_p

    out: Dict[str, Dict[str, object]] = {}
    for token in tokens:
        full = ("XC_" + token).encode("ascii")
        fid = lib.xc_functional_get_number(full)
        if fid <= 0:
            out[token] = {
                "needs_laplacian": False,
                "needs_tau": False,
                "is_development": False,
                "ext_param_names": [],
                "ext_param_descriptions": [],
            }
            continue

        func = lib.xc_func_alloc()
        try:
            rc = lib.xc_func_init(func, fid, 1)
            if rc != 0:
                out[token] = {
                    "needs_laplacian": False,
                    "needs_tau": False,
                    "is_development": False,
                    "ext_param_names": [],
                    "ext_param_descriptions": [],
                }
                continue

            info = lib.xc_func_get_info(func)
            flags = lib.xc_func_info_get_flags(info)
            npars = lib.xc_func_info_get_n_ext_params(info)
            ext_names: List[str] = []
            ext_descs: List[str] = []
            for i in range(npars):
                nraw = lib.xc_func_info_get_ext_params_name(info, i)
                draw = lib.xc_func_info_get_ext_params_description(info, i)
                ext_names.append((nraw or b"").decode("utf-8", errors="ignore"))
                ext_descs.append((draw or b"").decode("utf-8", errors="ignore"))
            out[token] = {
                "needs_laplacian": bool(flags & XC_FLAGS_NEEDS_LAPLACIAN),
                "needs_tau": bool(flags & XC_FLAGS_NEEDS_TAU),
                "is_development": bool(flags & XC_FLAGS_DEVELOPMENT),
                "ext_param_names": ext_names,
                "ext_param_descriptions": ext_descs,
            }
        finally:
            try:
                lib.xc_func_end(func)
            except Exception:
                pass
            lib.xc_func_free(func)

    return out


def line_for_entry(method_name: str, token: str) -> str:
    name_pad = max(0, 15 - len(method_name))
    token_pad = max(0, 26 - len(token))
    return (
        'funcs.append({"name": "'
        + method_name
        + '"'
        + (' ' * name_pad)
        + ' , "xc_functionals": {"'
        + token
        + '"'
        + (' ' * token_pad)
        + ' : {}}})\n'
    )


def runtime_probe_token(token: str, reference: str) -> Optional[str]:
    code = textwrap.dedent(
        f"""
        import math
        import psi4

        psi4.core.set_output_file('/tmp/psi4_libxc_probe.out', False)
        psi4.set_options({{
            'basis': 'sto-3g',
            'reference': '{reference}',
            'scf_type': 'df',
            'e_convergence': 1e-6,
            'd_convergence': 1e-6,
            'maxiter': 30,
        }})
        mol = psi4.geometry('He')
        fd = {{'name': 'probe', 'xc_functionals': {{'{token}': {{}}}}}}
        e = psi4.energy('scf', molecule=mol, dft_functional=fd)
        if not math.isfinite(float(e)):
            raise RuntimeError('Non-finite energy from runtime probe')
        print('OK', float(e))
        """
    )

    proc = subprocess.run(
        [sys.executable, "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd="/tmp",
    )

    if proc.returncode == 0:
        return None

    out = (proc.stderr or "").strip()
    if not out:
        out = (proc.stdout or "").strip()
    if out:
        first = out.splitlines()[-1]
        return f"runtime probe failed ({reference}, exit {proc.returncode}): {first}"
    return f"runtime probe failed ({reference}, exit {proc.returncode})"


def apply_insertions(libxc_file: Path, lines: List[str]) -> None:
    text_lines = libxc_file.read_text(encoding="utf-8").splitlines(keepends=True)
    marker_idx = None
    for i, ln in enumerate(text_lines):
        if "# yapf: enable" in ln:
            marker_idx = i
            break
    if marker_idx is None:
        raise RuntimeError(f"Could not find '# yapf: enable' marker in {libxc_file}")

    new_text = text_lines[:marker_idx] + lines + text_lines[marker_idx:]
    libxc_file.write_text("".join(new_text), encoding="utf-8")


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)

    header_path = args.header.resolve()
    dft_dir = args.dft_dir.resolve()
    libxc_file = args.libxc_file.resolve()

    if not header_path.exists():
        raise FileNotFoundError(f"Header path does not exist: {header_path}")
    if not dft_dir.exists():
        raise FileNotFoundError(f"DFT directory does not exist: {dft_dir}")
    if not libxc_file.exists():
        raise FileNotFoundError(f"libxc functionals file does not exist: {libxc_file}")

    header_tokens = parse_header_xc_tokens(header_path)
    registered_xc_keys, registered_names = collect_registered_data(dft_dir)

    report: Dict[str, object] = {
        "header": str(header_path),
        "dft_dir": str(dft_dir),
        "libxc_file": str(libxc_file),
        "header_xc_count": len(header_tokens),
        "registered_xc_count": len(registered_xc_keys),
        "already_registered": [],
        "candidates": [],
        "to_add": [],
        "unsupported": [],
        "laplacian_required": [],
        "external_param_required": [],
        "development_filtered": [],
        "runtime_unsupported": [],
        "name_collision": [],
        "batch_collision": [],
    }

    already_registered: List[str] = []
    candidates: List[str] = []
    for tok in header_tokens:
        if tok in registered_xc_keys:
            already_registered.append(tok)
        else:
            candidates.append(tok)

    report["already_registered"] = already_registered
    report["candidates"] = candidates

    safety_failures: Dict[str, Optional[str]] = {}
    flags_data: Dict[str, Dict[str, object]] = {}
    flag_probe_error: Optional[str] = None
    if not args.no_safety_check and candidates:
        try:
            flags_data = probe_libxc_flags(candidates)
        except Exception as exc:
            flag_probe_error = f"{type(exc).__name__}: {exc}"
        safety_failures = probe_supported(candidates)

    new_name_to_token: Dict[str, str] = {}
    to_add: List[Tuple[str, str]] = []
    unsupported: List[Dict[str, str]] = []
    laplacian_required: List[Dict[str, str]] = []
    external_param_required: List[Dict[str, str]] = []
    development_filtered: List[Dict[str, str]] = []
    runtime_unsupported: List[Dict[str, str]] = []
    name_collision: List[Dict[str, str]] = []
    batch_collision: List[Dict[str, str]] = []

    for tok in candidates:
        if tok in flags_data and flags_data[tok].get("needs_laplacian", False):
            reason = "requires Laplacian (XC_FLAGS_NEEDS_LAPLACIAN), which Psi4 does not support"
            if flags_data[tok].get("needs_tau", False):
                reason += "; functional also uses tau"
            laplacian_required.append({"token": tok, "reason": reason})
            continue

        if tok in flags_data:
            non_internal = [p for p in flags_data[tok].get("ext_param_names", []) if isinstance(p, str) and p and not p.startswith("_")]
            descs = [d for d in flags_data[tok].get("ext_param_descriptions", []) if isinstance(d, str)]
            desc_hint = next((d for d in descs if "must be calculated by the calling program" in d.lower()), "")
            if non_internal or desc_hint:
                reason = "requires external parameter(s) not auto-managed by Psi4"
                if non_internal:
                    reason += f": {', '.join(non_internal)}"
                if desc_hint:
                    reason += f" ({desc_hint})"
                external_param_required.append({"token": tok, "reason": reason})
                continue

            if flags_data[tok].get("is_development", False) and not args.allow_development:
                development_filtered.append({
                    "token": tok,
                    "reason": "marked XC_FLAGS_DEVELOPMENT in LibXC",
                })
                continue

        if tok in safety_failures and safety_failures[tok] is not None:
            unsupported.append({"token": tok, "reason": str(safety_failures[tok])})
            continue

        method_name = proposed_method_name(tok)
        lower = method_name.lower()

        if lower in registered_names:
            name_collision.append({
                "token": tok,
                "method_name": method_name,
                "reason": "method name or alias already registered",
            })
            continue

        if lower in new_name_to_token:
            batch_collision.append({
                "token": tok,
                "method_name": method_name,
                "reason": f"collides with candidate token {new_name_to_token[lower]}",
            })
            continue

        new_name_to_token[lower] = tok
        to_add.append((method_name, tok))

    if not args.no_safety_check and not args.no_runtime_check and to_add:
        kept: List[Tuple[str, str]] = []
        for method_name, tok in to_add:
            err = runtime_probe_token(tok, "rhf")
            if err is None and not args.no_runtime_check_uhf:
                err = runtime_probe_token(tok, "uhf")
            if err is None:
                kept.append((method_name, tok))
            else:
                runtime_unsupported.append({"token": tok, "method_name": method_name, "reason": err})
        to_add = kept

    report["to_add"] = [{"method_name": m, "token": t} for m, t in to_add]
    report["unsupported"] = unsupported
    report["laplacian_required"] = laplacian_required
    report["external_param_required"] = external_param_required
    report["development_filtered"] = development_filtered
    report["runtime_unsupported"] = runtime_unsupported
    report["name_collision"] = name_collision
    report["batch_collision"] = batch_collision
    if flag_probe_error:
        report["flag_probe_error"] = flag_probe_error

    print(f"Header _XC_ tokens:           {len(header_tokens)}")
    print(f"Already registered XC keys:   {len(already_registered)}")
    print(f"Candidate new XC keys:        {len(candidates)}")
    print(f"Need-Laplacian filtered:      {len(laplacian_required)}")
    print(f"Needs external params:        {len(external_param_required)}")
    print(f"Development filtered:         {len(development_filtered)}")
    print(f"Unsupported by XC_build:      {len(unsupported)}")
    print(f"Unsupported by runtime probe: {len(runtime_unsupported)}")
    print(f"Name collisions:              {len(name_collision)}")
    print(f"Batch collisions:             {len(batch_collision)}")
    print(f"New entries to add:           {len(to_add)}")
    if flag_probe_error:
        print(f"Flag probe warning:           {flag_probe_error}")

    if laplacian_required:
        print("\nRequires Laplacian (excluded):")
        for row in laplacian_required:
            print(f"  - {row['token']}: {row['reason']}")

    if external_param_required:
        print("\nNeeds external params (excluded):")
        for row in external_param_required:
            print(f"  - {row['token']}: {row['reason']}")

    if development_filtered:
        print("\nDevelopment functionals (excluded):")
        for row in development_filtered:
            print(f"  - {row['token']}: {row['reason']}")

    if unsupported:
        print("\nUnsupported functionals:")
        for row in unsupported:
            print(f"  - {row['token']}: {row['reason']}")

    if runtime_unsupported:
        print("\nRuntime probe excluded:")
        for row in runtime_unsupported:
            print(f"  - {row['token']} -> {row['method_name']}: {row['reason']}")

    if name_collision:
        print("\nName collisions:")
        for row in name_collision:
            print(f"  - {row['token']} -> {row['method_name']}: {row['reason']}")

    if batch_collision:
        print("\nBatch collisions:")
        for row in batch_collision:
            print(f"  - {row['token']} -> {row['method_name']}: {row['reason']}")

    if to_add:
        print("\nProposed additions:")
        for method_name, tok in to_add:
            print(f"  - {method_name} [{tok}]")

    if args.apply and to_add:
        insertion_lines = [line_for_entry(method_name, tok) for method_name, tok in to_add]
        apply_insertions(libxc_file, insertion_lines)
        print(f"\nUpdated {libxc_file} with {len(insertion_lines)} new entries.")
    elif args.apply:
        print("\nNo updates applied (no new entries to add).")

    if args.report_json:
        args.report_json.parent.mkdir(parents=True, exist_ok=True)
        args.report_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"JSON report written to {args.report_json}")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
