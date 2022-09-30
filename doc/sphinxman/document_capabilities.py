#!/usr/bin/env python

import ast
import argparse
import unicodedata
from pathlib import Path

from psi4.driver.procrouting.proc_data import method_governing_type_keywords
from psi4.driver.p4util.exceptions import sanitize_method

# TODO how to handle this?
method_governing_type_keywords["svwn"] = "scf_type"
method_governing_type_keywords["pbe"] = "scf_type"
method_governing_type_keywords["b3lyp"] = "scf_type"
method_governing_type_keywords["wb97x"] = "scf_type"
method_governing_type_keywords["b2plyp"] = "scf_type"

parser = argparse.ArgumentParser()
parser.add_argument("--sphinx", action="store_true", help="Sphinx format vs. eyes format")
parser.add_argument("--mode", default="details", choices=["details", "summary", "ccenergy", "fnocc", "dfmp2", "occ", "occ_oo", "occ_nonoo", "scf"], help="summary table vs. detailed (per-module) table vs. single-module table")
parser.add_argument("--stdsuite", default="stdsuite_psi4.txt", help="where is recorder file written by standard_suite_runner.py")
parser.add_argument("--writefile", action="append", help="where should table file be saved? defaults to autodoc_capabilities_{mode}.rst . repeat arg for multiple outputs")
parser.add_argument("--quiet", action="store_true", help="Don't also print to screen")
parser.add_argument("--second-in-rst", action="store_true", help="Suppress Sphinx substitution definitions to avoid warnings if multiple tables included in rst file")
parser.add_argument("--driver", default="eg", choices=["e", "eg", "egh"], help="maximum deriv column")
args = parser.parse_args()

blank = " "
equal = "="
dashh = "-"

notes = []

trans_cell = {
    "pass":   ('\u2713',       "|y|"),
    "(pass)": ('\u2713\u0332', "|d|"),
    "[pass]": ('\u2713\u0333', "|g|"),
    "fd":     ('\u2237',       "|e|"),
    "(fd)":   ('\u2237\u0332', "|b|"),
    "[fd]":   ('\u2237\u0333', "|c|"),
    "error":  ('',             ""   ),  # toggle for " " in visual mode
    # "error":  ('\u2717',       ""   ),  # toggle for "X" in visual mode
    # "wrong":  ('\u25a0')
}

# need indep letter for footnotes `[#f1]_` for each table so multiple can be included on one page, plus "h" for common headers. used are: a, c, d, e, f, m, n, o, r, s
fn_lbl = {"details": "d", "summary": "s", "ccenergy": "e", "fnocc": "f", "dfmp2": "m", "occ": "o", "occ_oo": "c", "occ_nonoo": "n", "scf": "r"}[args.mode]

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


def fcell(module, winnowed, contrast=False):
    global notes

    if winnowed:
        status = sorted(set([entry["status"] for entry in winnowed]), key=len)
        if len(winnowed) > 1 and len(status) > 1:
            if module.startswith("aaaa-") and status == ["fd", "pass"]:
                pass
            else:
                print("trouble:")
                for entry in winnowed:
                    print("\t", entry)

        cell = status[-1]  # prefer "pass" (grd1) to "fd" (grd0)
        note = sorted([entry["note"] for entry in winnowed], key=len)[-1]  # prefer "default" to ""
        if note == "defaultdefault":
            cell = "[" + cell + "]"
        elif note == "default":
            cell = "(" + cell + ")"

        if False:
        #if note:
            notes.append(note)
            notes = unique(notes)
            idx = notes.index(note) + 10
            cell += f" [#f{idx}]_"

        cell = trans_cell[cell][args.sphinx]
        if contrast:
            cell = cell.upper()

    else:
        cell = ""

    return cell


def fheader(head, width, contrast=False):
    trans = {
        # eghs
        "energy": ["|energy_fn|", "Energy", "E"],
        "gradient": [f"|gradient_fn|\ [#{fn_lbl}3]_", "Gradient", "G"],
        "hessian": [f"|hessian_fn|\ [#{fn_lbl}4]_", "Hessian", "H"],
        # fcae
        "ae": ["All-Electron", "AE", "A"],
        "fc": ["Frozen-Core", "FC", "F"],
        # refs
        "rhf": ["Restricted (RHF)", "RHF"],
        "uhf": ["Unrestricted (UHF)", "UHF"],
        "rohf": ["Restricted Open (ROHF)", "ROHF"],
        # types
        "conv": ["Conventional", "CV", "CONV", "CV"],
        "df": ["Density-Fit", "DF"],
        "cd": ["Cholesky-Decomp.", "CD"],
    }
    contrasted_trans = {k: v for k, v in trans.items()}
    contrasted_trans["ae"] = ["|A|"]
    contrasted_trans["fc"] = ["|F|"]

    if head in trans:
        source = contrasted_trans if contrast else trans
        for choice in source[head]:
            if len(choice) + 2 <= width:
                return choice
        return source[head][-1]
    else:
        return head


def visible_width(chars):
    return sum(not unicodedata.combining(ch) for ch in chars)


def compute_width(head, width):
    return width + len(head) - visible_width(head)


def fline_data(outer, mid, inner, guts, *, span=1, contrast=1):
    vertex = "|"
    eff_width = cwidth * span + span - 1

    sguts = [f"{fheader(item, eff_width):^{compute_width(item, eff_width)}}" for item in guts]

    if args.sphinx:
        csguts = [f"{fheader(item, eff_width, contrast=True):^{compute_width(item, eff_width)}}" for item in guts]
        lguts = [(csguts[i:i + contrast] if (b % 2) else sguts[i:i + contrast]) for b, i in enumerate(range(0, len(sguts), contrast))]
        lguts = [vertex.join(i) for i in lguts]
    else:
        lguts = [sguts[i:i + contrast] for i in range(0, len(sguts), contrast)]
        lguts = [blank.join(i) for i in lguts]

    datasection = vertex.join(lguts)

    outersection = vertex + f"{' ' + outer:<{l1width}}"
    midsection = (vertex + f"{' ' + mid:<{l2width}}") if l2width else ""  # mid header column is optional
    innersection = vertex + f"{' ' + inner:<{l3width}}"

    return blank * pwidth + outersection + midsection + innersection + vertex + datasection + vertex


def fline_fill(outer, mid, inner, guts, span=1):
    vertex = "+"
    eff_width = cwidth * span + span - 1

    sguts = [f"{item*eff_width}" for item in guts]
    datasection = vertex.join(sguts)

    outersection = vertex + f"{outer*l1width}"
    midsection = (vertex + f"{mid*l2width}") if l2width else ""  # mid header column is optional
    innersection = vertex + f"{inner*l3width}"

    return blank * pwidth + outersection + midsection + innersection + vertex + datasection + vertex


def extract_modules(winnowed):
    return sorted(list(set([entry["module"] for entry in winnowed])))


methods = [
"hf",
"mp2",
"mp2.5",
"mp3",
"mp4(sdq)",
"mp4",
"zapt2",
"cisd",
"qcisd",
"qcisd(t)",
"fci",
"remp2",
"lccd",
"lccsd",
"cepa(1)",
"cepa(3)",
"acpf",
"aqcc",
"ccd",
"bccd",
"cc2",
"ccsd",
"ccsd+t(ccsd)",
"ccsd(t)",
"a-ccsd(t)",
"bccd(t)",
"cc3",
"ccsdt-1a",
"ccsdt-1b",
"ccsdt-2",
"ccsdt-3",
"ccsdt",
"ccsdt(q)",
"ccsdtq",
"omp2",
"omp2.5",
"omp3",
"oremp2",
"olccd",
"svwn",
"pbe",
"b3lyp",
"b3lyp5",
"wb97x",
"b2plyp",
]

notes_holder = {
    ("lccsd", None): (", cepa(0)", None, None),
    ("a-ccsd(t)", None): ("FN", "a-CCSD(T) also known as CCSD(aT), Lambda-CCSD(T), and CCSD(T)_L", None),
    ("mp2", "DFMP2"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("omp2", "OCC"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("omp2.5", "OCC"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("omp3", "OCC"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("oremp2", "OCC"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("olccd", "OCC"): ("FN", "Also available for DFT references RKS/UKS", None),
    ("svwn", None): (", LSDA DFT", None, None),
    ("pbe", None): (", GGA DFT", None, None),
    ("b3lyp", None): (", Hybrid DFT", None, None),
    ("wb97x", None): (", LRC DFT", None, None),
    ("b2plyp", None): (", DH DFT", "DH-DFT only available with DF-MP2", None),
    ("cisd", None): (", ci\ *n*", "Arbitrary-order *n* through DETCI is inefficient byproduct of CI", ["fnocc"]),
    ("zapt2", None): (", zapt\ *n*", "Arbitrary-order *n* through DETCI is inefficient byproduct of CI", None),
    ("mp4", None): (", mp\ *n*", "Arbitrary-order *n* through DETCI is inefficient byproduct of CI", ["fnocc"]),
}


def method_title_append(mtd, mod):
    global notes

    keys = [(mtd, None)]
    if mod is not None:
        keys.append((mtd, mod))

    for key in keys:
        if key in notes_holder:
            append, note, skip = notes_holder[key]
            if skip and args.mode in skip:
                continue
            if append != "FN":
                mtd += append
            if note:
                notes.append(note)
                notes = unique(notes)
                idx = notes.index(note) + 10
                mtd += f"\ [#{fn_lbl}{idx}]_"

    return mtd


def module_title_append(mtd, mod):
    global notes

    keys = [(mtd, mod)]
    for key in keys:
        if key in notes_holder:
            append, note, skip = notes_holder[key]
            if skip and args.mode in skip:
                continue
            if append != "FN":
                mtd += append
            if note:
                notes.append(note)
                notes = unique(notes)
                idx = notes.index(note) + 10
                mod += f"\ [#{fn_lbl}{idx}]_"

    return mod


with open(args.stdsuite, "r") as fp:
    contents = fp.readlines()

stuff = [ast.literal_eval(ln) for ln in contents]


methods.extend([entry["method"] for entry in stuff])
methods = unique(methods)
fcaes = ["ae", "fc"]
if args.driver == "e":
    eghs = ["energy"]
elif args.driver == "eg":
    eghs = ["energy", "gradient"]
elif args.driver == "egh":
    eghs = ["energy", "gradient", "hessian"]
refs = ["rhf", "uhf", "rohf"]
corl_types = ["conv", "df", "cd"]
natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
natural_ref_rev = {v: k for k, v in natural_ref.items()}


# left margin, left outer header, left mid header (can be zero), left inner header, and body column widths
pwidth = 3
if args.mode == "summary":
    l1width, l2width, l3width = 26, 22, 26
elif args.mode == "details":
    l1width, l2width, l3width = 26, 0, 26
else:
    l1width, l2width, l3width = 26, 0, 25
cwidth = 3  # 13 minimum for footnotes


def table_builder__ref_driver_type_fcae():

    #  +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    #  |               Restricted (RHF)                |              Unrestricted (UHF)               |            Restricted Open (ROHF)             |
    #  +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    #  |     |psi4.energy|     |    |psi4.gradient|    |     |psi4.energy|     |    |psi4.gradient|    |     |psi4.energy|     |    |psi4.gradient|    |
    #  +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    #  |  CV   |  DF   |  CD   |  CV   |  DF   |  CD   |  CV   |  DF   |  CD   |  CV   |  DF   |  CD   |  CV   |  DF   |  CD   |  CV   |  DF   |  CD   |
    #  +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    #  | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F | A | F |
    #  +===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+

    place = '\u25fb'
    down_arrow = '\u2193'
    right_arrow = '\u2192'

    # row headers
    nameh = "name " + down_arrow + " " + right_arrow
    qcmoduleh = "|qc_module| " + down_arrow
    # column headers
    referenceh = "|scf__reference| " + right_arrow
    driverh = ""
    typeh = f"type\ [#{fn_lbl}1]_ " + down_arrow + " " + right_arrow
    fcaeh = f"|freeze_core|\ [#{fn_lbl}2]_\ " + right_arrow

    lines = []

    if not args.sphinx:
        lines.append(".. NOTE: this file is autogenerated for preview and not used in docs directly.")
        lines.append("")

    if args.sphinx and not args.second_in_rst:
        lines.append("")
        lines.append('.. role:: gbg')
        lines.append('.. raw:: html')
        lines.append('')
        lines.append('   <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>')
        lines.append('   <script>')
        lines.append('     $(document).ready(function() {')
        lines.append('       $(".gbg").parent().addClass("gbg-parent");')
        lines.append('     });')
        lines.append('   </script>')
        lines.append('   <style>')
        lines.append('      .gbg-parent {background-color:#cceecc;}')
        lines.append('   </style>')
        lines.append("")

        lines.append("")
        lines.append(".. |energy_fn|   replace:: :func:`~psi4.driver.energy()`")
        lines.append(".. |gradient_fn| replace:: :func:`~psi4.driver.gradient()`")
        lines.append(".. |hessian_fn|  replace:: :func:`~psi4.driver.hessian()`")
        lines.append(".. |freeze_core| replace:: |globals__freeze_core|")
        lines.append(".. |qc_module|   replace:: |globals__qc_module|")
        lines.append(".. |A| replace:: :gbg:`A`")
        lines.append(".. |F| replace:: :gbg:`F`")
        lines.append(".. |y| unicode:: U+2713        .. tick")
        lines.append(".. |d| unicode:: U+2713 U+0332 .. underlined tick")
        lines.append(".. |g| unicode:: U+2713 U+0333 .. double underlined tick")
        lines.append(".. |e| unicode:: U+2237        .. stencil")
        lines.append(".. |b| unicode:: U+2237 U+0332 .. underlined stencil")
        lines.append(".. |c| unicode:: U+2237 U+0333 .. double underlined stencil")
        lines.append(".. |Y| replace:: :gbg:`✓`")
        lines.append(".. |D| replace:: :gbg:`✓̲`")
        lines.append(".. |G| replace:: :gbg:`✓̳`")
        lines.append(".. |E| replace:: :gbg:`∷`")
        lines.append(".. |B| replace:: :gbg:`∷̲`")
        lines.append(".. |C| replace:: :gbg:`∷̳`")
        lines.append("")

    legend_lines = []

    if args.mode == "summary":
        lines.append(".. _`table:stdsuite`:")
        legend_lines.append(f""".. table:: Summary capabilities of |PSIfour|.""")
    elif args.mode == "details":
        lines.append(".. _`table:managedmethods`:")
        legend_lines.append(f""".. table:: Capabilities of |PSIfour|, including details of overlapping modules.""")
    elif args.mode == "ccenergy":
        lines.append(".. _`table:ccenergy_stdsuite`:")
        legend_lines.append(".. table:: Detailed capabilities of CCENERGY and related modules.""")
    elif args.mode == "fnocc":
        lines.append(".. _`table:fnocc_stdsuite`:")
        legend_lines.append(".. table:: Detailed capabilities of the FNOCC module.""")
    elif args.mode == "dfmp2":
        lines.append(".. _`table:dfmp2_stdsuite`:")
        legend_lines.append(".. table:: Detailed capabilities of the DFMP2 module.""")
    elif args.mode == "occ":
        lines.append(".. _`table:occ_stdsuite`:")
        legend_lines.append(".. table:: Detailed capabilities of the OCC module.""")
    elif args.mode == "occ_oo":
        lines.append(".. _`table:occ_stdsuite_oo`:")
        legend_lines.append(".. table:: Detailed orbital-optimized capabilities of the OCC module.""")
    elif args.mode == "occ_nonoo":
        lines.append(".. _`table:occ_stdsuite_nonoo`:")
        legend_lines.append(".. table:: Detailed non-orbital-optimized capabilities of the OCC module.""")
    elif args.mode == "scf":
        lines.append(".. _`table:scf_stdsuite`:")
        legend_lines.append(".. table:: Detailed capabilities of the SCF module.""")
    else:
        pass

    legend_lines.append(f'''"{trans_cell['pass'][args.sphinx]}" is runs analytically.''')

    if args.mode == "summary":
        legend_lines.append(f'''"{trans_cell['fd'][args.sphinx]}" is runs derivative with internal finite difference.''')

    if args.mode == "details":
        legend_lines.append(f"""Single underline "{trans_cell['(pass)'][args.sphinx]}" or "{trans_cell['(fd)'][args.sphinx]}" is default module when |globals__qc_module| unspecified.""")
    elif args.mode == "summary":
        pass
    else:
        legend_lines.append(f"""Single underline "{trans_cell['(pass)'][args.sphinx]}" is default module when |globals__qc_module| unspecified.""")

    if args.mode in ["details", "summary"]:
        legend_lines.append(f"""Double underline "{trans_cell['[pass]'][args.sphinx]}" or "{trans_cell['[fd]'][args.sphinx]}" is default algorithm type when type selector (e.g., |globals__cc_type|\ ) unspecified.""")
    else:
        legend_lines.append(f"""Double underline "{trans_cell['[pass]'][args.sphinx]}" is default algorithm type when type selector (e.g., |globals__cc_type|\ ) unspecified.""")

    lines.append("")
    lines.append(" ".join(legend_lines))

    lines.append("   :align: left")
    lines.append("")
    ncol = len(fcaes) * len(refs) * len(eghs) * len(corl_types)

    # from eghs span below
    contrast_span = 6

    if args.mode == "summary":

        #  +------------+----------------------+--------------------------+
        #  | ◻          |                      | |scf__reference| →       |
        #  +            +                      +                          +
        #  | name ↓ →   |                      | ◻                        |
        #  +            +                      +                          +
        #  | ◻          |                      | type select [#f1]_ ↓ →   |
        #  +            +                      +                          +
        #  | ◻          |                      | |freeze_core| [#f2]_ →   |
        #  +============+======================+==========================+

        lines.append(fline_fill(dashh, dashh, dashh, [dashh] * ncol))
        lines.append(fline_data(place, blank, referenceh, refs, span=len(fcaes)*len(corl_types)*len(eghs)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(nameh, blank, place, len(refs)*eghs, span=len(fcaes)*len(corl_types)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, typeh, len(refs)*len(eghs)*corl_types, span=len(fcaes)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, fcaeh, len(refs)*len(eghs)*len(corl_types)*fcaes, span=1, contrast=contrast_span))
        lines.append(fline_fill(equal, equal, equal, [equal] * ncol))

    elif args.mode == "details":

        #  +-------------------------+--------------------------+
        #  | ◻                       | |globals__qc_module| ↓   |
        #  +                         +                          +
        #  | ◻                       | |scf__reference| →       |
        #  +                         +                          +
        #  | name ↓ →                | ◻                        |
        #  +                         +                          +
        #  | type [#f1]_ ↓ →         | ◻                        |
        #  +                         +                          +
        #  | ◻                       | |freeze_core| [#f2]_ →   |
        #  +=========================+==========================+

        lines.append(fline_fill(dashh, blank, dashh, [dashh] * ncol))
        lines.append(fline_data(place, blank, qcmoduleh, ["|PSIfour| Capabilities"], span=ncol))  # watch out for replace delimiters lining up with cwidth cell delimiters. leads to bad formatting
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, referenceh, refs, span=len(fcaes)*len(corl_types)*len(eghs)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(nameh, blank, place, len(refs)*eghs, span=len(fcaes)*len(corl_types)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(typeh, blank, place, len(refs)*len(eghs)*corl_types, span=len(fcaes)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, fcaeh, len(refs)*len(eghs)*len(corl_types)*fcaes, span=1, contrast=contrast_span))
        lines.append(fline_fill(equal, blank, equal, [equal] * ncol))

    else:
        # single-module

        #  +---------------------+------------------------+
        #  | ◻                   | ◻                      |
        #  +                     +                        +
        #  | ◻                   | |scf__reference| →     |
        #  +                     +                        +
        #  | name ↓ →            | ◻                      |
        #  +                     +                        +
        #  | ◻                   | type [#f1]_ ↓ →        |
        #  +                     +                        +
        #  | ◻                   | |freeze_core| [#f2]_ → |
        #  +=====================+========================+

        module_caption = "OCC" if args.mode in ["occ_oo", "occ_nonoo"] else args.mode.upper()
        spacer = "  " if args.mode in ["dfmp2", "fnocc"] else ""  # avoid aligned delimiters

        lines.append(fline_fill(dashh, dashh, dashh, [dashh] * ncol))
        lines.append(fline_data(place, blank, place, [f"{spacer}{qcmoduleh[:-2]}\ ={module_caption} Capabilities"], span=ncol))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, referenceh, refs, span=len(fcaes)*len(corl_types)*len(eghs)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(nameh, blank, place, len(refs)*eghs, span=len(fcaes)*len(corl_types)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, typeh, len(refs)*len(eghs)*corl_types, span=len(fcaes)))
        lines.append(fline_fill(blank, blank, blank, [dashh] * ncol))
        lines.append(fline_data(place, blank, fcaeh, len(refs)*len(eghs)*len(corl_types)*fcaes, span=1, contrast=contrast_span))
        lines.append(fline_fill(equal, equal, equal, [equal] * ncol))

    for method in methods:
        if (args.mode == "occ_oo" and not method.startswith("o")) or (args.mode == "occ_nonoo" and method.startswith("o")):
            continue

        subset1 = [entry for entry in stuff if entry["method"] == method]
        if len(subset1) == 0:
            continue

        modules = extract_modules(subset1)
        if args.mode == "details":
            pass
        elif args.mode == "summary":
            modules = [modules[0]]
        else:
            if args.mode in ["occ_oo", "occ_nonoo"]:
                select_module = "occ"
            else:
                select_module = args.mode
            modules = [m for m in modules if m == f"psi4-{select_module}"]

            # skip method if module not represented
            if not modules:
                continue

        if args.sphinx:
            # to make Sphinx's alt-line coloring scheme highlight alt-methods, make each method block odd row count
            if (len(modules) % 2) != 1:
                modules.append("zzzz-")

        for module in modules:
            subset2 = [entry for entry in subset1 if entry["module"] == module]

            line = []
            column_block = True

            for ref in refs:
                for egh in eghs:
                    column_block = not column_block
                    for corl_type in corl_types:
                        for fcae in fcaes:
                            subset4 = [entry for entry in subset2 if (entry["fcae"] == fcae
                                                                      and entry["reference"] == ref
                                                                      and entry["driver"] == egh
                                                                      and (natural_ref_rev[entry["scf_type"]] == corl_type if (method_governing_type_keywords[method] == "scf_type") else entry["corl_type"] == corl_type)
                                                                      and entry["scf_type"] == natural_ref[corl_type]
                                                                      and (entry["status"] != "fd" or module.startswith("aaaa-"))
                                                                     )
                                      ]

                            line.append(fcell(module, subset4, column_block))

            module_caption = module[5:].upper()  # filter off "psi4-" in module
            module_title = module_title_append(method, module_caption)
            type_link = f"|globals__{method_governing_type_keywords[method]}|"
            method_title = method_title_append(method, module_caption)
            anchor = fn_lbl * 2
            method_anchor = f".. _{anchor}_{sanitize_method(method)}:"

            if args.mode == "summary":
                l1_col_header = method_title
                l2_col_header = method_anchor
                l3_col_header = type_link
            elif args.mode == "details":
                l2_col_header = ""
                if module == modules[0]:
                    l1_col_header = method_title
                    l3_col_header = method_anchor
                else:
                    if module == modules[1]:
                        l1_col_header = type_link
                    else:
                        l1_col_header = ""
                    l3_col_header = module_title
            else:
                l1_col_header = method_title
                l2_col_header = module_title
                l3_col_header = type_link

            lines.append(fline_data(l1_col_header, l2_col_header, l3_col_header, line, span=1, contrast=contrast_span))
            end_method = "-" if module == modules[-1] else " "
            if args.sphinx:
                lines.append(fline_fill(end_method, "-", "-", ["-"] * ncol))

        if not args.sphinx:
            lines.append(fline_fill("-", "-", "-", ["-"] * ncol))

    lines.append("")
    lines.append(f""".. [#{fn_lbl}1] Algorithm type selection keyword below. Values to the right: conventional ``CV``, density-fitted ``DF``, and Cholesky-decomposed ``CD``.""")
    lines.append(f""".. [#{fn_lbl}2] Active orbital values to the right: all-electron ``A`` and frozen-core ``F``.""")
    if "g" in args.driver:
        if args.mode == "summary":
            lines.append(f""".. [#{fn_lbl}3] Methods with no analytic gradients do not have finite difference explicitly marked by "{trans_cell['fd'][args.sphinx]}", but the capability can be gleaned from the energy availability.""")
        elif args.mode == "details":
            lines.append(f""".. [#{fn_lbl}3] Finite difference gradients are only marked explicitly by "{trans_cell['fd'][args.sphinx]}" for overall (not per-method) lines and when at least one case has analytic gradients implemented, but the capability can be gleaned from the energy availability.""")
        else:
            lines.append(f""".. [#{fn_lbl}3] Finite difference gradients are not marked explicitly by "{trans_cell['fd'][args.sphinx]}", but the capability can be gleaned from the energy availability.""")
    if "h" in args.driver:
        if args.mode == "summary":
            lines.append(f""".. [#{fn_lbl}4] Methods with no analytic Hessians do not have finite difference explicitly marked by "{trans_cell['fd'][args.sphinx]}", but the capability can be gleaned from the gradient or energy availability.""")
        elif args.mode == "details":
            lines.append(f""".. [#{fn_lbl}4] Finite difference Hessians are only marked explicitly by "{trans_cell['fd'][args.sphinx]}" for overall (not per-method) lines and when at least one case has analytic Hessians implemented, but the capability can be gleaned from the gradient or energy availability.""")
        else:
            lines.append(f""".. [#{fn_lbl}4] Finite difference Hessians are not marked explicitly by "{trans_cell['fd'][args.sphinx]}", but the capability can be gleaned from the gradient or energy availability.""")
    for idx, note in enumerate(notes):
        lines.append(f".. [#{fn_lbl}{idx+10}] {note}")

    if not args.quiet:
        print("\n".join(lines))

    writefiles = args.writefile or [""]
    for wf in writefiles:
        if wf:
            tablefile = wf
        else:
            tablefile = f"autodoc_capabilities_{args.mode}.rst"

        tablefile = Path(tablefile)
        tablefile.write_text("\n".join(lines))


def pts(category, pyfile):
    print('Auto-documenting %s file %s' % (category, pyfile))


table_builder__ref_driver_type_fcae()
pts('capabilities', args.mode)

