#
#   Exact verification of generated kernels
#
#   The loop bodies of a generated .C file are executed statement by statement in
#   exact rational arithmetic for random inputs and compared with the recurrence
#   relations of recurrences.py, which are unrolled symbolically by SymPy without
#   any of the path optimizations of dag.py.  Any difference is an error in the
#   generated code, not roundoff.
#
#   Usage: python verify.py            (a default selection of kernels)
#          python verify.py vrr:D:P:0 hrr:G:F:0 vrr:F:F:1 ...   (kind:MSA:MSB:same_center)
#

import codegen
from multishells import MULTISHELLS, get_hrr_blist
from recurrences import do_VRR_multishell, do_HRR_multishell

import random
import re
import sys
import tempfile
from fractions import Fraction
from sympy import Rational

def random_fraction(rng):
    return Fraction(rng.randint(-40, 40) or 7, rng.randint(7, 40))

def parse_loops(text):
    """Returns (scalar statements, loops, copy); a loop is a list of (lhs, accumulate, rhs).

    copy = (buffer name, row of every output) describes the copy of a block buffer to the
    output arrays of a blocked kernel, and is None for kernels that work in place."""
    scalars, loops, cur = [], [], None
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith("for (int g ="):
            cur = []
        elif line == "} // for g":
            loops.append(cur)
            cur = None
        elif cur is not None and "=" in line:
            acc = "+=" in line
            lhs, rhs = line.rstrip(";").split("+=" if acc else "=", 1)
            cur.append((lhs.replace("const T", "").strip(), acc, rhs.strip()))
        elif cur is None:
            m = re.match(r"const T (\w+) = ([^;]*);", line)
            if m and "v_" not in m.group(2):
                scalars.append((m.group(1), m.group(2)))
    copy = None
    m = re.search(r"const T\* __restrict src = &(\w+)\[(row\[i\]|i)\*", text)
    if m:
        n = int(re.search(r"for \(int i = 0; i < (\d+); \+\+i\)", text).group(1))
        rows = list(range(n))
        if m.group(2) != "i":
            table = re.search(r"static const short row\[\d+\] = \{([^}]*)\}", text).group(1)
            rows = [int(r) for r in table.replace("\n", " ").split(",") if r.strip()]
        copy = (m.group(1), rows)
    return scalars, loops, copy

def layouts(path):
    """The function body of a kernel, once per layout that it contains."""
    text = open(path).read()
    if "*sizeof(T) >= " not in text:
        return [text]
    # blocked and in place variant, selected by the size of T
    head, rest = text.split("*sizeof(T) >= ", 1)
    blocked, in_place = rest.split("\n    } else {\n", 1)
    return [blocked, in_place]

def pythonize(s):
    s = re.sub(r"\[g\]", "", s)
    return re.sub(r"\[(\d+)\]", r"_\1", s)

def run_loop(stmts, env, out):
    """Execute one grid loop body; out collects the stored arrays (accumulating into them)."""
    env = dict(env)
    for (lhs, acc, rhs) in stmts:
        val = eval(pythonize(rhs), {"__builtins__": {}}, env)
        if "[g]" in lhs:
            name = pythonize(lhs)
            out[name] = (out.get(name, 0) + val) if acc else val
        else:
            env[lhs] = val
    return out

def sympy_value(expr, env):
    return expr.xreplace({s: Rational(env[pythonize(s.name)].numerator, env[pythonize(s.name)].denominator)
                          for s in expr.free_symbols})

def verify(kind, MSA, MSB, same_center, nprim=2, seed=1, dname=None, verbose=True):
    vLA, vLB = MULTISHELLS[MSA], MULTISHELLS[MSB]
    if kind == "hrr" and (max(vLA) == 0 or max(vLB) == 0):
        return True
    name = f"{kind.upper()}_{MSA}_{MSB}" + ("_SAMECENTER" if same_center else "")

    with tempfile.TemporaryDirectory() as tmp:
        d = dname or tmp
        getattr(codegen, f"generate_{kind}")(MSA=MSA, MSB=MSB, same_center=same_center, dname=d, verbose=False)
        texts = layouts(f"{d}/{name}.C")

    # reference expressions
    if kind == "vrr":
        ref_prim, ref_hrr = do_VRR_multishell(vLA, vLB, same_center), None
        prim_ratio = [(iLA, iLB) for iLA, LA in enumerate(vLA) for _ in range((LA+1)*(LA+2)//2)
                      for iLB, LB in enumerate(vLB) for _ in range((LB+1)*(LB+2)//2)]
    else:
        blist = get_hrr_blist(vLA, vLB, same_center)
        ref_prim, ref_hrr = do_HRR_multishell(vLA, vLB, same_center, blist)
        prim_ratio = [(vLA.index(k[6]), vLB.index(k[7])) for k in blist]

    nbad = ntot = 0
    for text in texts:
        bad, tot = verify_layout(text, kind, vLA, vLB, ref_prim, ref_hrr, prim_ratio, nprim, seed)
        nbad, ntot = nbad + bad, ntot + tot

    if verbose:
        note = " (blocked and in place)" if len(texts) > 1 else ""
        print(f"{name:28s} {ntot:5d} integrals checked exactly{note}: {'OK' if nbad == 0 else f'{nbad} WRONG'}")
    return nbad == 0

def verify_layout(text, kind, vLA, vLB, ref_prim, ref_hrr, prim_ratio, nprim, seed):
    """Returns (number of wrong integrals, number of integrals) of one layout of a kernel."""
    nLA, nLB = len(vLA), len(vLB)
    Nmax = max(vLA) + max(vLB)
    scalars, loops, copy = parse_loops(text)

    rng = random.Random(seed)
    AB = {f"AB{c}": random_fraction(rng) for c in "xyz"}

    # primitive loop: the first primitive assigns (loop 0), the others accumulate (loop 1)
    out, expected = {}, [0]*len(ref_prim)
    for prim in range(nprim):
        env = dict(AB)
        env["tpi"] = random_fraction(rng)
        for c in "xyz":
            env[f"PA{c}"] = random_fraction(rng)
            env[f"PB{c}"] = env[f"PA{c}"] + AB[f"AB{c}"]
            env[f"PG{c}"] = random_fraction(rng)
        for N in range(Nmax+1):
            env[f"scaled_boys_{N}"] = random_fraction(rng)
        for r in range(nLA*nLB-1):
            env[f"ratios_{r}"] = random_fraction(rng)
        for (sname, sexpr) in scalars:
            if sname not in env:
                env[sname] = eval(pythonize(sexpr), {"__builtins__": {}}, env)
        run_loop(loops[0 if prim == 0 else 1], env, out)
        for i, e in enumerate(ref_prim):
            iLA, iLB = prim_ratio[i]
            ratio = 1 if (iLA == nLA-1 and iLB == nLB-1) else env[f"ratios_{iLA*nLB+iLB}"]
            expected[i] += ratio * sympy_value(e, env)

    def result(values, i):
        # integral i: in its own array, or in a row of the block buffer that is copied out
        if copy is None:
            return values.get(f"integral_{i}")
        return values.get(f"{copy[0]}_{copy[1][i]}")

    if kind == "vrr":
        return sum(1 for i, e in enumerate(expected) if result(out, i) != e), len(expected)

    # the base integrals are always computed into base_integral_*
    nbad = sum(1 for i, e in enumerate(expected) if out.get(f"base_integral_{i}") != e)
    env = dict(AB)
    env.update(out)
    final = run_loop(loops[2], env, {}) if len(loops) > 2 else {}
    final.update(out)
    nbad += sum(1 for i, e in enumerate(ref_hrr) if result(final, i) != sympy_value(e, env))
    return nbad, len(expected) + len(ref_hrr)

DEFAULT = ["vrr:S:S:0", "vrr:P:S:0", "vrr:S:P:0", "vrr:P:P:0", "vrr:D:P:0", "vrr:P:D:0", "vrr:D:D:0", "vrr:F:D:0",
           "vrr:P:P:1", "vrr:D:P:1", "vrr:D:D:1", "vrr:F:F:1",
           "hrr:P:P:0", "hrr:D:P:0", "hrr:P:D:0", "hrr:D:D:0", "hrr:F:D:0", "hrr:D:F:0",
           "hrr:P:P:1", "hrr:D:P:1", "hrr:P:D:1", "hrr:D:D:1", "hrr:F:F:1"]

if __name__ == "__main__":
    ok = True
    for spec in (sys.argv[1:] or DEFAULT):
        kind, MSA, MSB, sc = spec.split(":")
        ok = verify(kind, MSA, MSB, bool(int(sc))) and ok
    sys.exit(0 if ok else 1)
