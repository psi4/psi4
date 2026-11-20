#
#   For generating integral kernel C++ code
#

from multishells import MULTISHELLS, num_cartesian_integrals, get_hrr_blist
from dag import Graph, VRRBuilder, HRRBuilder
from codestrings import get_codestrings

import io
import textwrap
from pathlib import Path

# Number of random restarts ("kicks") of the recurrence path search after its local
# search has converged.  More kicks give slightly cheaper kernels and a slower codegen.
SEARCH_KICKS = 20

# Layout of the kernels
#
# The arrays of the interface are gb_size elements apart.  For a gb_size like 128 this is a
# divisor of the 4096 bytes spanned by one way of the L1 cache: the same gridpoint of all
# arrays then maps to the same few cache sets (and its loads and stores alias in their low
# 12 address bits), so that a kernel which touches many arrays per gridpoint keeps evicting
# its own data.  Such kernels are "blocked": the gridpoints are processed in blocks of a
# few cache lines, the integrals of a block are contracted in a contiguous buffer on the
# stack, and the buffer is copied to the output arrays one array after the other.  The
# copy costs a load and a store per integral, which does not pay off for small kernels;
# those work on the arrays in place.
#
# A kernel is blocked if it touches at least this many bytes per gridpoint (this depends
# on the precision, so that some kernels contain both variants) ...
BLOCK_MIN_BYTES = 280
# ... a block has this many gridpoints ...
BLOCK_POINTS = 32
# ... unless the buffers of the block would then exceed this many bytes.
BLOCK_MAX_BYTES = 24*1024

def cartesian_components(L):
    return [(lx, ly, L-lx-ly) for lx in range(L,-1,-1) for ly in range(L-lx,-1,-1)]

def ratio_expression(iLA, iLB, nLA, nLB):
    # the contraction coefficients of the highest L pair are part of the scaled Boys function
    if (iLA == nLA-1) and (iLB == nLB-1):
        return None
    return f"ratios[{iLA*nLB+iLB}]"

def build_vrr_graph(same_center, lqns, verbose=False):
    """Find a cheap VRR evaluation of the integrals lqns = [(a, b), ...]; returns the graph and its nodes."""
    graph = Graph()
    vrr = VRRBuilder(same_center, graph, kicks=SEARCH_KICKS)
    nodes = graph.factor(vrr.build([vrr.target(a, b) for (a, b) in lqns], verbose))
    if verbose:
        nmul, nadd = graph.count(nodes)
        print(f"  VRR: {nmul} multiplications and {nadd} additions per gridpoint and primitive.")
    return graph, nodes

def unique_rows(values):
    """Number the distinct values; returns (distinct values, row of every value)."""
    rows, row_of = [], []
    for v in values:
        if v not in rows:
            rows.append(v)
        row_of.append(rows.index(v))
    return rows, row_of

def write_grid_loop(f, graph, stores, indent, assign, npoints, stores_last):
    """Write one SIMD loop over gridpoints.

    stores: list of (output array name, graph node, ratio expression or None)"""
    op = "=" if assign else "+="
    order, nuses = graph.schedule([node for (_, node, _) in stores])
    outputs = {}
    for (lhs, node, ratio) in stores:
        outputs.setdefault(node, []).append((lhs, ratio))

    body, final = [], {}
    names = {}
    nvar = 0
    for node in order:
        if graph.defs[node][0] == 'leaf':
            names[node] = graph.defs[node][1]
            if node not in outputs:
                continue
        expr = graph.expr(node, names)
        # a value that only feeds a single output is written straight into the store
        if graph.defs[node][0] != 'leaf' and not (node in outputs and nuses[node] == 1):
            names[node] = f"x{nvar}"
            nvar += 1
            body.append(f"const T {names[node]} = {expr};")
            expr = names[node]
        for (lhs, ratio) in outputs[node] if node in outputs else []:
            if ratio is None:
                line = f"{lhs}[g] {op} {expr};"
            elif expr == names.get(node):
                line = f"{lhs}[g] {op} {ratio}*{expr};"
            else:
                line = f"{lhs}[g] {op} {ratio}*({expr});"
            if stores_last:
                final[lhs] = line
            else:
                body.append(line)

    f.write(f"{indent}#pragma omp simd\n")
    f.write(f"{indent}for (int g = 0; g < {npoints}; ++g) {{\n\n")
    for line in body:
        f.write(f"{indent}    {line}\n")
    if stores_last:
        f.write("\n")
        for (lhs, _, _) in stores:
            f.write(f"{indent}    {final[lhs]}\n")
    f.write(f"\n{indent}}} // for g\n")

def write_primitive_loop(f, graph, stores, Nmax, nratios, blocked):
    """Write the loop over primitives that evaluates and contracts the VRR.

    In place, all stores are at the end of the loop body: the compiler cannot prove that the
    outputs do not overlap the inputs, and reloads the inputs after every store.  A blocked
    kernel stores into its own buffer, as soon as a value is complete (fewer live values)."""
    npoints = "g_block" if blocked else "gb_size"
    # the loop bodies are generated first: they decide which scalars are needed
    body_assign, body_accumulate = io.StringIO(), io.StringIO()
    write_grid_loop(body_assign, graph, stores, " "*12, True, npoints, not blocked)
    write_grid_loop(body_accumulate, graph, stores, " "*12, False, npoints, not blocked)

    inputs = [("PGx", "v_PGx[prim*gb_size"), ("PGy", "v_PGy[prim*gb_size"), ("PGz", "v_PGz[prim*gb_size")]
    inputs += [(f"scaled_boys_{N}", f"v_scaled_boys[prim*{(Nmax+1)}*gb_size+{N}*gb_size") for N in range(0, Nmax+1)]

    f.write("\n    for (int prim = 0; prim < nprim; ++prim) {\n\n")

    if blocked:
        f.write("        if (padded) {\n")
        f.write("            for (int k = 0; k < gb_size; ++k) {\n")
        for i, (name, array) in enumerate(inputs):
            f.write(f"                pad[{i}*g_block+k] = {array}+k];\n")
        f.write("            }\n")
        f.write("        }\n")
        for i, (name, array) in enumerate(inputs):
            f.write(f"        const T* __restrict {name} = padded ? &pad[{i}*g_block] : &{array} + g_lo];\n")
    else:
        for (name, array) in inputs:
            f.write(f"        const T* __restrict {name} = &{array}];\n")
    if nratios > 0:
        f.write(f"        const T* __restrict ratios = &v_ratios[prim*{nratios}];\n")
    f.write("\n")

    for line in graph.scalar_definitions():
        f.write(f"        {line}\n")

    f.write("\n        // overwrite integrals for first primitive, accumulate otherwise\n")
    f.write("        if (prim==0) {\n")
    f.write(body_assign.getvalue())
    f.write("        } else {\n\n")
    f.write(body_accumulate.getvalue())
    f.write("        } // if/else prim\n")
    f.write("    } // for prim\n")

def block_points(nrows, size):
    """Gridpoints per block for buffers of nrows arrays of a type with the given size."""
    points = BLOCK_POINTS
    while points*size*nrows > BLOCK_MAX_BYTES and points*size > 64:
        points //= 2
    return points

def write_block_buffers(f, buffers, Nmax):
    """Declare the buffers of a blocked kernel; buffers = [(name, number of rows), ...]."""
    nrows = sum(n for (_, n) in buffers)
    f.write("    if (nprim < 1) return;\n\n")
    f.write("    // The gridpoints are processed in blocks.  The integrals of a block are contracted in\n")
    f.write("    // a buffer that stays in the L1 cache and are then copied to the output arrays.\n")
    f.write(f"    enum {{ g_block = (sizeof(T) <= 4) ? {block_points(nrows, 4)} : {block_points(nrows, 8)} }};\n")
    for (name, n) in buffers:
        f.write(f"    alignas(64) T {name}[{n}*g_block];\n")
        for i in range(0, n):
            f.write(f"    T* __restrict {name}_{i} = &{name}[{i}*g_block];\n")
    f.write("\n    // fewer gridpoints than one block: work on a zero padded copy of the inputs\n")
    f.write("    const bool padded = (gb_size < g_block);\n")
    f.write(f"    alignas(64) T pad[{Nmax+4}*g_block];\n")
    f.write("    if (padded) {\n")
    f.write(f"        for (int k = 0; k < {Nmax+4}*g_block; ++k) pad[k] = 0;\n")
    f.write("    }\n")

def write_block_loop(f, body):
    f.write("\n    for (int g_lo = 0; g_lo < gb_size; g_lo += g_block) {\n\n")
    f.write("        // all blocks are complete: the last one overlaps its predecessor\n")
    f.write("        if (g_lo + g_block > gb_size && !padded) g_lo = gb_size - g_block;\n")
    f.write(textwrap.indent(body, "    "))
    f.write("    } // for g_lo\n")

def write_copy_out(f, source, row_of):
    """Copy the rows of a block buffer to the output arrays, one array after the other."""
    f.write("\n    // copy the block to the output arrays\n")
    row = "i"
    if row_of != list(range(len(row_of))):
        f.write(f"    static const short row[{len(row_of)}] = {{")
        for i in range(0, len(row_of), 20):
            f.write("\n        " + " ".join(f"{r}," for r in row_of[i:i+20]))
        f.write("\n    };\n")
        row = "row[i]"
    f.write(f"    for (int i = 0; i < {len(row_of)}; ++i) {{\n")
    f.write(f"        const T* __restrict src = &{source}[{row}*g_block];\n")
    f.write("        T* __restrict dst = &integrals[i*gb_size + g_lo];\n")
    f.write("        if (padded) {\n")
    f.write("            SNK_COPY(dst, src, gb_size);\n")
    f.write("        } else {\n")
    f.write("            SNK_COPY(dst, src, g_block);\n")
    f.write("        }\n")
    f.write("    }\n")

# Compilers do not reliably turn a short copy loop into vector moves; the builtin needs no header.
COPY_MACRO = """#ifndef SNK_COPY
#if defined(__GNUC__) || defined(__clang__)
#define SNK_COPY(dst, src, n) __builtin_memcpy(dst, src, (n)*sizeof(*(dst)))
#else
#define SNK_COPY(dst, src, n) for (int k_ = 0; k_ < (n); ++k_) (dst)[k_] = (src)[k_]
#endif
#endif

"""

def write_kernel(dname, func_name, hrr, need_coeff_ratios, same_center, narrays, body):
    """Write header and implementation; body(blocked) returns the function body of a layout."""
    args_str_T, args_str_f, args_str_d, header_str = get_codestrings(hrr, need_coeff_ratios, same_center)
    func_proto_h = f"void {func_name}(\n{header_str}"
    func_proto_T = f"void {func_name}(\n{args_str_T}"
    func_proto_f = f"void {func_name}<float>(\n{args_str_f}"
    func_proto_d = f"void {func_name}<double>(\n{args_str_d}"

    # which layouts are needed for single and double precision?
    blocked = {narrays*size >= BLOCK_MIN_BYTES for size in (4, 8)}
    if len(blocked) == 1:
        code = body(blocked.pop())
    else:
        code = f"    if ({narrays}*sizeof(T) >= {BLOCK_MIN_BYTES}) {{\n\n"
        code += textwrap.indent(body(True), "    ")
        code += "\n    } else {\n\n"
        code += textwrap.indent(body(False), "    ")
        code += "\n    }\n"

    code_dir = Path(dname)
    code_dir.mkdir(exist_ok=True)

    # Header file
    with open(f"{dname}/{func_name}.h", 'w') as f:
        f.write(f"#ifndef SNK_{func_name}\n")
        f.write(f"#define SNK_{func_name}\n\n")
        f.write("template <typename T>\n")
        f.write(f"{func_proto_h}\n\n")
        f.write(f"#endif // SNK_{func_name}")

    # Implementation file
    with open(f"{dname}/{func_name}.C", 'w') as f:
        f.write(f'#include "{func_name}.h"\n\n')
        if "SNK_COPY" in code:
            f.write(COPY_MACRO)
        f.write("template <typename T>\n")
        f.write(f"{func_proto_T} \n{{\n")
        f.write(code)
        f.write(f"}} // {func_name}\n\n")
        f.write(f"template {func_proto_f};\n\n")
        f.write(f"template {func_proto_d};\n\n")

#####################################
#                                   #
#   VRR Automatic Code Generation   #
#                                   #
#####################################

def generate_vrr(MSA, MSB, same_center, dname="SNK_CODEGEN", verbose=False):
    vLA, vLB = MULTISHELLS[MSA], MULTISHELLS[MSB]
    nci = num_cartesian_integrals(vLA, vLB)

    nLA, nLB = len(vLA), len(vLB)
    Nmax = max(vLA) + max(vLB)

    if verbose:
        sc_str = "same center," if same_center else ""
        print(f"Automatic integral code generation, {MSA}-{MSB} (multi)shellpair,{sc_str} {nci} cartesian integrals.")
        print(f"Searching for the VRR path...")

    lqns, ratios = [], []
    for iLA, LA in enumerate(vLA):
        for a in cartesian_components(LA):
            for iLB, LB in enumerate(vLB):
                for b in cartesian_components(LB):
                    lqns.append((a, b))
                    ratios.append(ratio_expression(iLA, iLB, nLA, nLB))
    if (len(lqns) != nci):
        raise RuntimeError(f"codegen.generate_vrr(): Expected {nci} integrals, made {len(lqns)} integrals!")
    graph, nodes = build_vrr_graph(same_center, lqns, verbose)
    values = list(zip(nodes, ratios))

    if verbose:
        print(f"Writing generated code...")

    def body(blocked):
        f = io.StringIO()
        if blocked:
            # same center integrals are not all different: only the distinct ones are contracted
            rows, row_of = unique_rows(values)
            write_block_buffers(f, [("block", len(rows))], Nmax)
            stores = [(f"block_{i}", node, ratio) for i, (node, ratio) in enumerate(rows)]
            loops = io.StringIO()
            write_primitive_loop(loops, graph, stores, Nmax, nLA*nLB-1, True)
            write_copy_out(loops, "block", row_of)
            write_block_loop(f, loops.getvalue())
        else:
            for i in range(0, nci):
                f.write(f"    T* __restrict integral_{i} = &integrals[{i}*gb_size];\n")
            stores = [(f"integral_{i}", node, ratio) for i, (node, ratio) in enumerate(values)]
            write_primitive_loop(f, graph, stores, Nmax, nLA*nLB-1, False)
        return f.getvalue()

    func_name = f"VRR_{MSA}_{MSB}"
    if same_center:
        func_name += "_SAMECENTER"
    # arrays per gridpoint: the integrals, the Boys functions and PG of one primitive
    write_kernel(dname, func_name, False, nLA*nLB > 1, same_center, nci + Nmax+1 + 3, body)

    if verbose:
        print(f"Finished!")

#####################################
#                                   #
#   HRR Automatic Code Generation   #
#                                   #
#####################################

def generate_hrr(MSA, MSB, same_center, dname="SNK_CODEGEN", verbose=False):
    vLA, vLB = MULTISHELLS[MSA], MULTISHELLS[MSB]
    nci = num_cartesian_integrals(vLA, vLB)
    hrr_blist = get_hrr_blist(vLA, vLB, same_center)

    nLA, nLB = len(vLA), len(vLB)
    LAmax, LBmax = max(vLA), max(vLB)
    Nmax = LAmax + LBmax
    if (LAmax==0) or (LBmax==0):
        if verbose:
            print(f"Skipping this {MSA}-{MSB} (multi)shellpair, no HRRs to do!")
        return

    if verbose:
        sc_str = "same center," if same_center else ""
        print(f"Automatic integral code generation, {MSA}-{MSB} (multi)shellpair,{sc_str} {nci} cartesian integrals.")
        print(f"Searching for the VRR path ({len(hrr_blist)} base integrals)...")

    #
    #   Base integrals
    #

    bi_lqns, bi_ratios = [], []
    for lqns in hrr_blist:   # insertion order is the base integral index
        bi_lqns.append((lqns[0:3], lqns[3:6]))
        bi_ratios.append(ratio_expression(vLA.index(lqns[6]), vLB.index(lqns[7]), nLA, nLB))
    bi_graph, bi_nodes = build_vrr_graph(same_center, bi_lqns, verbose)
    bi_stores = [(f"base_integral_{i}", node, ratio) for i, (node, ratio) in enumerate(zip(bi_nodes, bi_ratios))]

    #
    #   Final HRR integrals
    #

    if verbose:
        print(f"Searching for the HRR path ({nci} final integrals)...")

    graph = Graph()
    hrr = HRRBuilder(same_center, graph, hrr_blist)
    targets = []
    for LA in vLA:
        for a in cartesian_components(LA):
            for LB in vLB:
                for b in cartesian_components(LB):
                    targets.append(a + b + (LA, LB))
    if (len(targets) != nci):
        raise RuntimeError(f"codegen.generate_hrr(): Expected {nci} integrals, made {len(targets)} integrals!")
    nodes = hrr.build(targets, verbose)

    if verbose:
        nmul, nadd = graph.count(nodes)
        print(f"  HRR: {nmul} multiplications and {nadd} additions per gridpoint.")
        print(f"Writing generated code...")

    def body(blocked):
        f = io.StringIO()
        if blocked:
            # The base integrals of a block live in a local buffer: base_integrals is not used.
            # Same center integrals are base integrals: there is nothing to compute in the HRR.
            copies = all(graph.defs[node][0] == 'leaf' for node in nodes)
            rows, row_of = unique_rows(nodes)
            write_block_buffers(f, [("base_integral", len(bi_stores))] + ([] if copies else [("block", len(rows))]), Nmax)
            loops = io.StringIO()
            write_primitive_loop(loops, bi_graph, bi_stores, Nmax, nLA*nLB-1, True)
            if copies:
                write_copy_out(loops, "base_integral", [hrr.base_index(node) for node in nodes])
            else:
                stores = [(f"block_{i}", node, None) for i, node in enumerate(rows)]
                loops.write("\n")
                write_grid_loop(loops, graph, stores, " "*4, True, "g_block", False)
                write_copy_out(loops, "block", row_of)
            write_block_loop(f, loops.getvalue())
        else:
            for i in range(0, len(bi_stores)):
                f.write(f"    T* __restrict base_integral_{i} = &base_integrals[{i}*gb_size];\n")
            write_primitive_loop(f, bi_graph, bi_stores, Nmax, nLA*nLB-1, False)
            f.write("\n")
            for i in range(0, nci):
                f.write(f"    T* __restrict integral_{i} = &integrals[{i}*gb_size];\n")
            f.write("\n")
            stores = [(f"integral_{i}", node, None) for i, node in enumerate(nodes)]
            write_grid_loop(f, graph, stores, " "*4, True, "gb_size", True)
        return f.getvalue()

    func_name = f"HRR_{MSA}_{MSB}"
    if same_center:
        func_name += "_SAMECENTER"
    # arrays per gridpoint: base integrals, integrals, the Boys functions and PG of one primitive
    write_kernel(dname, func_name, True, nLA*nLB > 1, same_center, len(bi_stores) + nci + Nmax+1 + 3, body)

    if verbose:
        print(f"Finished!")
