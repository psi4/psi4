#
#   For generating integral kernel C++ code
#

from multishells import MULTISHELLS, num_cartesian_integrals
from recurrences import do_VRR_multishell, get_hrr_blist, do_HRR_multishell
from codestrings import get_codestrings

from pathlib import Path
from sympy import cse
from sympy.printing.c import C99CodePrinter

class CustomCPrinter(C99CodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        if exp.is_Integer and exp.is_positive:
            base_str = self._print(base)
            return "*".join([base_str] * int(exp))
        return super()._print_Pow(expr)

def ccode_nopow(expr, **settings):
    return CustomCPrinter(settings).doprint(expr)

#####################################
#                                   #
#   VRR Automatic Code Generation   #
#                                   #
#####################################

def generate_vrr(MSA, MSB, same_center, dname="SNK_CODEGEN", verbose=False):
    vLA, vLB = MULTISHELLS[MSA], MULTISHELLS[MSB]
    nci = num_cartesian_integrals(vLA, vLB)

    nLA, nLB = len(vLA), len(vLB)
    LAmax, LBmax = max(vLA), max(vLB)
    Nmax = LAmax + LBmax
    
    if verbose:
        sc_str = "same center," if same_center else ""
        print(f"Automatic integral code generation, {MSA}-{MSB} (multi)shellpair,{sc_str} {nci} cartesian integrals.")

    vrr_exprs = do_VRR_multishell(vLA, vLB, same_center, verbose)

    if verbose:
        print(f"Performing CSE...")

    cse_repls, cse_exprs = cse(vrr_exprs, optimizations='basic')
    
    if verbose:
        print(f"Writing generated code...")

    need_coeff_ratios = (len(vLA)>1) or (len(vLB)>1)
    args_str_T, args_str_f, args_str_d, header_str = get_codestrings(False, need_coeff_ratios, same_center)
    
    func_name = f"VRR_{MSA}_{MSB}"
    if same_center:
        func_name += "_SAMECENTER"
    func_proto_h = f"void {func_name}(\n{header_str}"
    func_proto_T = f"void {func_name}(\n{args_str_T}"
    func_proto_f = f"void {func_name}<float>(\n{args_str_f}"
    func_proto_d = f"void {func_name}<double>(\n{args_str_d}"

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
        f.write("template <typename T>\n")
        f.write(f"{func_proto_T} \n{{\n")

        for i in range(0, len(cse_exprs)):
            f.write(f"    T* __restrict integral_{i} = &integrals[{i}*gb_size];\n")
         
        f.write("\n    for (int prim = 0; prim < nprim; ++prim) {\n\n")
            
        f.write("        const T tpi = 0.5 / v_p[prim];\n") 
        if not same_center:
            if (LAmax > 0):
                f.write("        const T PAx = v_PAx[prim];\n")
                f.write("        const T PAy = v_PAy[prim];\n")
                f.write("        const T PAz = v_PAz[prim];\n")
            if (LBmax > 0):
                f.write("        const T PBx = PAx + ABx;\n")
                f.write("        const T PBy = PAy + ABy;\n")
                f.write("        const T PBz = PAz + ABz;\n") 
        
        f.write("        const T* __restrict PGx = &v_PGx[prim*gb_size];\n")
        f.write("        const T* __restrict PGy = &v_PGy[prim*gb_size];\n")
        f.write("        const T* __restrict PGz = &v_PGz[prim*gb_size];\n")
        if need_coeff_ratios:
            f.write(f"        const T* __restrict ratios = &v_ratios[prim*{nLA*nLB-1}];\n\n")
        else:
            f.write("\n")

        for N in range(0, Nmax+1):
            f.write(f"        const T* __restrict scaled_boys_{N} = &v_scaled_boys[prim*{(Nmax+1)}*gb_size+{N}*gb_size];\n")
        
        f.write("\n        // overwrite integrals for first primitive, accumulate otherwise\n")
        f.write("        if (prim==0) {\n")

        f.write(f"            #pragma omp simd\n")
        f.write(f"            for (int g = 0; g < gb_size; ++g) {{\n\n")
        for r in cse_repls:
            f.write(f"                const T {r[0]} = {ccode_nopow(r[1])};\n")
        if (len(cse_repls) > 0):
            f.write("\n")

        if need_coeff_ratios:
            integral_count = 0
            for iLA in range(0, nLA):
                LA = vLA[iLA]
                for lax in range(LA,-1,-1):
                    for lay in range(LA-lax,-1,-1):
                        for iLB in range(0, nLB):
                            LB = vLB[iLB]
                            for lbx in range(LB,-1,-1):
                                for lby in range(LB-lbx,-1,-1):
                                    cse_expr = ccode_nopow(cse_exprs[integral_count]);
                                    if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                        f.write(f"                integral_{integral_count}[g] = {cse_expr};\n")
                                    else:
                                        f.write(f"                integral_{integral_count}[g] = ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                    integral_count += 1
        else: 
            for i, e in enumerate(cse_exprs):
                f.write(f"                integral_{i}[g] = {ccode_nopow(e)};\n")

        f.write("\n            } // for g\n")
        f.write("        } else {\n")

        f.write(f"\n            #pragma omp simd\n")
        f.write(f"            for (int g = 0; g < gb_size; ++g) {{\n\n")
        for r in cse_repls:
            f.write(f"                const T {r[0]} = {ccode_nopow(r[1])};\n")
        if (len(cse_repls) > 0):
            f.write("\n")

        if need_coeff_ratios:
            integral_count = 0
            for iLA in range(0, nLA):
                LA = vLA[iLA]
                for lax in range(LA,-1,-1):
                    for lay in range(LA-lax,-1,-1):
                        for iLB in range(0, nLB):
                            LB = vLB[iLB]
                            for lbx in range(LB,-1,-1):
                                for lby in range(LB-lbx,-1,-1):
                                    cse_expr = ccode_nopow(cse_exprs[integral_count]);
                                    if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                        f.write(f"                integral_{integral_count}[g] += {cse_expr};\n")
                                    else:
                                        f.write(f"                integral_{integral_count}[g] += ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                    integral_count += 1
        else: 
            for i, e in enumerate(cse_exprs):
                f.write(f"                integral_{i}[g] += {ccode_nopow(e)};\n")

        f.write("\n            } // for g\n")

        f.write("        } // if/else prim\n")

        f.write("    } // for prim\n")
        f.write(f"}} // {func_name}\n\n")
        f.write(f"template {func_proto_f};\n\n")
        f.write(f"template {func_proto_d};\n\n")

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
    if (LAmax==0) or (LBmax==0):
        if verbose:
            print(f"Skipping this {MSA}-{MSB} (multi)shellpair, no HRRs to do!")
        return

    Nmax = LAmax + LBmax

    if verbose:
        sc_str = "same center," if same_center else ""
        print(f"Automatic integral code generation, {MSA}-{MSB} (multi)shellpair,{sc_str} {nci} cartesian integrals.")

    hrr_bi_exprs, hrr_exprs = do_HRR_multishell(vLA, vLB, same_center, hrr_blist, verbose)

    if verbose:
        print(f"Performing CSE...")

    cse_bi_repls, cse_bi_exprs = cse(hrr_bi_exprs, optimizations='basic')
    cse_repls, cse_exprs = cse(hrr_exprs, optimizations='basic')

    need_coeff_ratios = (len(vLA)>1) or (len(vLB)>1)
    args_str_T, args_str_f, args_str_d, header_str = get_codestrings(True, need_coeff_ratios, same_center)
   
    code_dir = Path(dname)
    code_dir.mkdir(exist_ok=True)
     
    if verbose:
        print(f"Writing generated code...")

    func_name = f"HRR_{MSA}_{MSB}"
    if same_center:
        func_name += "_SAMECENTER"
    func_proto_h = f"void {func_name}(\n{header_str}"
    func_proto_T = f"void {func_name}(\n{args_str_T}"
    func_proto_f = f"void {func_name}<float>(\n{args_str_f}"
    func_proto_d = f"void {func_name}<double>(\n{args_str_d}"
    
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
        f.write("template <typename T>\n")
        f.write(f"{func_proto_T} \n{{\n")

        #
        #   Base integrals
        #
        
        for i in range(0, len(cse_bi_exprs)):
            f.write(f"    T* __restrict base_integral_{i} = &base_integrals[{i}*gb_size];\n")
         
        f.write("\n    for (int prim = 0; prim < nprim; ++prim) {\n\n")
            
        f.write("        const T tpi = 0.5 / v_p[prim];\n") 
        if not same_center:
            if (LAmax > 0):
                f.write("        const T PAx = v_PAx[prim];\n")
                f.write("        const T PAy = v_PAy[prim];\n")
                f.write("        const T PAz = v_PAz[prim];\n")
            if (LBmax > 0):
                f.write("        const T PBx = PAx + ABx;\n")
                f.write("        const T PBy = PAy + ABy;\n")
                f.write("        const T PBz = PAz + ABz;\n") 
        
        f.write("        const T* __restrict PGx = &v_PGx[prim*gb_size];\n")
        f.write("        const T* __restrict PGy = &v_PGy[prim*gb_size];\n")
        f.write("        const T* __restrict PGz = &v_PGz[prim*gb_size];\n")
        if need_coeff_ratios:
            f.write(f"        const T* __restrict ratios = &v_ratios[prim*{nLA*nLB-1}];\n\n")
        else:
            f.write("\n")

        for N in range(0, Nmax+1):
            f.write(f"        const T* __restrict scaled_boys_{N} = &v_scaled_boys[prim*{(Nmax+1)}*gb_size+{N}*gb_size];\n")
        
        f.write("\n        // overwrite integrals for first primitive, accumulate otherwise\n")
        f.write("        if (prim==0) {\n")
        
        f.write(f"            #pragma omp simd\n")
        f.write(f"            for (int g = 0; g < gb_size; ++g) {{\n\n")
        for r in cse_bi_repls:
            f.write(f"                const T {r[0]} = {ccode_nopow(r[1])};\n")
        if (len(cse_bi_repls) > 0):
            f.write("\n")
        
        if need_coeff_ratios:
            integral_count = 0
            for iLA in range(0, nLA):
                LA = vLA[iLA]
                for iLB in range(0, nLB):
                    LB = vLB[iLB]
                    if same_center:
                        L = LA + LB
                        for lx in range(L,-1,-1):
                            for ly in range(L-lx,-1,-1):
                                cse_expr = ccode_nopow(cse_bi_exprs[integral_count]);
                                if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                    f.write(f"                base_integral_{integral_count}[g] = {cse_expr};\n")
                                else:
                                    f.write(f"                base_integral_{integral_count}[g] = ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                integral_count += 1
                    else:
                        for L in range(max(LA,LB), LA+LB+1):
                            for lx in range(L,-1,-1):
                                for ly in range(L-lx,-1,-1):
                                    cse_expr = ccode_nopow(cse_bi_exprs[integral_count]);
                                    if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                        f.write(f"                base_integral_{integral_count}[g] = {cse_expr};\n")
                                    else:
                                        f.write(f"                base_integral_{integral_count}[g] = ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                    integral_count += 1
        else: 
            for i, e in enumerate(cse_bi_exprs):
                f.write(f"                base_integral_{i}[g] = {ccode_nopow(e)};\n")

        f.write("\n            } // for g\n")
        f.write("        } else {\n")
        
        f.write(f"\n            #pragma omp simd\n")
        f.write(f"            for (int g = 0; g < gb_size; ++g) {{\n\n")
        for r in cse_bi_repls:
            f.write(f"                const T {r[0]} = {ccode_nopow(r[1])};\n")
        if (len(cse_bi_repls) > 0):
            f.write("\n")
        
        if need_coeff_ratios:
            integral_count = 0
            for iLA in range(0, nLA):
                LA = vLA[iLA]
                for iLB in range(0, nLB):
                    LB = vLB[iLB]
                    if same_center: 
                        L = LA + LB
                        for lx in range(L,-1,-1):
                            for ly in range(L-lx,-1,-1):
                                cse_expr = ccode_nopow(cse_bi_exprs[integral_count]);
                                if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                    f.write(f"                base_integral_{integral_count}[g] += {cse_expr};\n")
                                else:
                                    f.write(f"                base_integral_{integral_count}[g] += ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                integral_count += 1
                    else:
                        for L in range(max(LA,LB), LA+LB+1):
                            for lx in range(L,-1,-1):
                                for ly in range(L-lx,-1,-1):
                                    cse_expr = ccode_nopow(cse_bi_exprs[integral_count]);
                                    if (iLA==(nLA-1)) and (iLB==(nLB-1)):
                                        f.write(f"                base_integral_{integral_count}[g] += {cse_expr};\n")
                                    else:
                                        f.write(f"                base_integral_{integral_count}[g] += ratios[{iLA*nLB+iLB}] * ({cse_expr});\n")
                                    integral_count += 1
        else: 
            for i, e in enumerate(cse_bi_exprs):
                f.write(f"                base_integral_{i}[g] += {ccode_nopow(e)};\n")
        
        f.write("\n            } // for g\n")
        f.write("        } // if/else prim\n")
        
        f.write("    } // for prim\n\n")

        #
        #   Final HRR integrals
        #
        
        for i in range(0, len(cse_exprs)):
            f.write(f"    T* __restrict integral_{i} = &integrals[{i}*gb_size];\n") 
        f.write("\n    #pragma omp simd\n")
        f.write("    for (int g = 0; g < gb_size; ++g) {\n\n")
        for r in cse_repls:
            f.write(f"        const T {r[0]} = {ccode_nopow(r[1])};\n")
        if (len(cse_repls) > 0):
            f.write("\n")

        for i, e in enumerate(cse_exprs):
            f.write(f"        integral_{i}[g] = {ccode_nopow(e)};\n")
        f.write("\n    } // for g\n")

        f.write(f"}} // {func_name}\n\n")
        f.write(f"template {func_proto_f};\n\n")
        f.write(f"template {func_proto_d};\n\n")
    
    if verbose:
        print(f"Finished!")
