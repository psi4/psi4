#!/usr/bin/env python

"""
isort:skip_file
"""
import collections
import importlib
import math
import os
import subprocess
import sys
import sysconfig
import time

if sys.version_info <= (3, 0):
    print('Much of this script needs py3')
    sys.exit()

# good
#import numpy as np
#import psi4

# bad
import psi4
import numpy as np


def test_threaded_blas(args):
    threads = int(args.nthread)

    times = {}

    size = [200, 500, 2000, 4000]
    threads = [1, threads]

    for th in threads:
        psi4.set_num_threads(th)

        for sz in size:
            nruns = max(1, int(1.e10 / (sz ** 3)))

            a = psi4.core.Matrix(sz, sz)
            b = psi4.core.Matrix(sz, sz)
            c = psi4.core.Matrix(sz, sz)

            tp4 = time.time()
            for n in range(nruns):
                c.gemm(False, False, 1.0, a, b, 0.0)

            retp4 = (time.time() - tp4) / nruns

            tnp = time.time()
            for n in range(nruns):
                np.dot(a, b, out=np.asarray(c))

            retnp = (time.time() - tnp) / nruns
            print("Time for threads %2d, size %5d: Psi4: %12.6f  NumPy: %12.6f" % (th, sz, retp4, retnp))
            if sz == 4000:
                times[f"p4-n{th}"] = retp4
                times[f"np-n{th}"] = retnp
                assert psi4.get_num_threads() == th

    rat1 = times["np-n" + str(threads[-1])] / times["p4-n" + str(threads[-1])]
    rat2 = times["p4-n" + str(threads[0])] / times["p4-n" + str(threads[-1])]
    print("  NumPy@n%d : Psi4@n%d ratio (want ~1): %.2f" % (threads[-1], threads[-1], rat1))
    print("   Psi4@n%d : Psi4@n%d ratio (want ~%d): %.2f" % (threads[0], threads[-1], threads[-1], rat2))
    if args.passfail:
        assert math.isclose(rat1, 1.0, rel_tol=0.4), f'PsiAPI:NumPy speedup {rat1} !~= 1.0'
        assert math.isclose(rat2, threads[-1], rel_tol=0.4), f'PsiAPI speedup {rat2} !~= {threads[-1]}'


def test_psithon(args):

    inputdat = """
molecule dimer {
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
}

set {
    BASIS jun-cc-pVQZ
    SCF_TYPE DF
    FREEZE_CORE True
    DF_BASIS_ELST jun-cc-pVQZ-RI
}

energy('sapt0')

compare_values(85.189064531275775, dimer.nuclear_repulsion_energy(), 9, "Nuclear Repulsion Energy")
compare_values(-0.00343130969, psi4.variable("SSAPT0 ELST ENERGY"), 6, "sSAPT0 elst")
compare_values( 0.00368418323, psi4.variable("SSAPT0 EXCH ENERGY"), 6, "sSAPT0 exch")
compare_values(-0.00093297498, psi4.variable("SSAPT0 IND ENERGY"), 6, "sSAPT0 ind")
compare_values(-0.00231534918, psi4.variable("SSAPT0 DISP ENERGY"), 6, "sSAPT0 disp")
compare_values(-0.00299545062, psi4.variable("SSAPT0 TOTAL ENERGY"), 6, "sSAPT0")
"""

    tfn = '_thread_test_input_psi4_yo'
    with open(tfn + '.in', 'w') as fp:
        fp.write(inputdat)

    run_psithon_inputs(args, tfn=tfn, label='Psi4')


def run_psithon_inputs(args, tfn, label):

    times = {}
    threads = [1, int(args.nthread)]

    for nt in threads:
        t0 = time.time()
        cmd = f"""psi4 -i {tfn}.in -o {tfn}_n{nt}.out -n{nt}"""
        print(f'Running {cmd} ...')
        subprocess.call(cmd.split())

        t1 = time.time()
        times[nt] = t1 - t0
        print("Time for threads %2d: Psi4: %12.6f" % (nt, times[nt]))

    rat1 = times[threads[0]] / times[threads[-1]]
    print("   Psi4@n%d : Psi4@n%d ratio (want ~%d): %.2f" % (threads[0], threads[-1], threads[-1], rat1))
    if args.passfail:
        #assert math.isclose(rat1, threads[-1], rel_tol=0.6), 'Psithon speedup {} !~= {}'.format(rat1, threads[-1])
        assert rat1 > 1.25, f'{label} Psithon speedup {rat1} !~= {threads[-1]}'


def test_plugin_dfmp2(args):

    inputdat = """
import %s

memory 2 gb

molecule {
0 1
C     0.0000000    0.0000000    1.0590353
C     0.0000000   -1.2060084    1.7576742
C     0.0000000   -1.2071767    3.1515905
C     0.0000000    0.0000000    3.8485751
C     0.0000000    1.2071767    3.1515905
C     0.0000000    1.2060084    1.7576742
H     0.0000000    0.0000000   -0.0215805
H     0.0000000   -2.1416387    1.2144217
H     0.0000000   -2.1435657    3.6929953
H     0.0000000    0.0000000    4.9301499
H     0.0000000    2.1435657    3.6929953
H     0.0000000    2.1416387    1.2144217
--
0 1
C    -1.3940633    0.0000000   -2.4541524
C    -0.6970468    1.2072378   -2.4546277
C     0.6970468    1.2072378   -2.4546277
C     1.3940633    0.0000000   -2.4541524
C     0.6970468   -1.2072378   -2.4546277
C    -0.6970468   -1.2072378   -2.4546277
H    -2.4753995    0.0000000   -2.4503221
H    -1.2382321    2.1435655   -2.4536764
H     1.2382321    2.1435655   -2.4536764
H     2.4753995    0.0000000   -2.4503221
H     1.2382321   -2.1435655   -2.4536764
H    -1.2382321   -2.1435655   -2.4536764
}

set {
  scf_type df
  mp2_type df
  basis aug-cc-pvdz
  freeze_core true
}

e, wfn = energy('plugdfmp2', return_wfn=True)
compare_values(-1.6309450762271729, wfn.variable('MP2 CORRELATION ENERGY'), 5, 'df-mp2 energy')  # aug-cc-pvdz
#compare_values(-1.5720781831194317, wfn.variable('MP2 CORRELATION ENERGY'), 5, 'df-mp2 energy')  # cc-pvdz
""" % (args.module)

    tfn = '_dfmp2_plugin_thread_test_input_psi4_yo'
    with open(tfn + '.in', 'w') as fp:
        fp.write(inputdat)

    run_psithon_inputs(args, tfn=tfn, label='Plugin dfmp2')


def print_math_ldd(args):

    module, sharedlibrary_woext = args.module.split('/')
    mod = importlib.import_module(module)
    exts = [sysconfig.get_config_var("EXT_SUFFIX"), '.so']
    for ext in exts:
        modcore = os.path.dirname(os.path.abspath(mod.__file__)) + os.path.sep + sharedlibrary_woext + ext
        if os.path.isfile(modcore):
            break

    if sys.platform.startswith('linux'):
        lddish = 'ldd -v'
    elif sys.platform.startswith('darwin'):
        lddish = 'otool -L'
    else:
        print('Not available w/o `ldd` or `otool`')
        return True

    cmd = """{} {} | grep -e ':' -e 'mkl' -e 'openblas' -e 'iomp5' -e 'gomp' -e 'libomp'""".format(lddish, modcore)
    print(f'Running {cmd} ...')
    subprocess.call(cmd, shell=True)
    lddout = subprocess.getoutput(cmd)
    report = {'mkl': lddout.count('libmkl'),
              'iomp5': lddout.count('libiomp5'),
              'openblas': lddout.count('libopenblas'),
              'omp': lddout.count('libomp'),
              'gomp': lddout.count('libgomp')}
    print(report)
    
    if sys.platform.startswith('linux'):
        slddout = collections.defaultdict(list)
        key = ''
        for ln in lddout.splitlines():
            if ':' in ln:
                key = ln.strip()
            else:
                slddout[key].append(ln.strip())

        for k, v in slddout.items():
            if modcore in k:
                tlddout = '\n'.join(v)
        treport = {'mkl': tlddout.count('libmkl'),
                   'iomp5': tlddout.count('libiomp5'),
                   'openblas': tlddout.count('libopenblas'),
                   'omp': tlddout.count('libomp'),
                   'gomp': tlddout.count('libgomp')}
        print(treport)
        if args.passfail:
            if sys.platform.startswith('linux'):
                assert (not treport['iomp5'] and not treport['omp'] and not treport['gomp']) is False

    report = {k : bool(v) for k, v in report.items()}
    okmkl = report['mkl'] and report['iomp5'] and not report['openblas'] and not report['gomp']
    okiomp5 = not report['mkl'] and report['iomp5'] and not report['openblas'] and not report['gomp']
    okopenblas = not report['mkl'] and not report['iomp5'] and report['openblas'] and report['gomp']
    
    omplike = (report['iomp5'] or report['omp']) and not report['gomp']
    okmkl2 = report['mkl'] and omplike and not report['openblas']
    if args.passfail:
        if sys.platform.startswith('linux'):
            assert okmkl2 != okopenblas
        elif sys.platform.startswith('darwin'):
            # plugins on Mac won't show mkl through otool (linked to psi4.core)
            assert (okmkl2 != okopenblas) or (omplike != okopenblas)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Psi4 threading tester. `psi4` and `import psi4` expected in *PATHs")
    parser.add_argument("-n", "--nthread", default=4,
                        help="""Number of threads to use. Psi4 disregards OMP_NUM_THREADS/MKL_NUM_THREADS.""")
    parser.add_argument("--passfail", action='store_true',
                        help="""Instead of just printing, run as tests.""")
    parser.add_argument("--module", default='psi4/core',
                        help="""In --ldd mode, module and shared library (w/o extension) to analyze, e.g., 'greatplugin/cxxcode.so' or 'psi4/core.cpython-36m-x86_64-linux-gnu.so'.
In --plugin-dfmp2 mode, name of dfmp2 module to load, e.g., 'plugdfmp2'.""")

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--psiapi', action='store_true',
                        help="""Test Psi4 (PsiAPI) vs NumPy in threaded matrix multiply.""")
    group.add_argument('--psithon', action='store_true',
                        help="""Test Psi4 (PSIthon) in threaded SAPT calc.""")
    group.add_argument('--ldd', action='store_true',
                        help="""Run ldd to examine BLAS and OMP linking of psi4/core.so (or whatever supplied to --module).""")
    group.add_argument('--plugin-dfmp2', action='store_true',
                        help="""Test dfmp2 plugin template (PSIthon) threading.""")

    args, unknown = parser.parse_known_args()

    if args.psiapi:
        test_threaded_blas(args)
    elif args.psithon:
        test_psithon(args)
    elif args.ldd:
        print_math_ldd(args)
    elif args.plugin_dfmp2:
        test_plugin_dfmp2(args)
    else:
        test_threaded_blas(args)
        test_psithon(args)
        print_math_ldd(args)


"""
PLUG="plugdfmp2"
THD=8
# * build psi4 and test its threading
PYTHONPATH=stage/lib/ python stage/share/psi4/scripts/test_threading.py --passfail --ldd
PATH=stage/bin/:$PATH PYTHONPATH=stage/lib/ python stage/share/psi4/scripts/test_threading.py --passfail -n$THD
# * build an OpenMP plugin and test its threading
stage/bin/psi4 --plugin-name $PLUG --plugin-template dfmp2
cd $PLUG && `../stage/bin/psi4 --plugin-compile` && make && cd ..
PYTHONPATH=stage/lib/:. python stage/share/psi4/scripts/test_threading.py --passfail --ldd --module="$PLUG/$PLUG"
PATH=stage/bin/:$PATH PYTHONPATH=stage/lib/:. python stage/share/psi4/scripts/test_threading.py --passfail --plugin-dfmp2 --module="$PLUG" -n$THD
"""
