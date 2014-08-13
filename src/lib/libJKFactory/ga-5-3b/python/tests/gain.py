import mpi4py.MPI as MPI
from ga4py import ga
from ga4py import gain
import numpy as np
from getopt import getopt
import sys
import traceback

def me():
    return ga.pgroup_nodeid(ga.pgroup_get_default())
def nproc():
    return ga.pgroup_nnodes(ga.pgroup_get_default())
def sync():
    ga.pgroup_sync(ga.pgroup_get_default())

PREPRINT = False
PRINTERS = [0]
#PRINTERS = [i for i in range(nproc())]

try:
    import colorama
    colorama.init()
    RED = colorama.Fore.RED
    YELLOW = colorama.Fore.YELLOW
    GREEN = colorama.Fore.GREEN
    RESET = colorama.Fore.RESET
except Exception:
    RED = ""
    YELLOW = ""
    GREEN = ""
    RESET = ""

# each test is exec'd by twice, once for numpy and once for gain
# 'm' is the current module
tests = [
    "result = m.arange(100, dtype=m.float32)",
    "result = m.arange(100, dtype=m.float32).copy()",
    "result = m.arange(100, dtype=m.float32)[9:17].copy()",
    """a = m.arange(100, dtype=m.float32)
       result = m.sin(a[49::-2],a[50::2])""",
    """a = m.arange(100, dtype=m.float32)
       result = m.sin(a[49::-2],a[50::2])""",
    """a = m.arange(100, dtype=m.float32)
       a[5:15] = m.arange(10, dtype=m.float32)
       result = a""",
    """a = m.arange(100, dtype=m.float32)
       a[5:15] = np.arange(10, dtype=m.float32)
       result = a""",
    "result = m.sin(1)",
    "result = m.sin([1,2,3])",
    "result = m.sin([1,2,3], m.asarray([0,0,0], dtype=m.float32))""",
    """c = m.ones(10, dtype=m.int16)
       d = m.ones(10, dtype=m.int16)
       result = m.sin(c,d)""",
    "result = m.sin(m.ones(10, dtype=m.int16))",
    """a = m.arange(100, dtype=m.float32)
       b = a[:50:2]
       c = a[50::2]
       result = m.add(b,c)""",
    # commented out because numpy produces inconsistent results
    # whereas gain produces the 'correct' results due to internal copies
    #"""a = m.arange(100, dtype=m.float32)
    #   b = a[:50:2]
    #   c = a[50::2]
    #   d = a[::4]
    #   result = m.add(b,c,d)""",
    """a = m.arange(100, dtype=m.float32)
       b = a[:50:2]
       result = m.add(b,5)""",
    """a = m.arange(100, dtype=m.float32)
       b = a[50::2]
       result = m.add(5,b)""",
    "result = m.ones((3,4,5), dtype=m.float32)",
    """s = m.ones((3,4,5), dtype=m.float32)
       t = m.ones((3,4,5), dtype=m.float32)
       result = m.sin(s,t)""",
    """x = m.ones((2,3,4))
       result = m.add(x,x)""",
    """x = m.ones((2,3,4))
       y = m.ones((3,4))
       result = m.add(x,y)""",
    "result = m.linspace(2.0,3.0,num=5)",
    "result = m.linspace(2.0,3.0,num=5,endpoint=False)",
    "result = m.linspace(2.0,3.0,num=5,retstep=True)",
    "result = m.logspace(2.0,3.0,num=4)",
    "result = m.logspace(2.0,3.0,num=4,endpoint=False)",
    "result = m.logspace(2.0,3.0,num=4,base=2.0)",
    """f = m.ones((100,200), dtype=float)
       g = m.ones((200,300), dtype=float)
       result = m.dot(f,g)""",
    """h = m.arange(100, dtype=m.float32)
       result = m.dot(h,h)""",
    "result = m.dot(6,10)",
    """h = m.arange(100, dtype=m.float32)
       result = m.dot(6,h)""",
    "result = m.eye(24,25)",
    "result = m.eye(24,25,4)",
    "result = m.eye(24,25,-8)",
    "result = m.identity(11)",
    "result = m.add.reduce([1,2,3,4])",
    "result = m.add.reduce(m.arange(100))",
    "result = m.add.reduce(m.ones((100,200)))",
    "result = m.add.reduce(m.ones((100,200,300)))",
    "result = m.subtract.reduce([1,2,3,4])",
    "result = m.subtract.reduce(m.arange(100))",
    """i = m.arange(100*200).reshape((100,200))
       result = m.subtract.reduce(i)""",
    """i = m.arange(100*200*300).reshape((100,200,300))
       result = m.subtract.reduce(i)""",
    "result = m.add.accumulate(m.ones(7))",
    "result = m.add.accumulate(m.ones((7,7)))",
    "result = m.add.accumulate(m.ones((7,7,7)), axis=0)",
    "result = m.add.accumulate(m.ones((7,7,7)), axis=1)",
    "result = m.add.accumulate(m.ones((7,7,7)), axis=2)",
    "result = m.alen((1,2,3))",
    "result = m.alen(m.zeros((4,5,6)))",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       result = i""",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       result = i.flat""",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       result = i.flat[2]""",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       result = i.flat[2:19]""",
    """i = m.zeros((4,25))
       j = i.flat
       j[:] = 6
       result = j""",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       j = i.flat
       j[:] = 6
       j[3:19] = 7
       result = j""",
    """foo = m.zeros((3,4,5))
       bar = m.arange(3*4*5)
       foo.flat = bar
       result = foo""",
    """foo = m.zeros((3,4,5))
       foo.flat = 6
       result = foo""",
    """foo = m.zeros((3,4,5))
       bar = m.arange(3*4*5)
       i = m.zeros((3,4,5))
       result = m.add(foo.flat,bar,i.flat)""",
    """foo = m.zeros((3,4,5))
       bar = m.arange(3*4*5)
       foo.flat = bar
       baz = foo[1:,::2,2:5]
       result = baz.flat[2:7]""",
    """foo = m.zeros((3,4,5))
       bar = m.arange(3*4*5)
       foo.flat = bar
       baz = foo[1,1:,2:]
       result = baz.flat[2:7]""",
    """foo = m.zeros((3,4,5))
       bar = m.arange(3*4*5)
       foo.flat = bar
       baz = foo[1,1:,None,2:]
       result = baz.flat[2:7]""",
    "result = m.clip(m.arange(10), 1, 8)",
    "result = m.clip(m.arange(100), 10, 80)",
    "result = m.clip(m.arange(10), [3,4,1,1,1,4,4,4,4,4], 8)",
    """foo = np.arange(4*25).reshape(4,25)
       i = m.zeros((4,25))
       i[:] = foo
       result = i.transpose()""",
    """foo = np.arange(4*5*77).reshape(4,5,77)
       k = m.zeros((4,5,77))
       k[:] = foo
       result = k.transpose()""",
    """foo = np.arange(4*5*77).reshape(4,5,77)
       k = m.zeros((4,5,77))
       k[:] = foo
       result = k.transpose(1,2,0).shape""",
    """foo = np.arange(4*5*77).reshape(4,5,77)
       k = m.zeros((4,5,77))
       k[:] = foo
       result = k.transpose(1,2,0)""",
    "result = m.arange(3*4*5).reshape((3,4,5))",
    "result = m.arange(3*4*5).reshape((3,4,5))[1:,2:,::2].reshape((2,3,2))",
    """a = m.zeros((10,2))
       a.shape = (20)
       result = a""",
    """a = m.empty((4,2,6,92), dtype=m.float32)
       a.fill(6789.0)
       result = a""",
    """a = m.empty((4,2,6,92), dtype=m.int32)
       a.fill(6789.0)
       result = a""",
    "result = m.indices((3,4,5))",
    "result = m.indices((3,4,5)).max()",
    "result = m.indices((3,4,5)).max(0)",
    "result = m.arange(3*4*5).reshape((5,4,3))",
    "result = m.arange(3*4*5).reshape((5,4,3)).T",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[1:,2:,3:]",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[0,2:,3:]",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[1:,0,3:]",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[1:,2:,0]",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[1:,None,2:,0]",
    "result = m.arange(3*4*5).reshape((5,4,3)).T[None,1:,None,2:,None,0,None]",
    "result = m.add.reduce(m.arange(3*4*5).reshape((5,4,3)).T)",
    """result = m.ones((15,2))
       result = m.multiply(result,5,result)""",
    """a = m.ones((15,2))
       b = a[:,0]
       result = m.multiply(b,5,b)""",
    """result = m.ones((15,2)) * 5""",
    """result = m.ones((15,2))
       result *= 5""",
    """result = m.ones((15,2))[:,0] * 5""",
    """result = m.ones((15,2))
       result[:,0] *= 5""",
    "result = m.arange(3*4).reshape(3,4).diagonal()",
    "result = m.arange(3*4).reshape(3,4).diagonal(1)",
    "result = m.arange(3*4).reshape(3,4).diagonal(-1)",
]

# the current module, either numpy or gain
m = None

# count test results
passes = 0
x_failures = 0
np_failures = 0
gain_failures = 0
epic_failures = 0

class PrintZero(object):
    def __init__(self):
        self.stdout = sys.stdout
    def write(self, something):
        if me() in PRINTERS:
            self.stdout.write(something)
    def flush(self):
        if me() in PRINTERS:
            self.stdout.flush()
sys.stdout = PrintZero()

def _dtype(thing):
    try:
        return thing.dtype
    except:
        try:
            return thing.base.dtype
        except:
            return None
def _shape(thing):
    try:
        return thing.shape
    except:
        try:
            return thing.base.size # it's a flatiter
        except:
            return None

def print_result(result_np,result_gain,diff):
    if isinstance(result_np, np.flatiter):
        _result_np = result_np.copy()
    else:
        _result_np = result_np
    print """%s
---------------------- numpy ---------------------------------
%s
%s
%s
%s
---------------------- gain ----------------------------------
%s
%s
%s
%s
---------------------- difference ----------------------------
%s
%s
""" % (RED,
    type(result_np),   _dtype(result_np),   _shape(result_np),   _result_np,
    type(result_gain), _dtype(result_gain), _shape(result_gain),  result_gain,
    diff, RESET)

def run_tests():
    global m
    for test in tests:
        run_test(test)

def run_test(test):
    global passes,x_failures,np_failures,gain_failures,epic_failures
    # sanity check that the test is correctly written
    if test.count("result") < 1:
        raise SyntaxError, "TEST ERROR: %s" % test
    # clean up whitespace and pretty print the test string
    test_lines = [line.strip() for line in test.splitlines()]
    test = '\n'.join(test_lines)
    if PREPRINT:
        print " TESTING:"
        print '    ' + '\n    '.join(test_lines)
    any_err = False
    e_np = None
    e_gain = None
    result = None
    result_np = None
    result_gain = None
    m = np
    # some temporary labels for results
    try:
        exec test
        result_np = result
    except Exception,e:
        e_np = e
        tb_np = sys.exc_info()
    m = gain
    try:
        exec test
        result_gain = result
    except Exception,e:
        e_gain = e
        tb_gain = sys.exc_info()
    # sync now since most operations sync on the way in, not on the way out
    sync()
    if e_np is None and e_gain is None:
        err = False
        hard_err = False
        hard_err_tb = None
        try:
            diff = None
            if isinstance(result_gain, (gain.ndarray,gain.flatiter)):
                diff = result_np-result_gain.get()
            elif isinstance(result_gain, tuple):
                if (len(result_gain) > 1
                    and type(result_gain[0]) is type(result_gain[1])):
                    diff = np.asarray(result_np)-np.asarray(result_gain)
                else:
                    result_np = result_np[0]
                    result_gain = result_gain[0]
                    diff = result_np-result_gain.get()
            else:
                diff = result_np-result_gain
            if not np.all(diff == 0):
                print_result(result_np,result_gain,diff)
                err = True
            if not (_dtype(result_np) == _dtype(result_gain)):
                print RED + "different types np=%s gain=%s" % (
                        _dtype(result_np), _dtype(result_gain)) + RESET
                err = True
        except Exception,e:
            print "%scaught exception: %s%s" % (RED,e,RESET)
            hard_err = True
            hard_err_tb = sys.exc_info()
        if hard_err:
            print " RESULT: %sEPIC FAIL%s" % (RED,RESET)
            print "".join(traceback.format_exception(*hard_err_tb))
            epic_failures += 1
            any_err = True
        elif err:
            print " RESULT: %sFAIL%s" % (RED,RESET)
            gain_failures += 1
            any_err = True
        else:
            #print " RESULT: %sPASS%s" % (GREEN,RESET)
            passes += 1
    elif e_np is None:
        print " RESULT: %sFAIL -- gain exception only: %s%s" % (
                RED,e_gain,RESET)
        print "".join(traceback.format_exception(*tb_gain))
        gain_failures += 1
        any_err = True
    elif e_gain is None:
        print " RESULT: %sFAIL -- numpy exception only: %s%s" % (
                RED,e_np,RESET)
        print "".join(traceback.format_exception(*tb_np))
        np_failures += 1
        any_err = True
    else: # both errors are set
        if str(e_np) != str(e_gain):
            print " RESULT: %sFAIL (dffering exceptions)%s" % (RED,RESET)
            print "   %snp:%s'%s'%s" % (RED,type(e_np),e_np,RESET)
            print " %sgain:%s'%s'%s" % (RED,type(e_gain),e_gain,RESET)
            print " --- NumPy traceback ---"
            print "".join(traceback.format_exception(*tb_np))
            print " --- GAiN traceback ---"
            print "".join(traceback.format_exception(*tb_gain))
            gain_failures += 1
        else:
            #print " RESULT: %sXFAIL%s" % (GREEN,RESET)
            xfailures += 1
        any_err = True
    if any_err:
        print "TEST WAS: %s" % test_lines[0]
        for line in test_lines[1:]:
            print "         %s" % line

if __name__ == '__main__':
    profile = False
    use_groups = False
    use_color = True
    (optsvals,args) = getopt(sys.argv[1:],'pgc')
    for (opt,val) in optsvals:
        if opt == '-p':
            profile = True
        elif opt == '-g':
            use_groups= True
        elif opt == '-c':
            use_color = False
    if not use_color:
        RED = ""
        YELLOW = ""
        GREEN = ""
        RESET = ""
    if profile:
        import cProfile
        print "Profiling enabled"
        cProfile.run("run_tests()", "gaintest.%s.prof" % str(me()))
    elif use_groups:
        midproc = nproc()//2
        proclist_first = range(0,midproc)
        proclist_last  = range(midproc,nproc())
        group_id_first = ga.pgroup_create(proclist_first)
        group_id_last  = ga.pgroup_create(proclist_last)
        if me() in proclist_first:
            ga.pgroup_set_default(group_id_first)
            run_tests()
        ga.pgroup_set_default(ga.pgroup_get_world())
        sync()
        if me() in proclist_last:
            ga.pgroup_set_default(group_id_last)
            run_tests()
        ga.pgroup_set_default(ga.pgroup_get_world())
        sync()
        print "All done with groups"
    else:
        run_tests()
    print ""
    print "%s           Passed: %s%s" % (GREEN,passes,RESET)
    print "%sExpected Failures: %s%s" % (GREEN,x_failures,RESET)
    print "%s   NumPy Failures: %s%s" % (YELLOW,np_failures,RESET)
    print "%s    GAiN Failures: %s%s" % (RED,gain_failures,RESET)
    print "%s    Epic Failures: %s%s" % (RED,epic_failures,RESET)
    if not me():
        ga.print_stats()
