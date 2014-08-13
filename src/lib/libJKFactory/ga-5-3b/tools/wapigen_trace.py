#!/usr/bin/env python

'''Generate the wapi_trace.c source from the ga-papi.h header.'''

import sys

def get_signatures(header):
    # first, gather all function signatures from ga-papi.h aka argv[1]
    accumulating = False
    signatures = []
    current_signature = ''
    EXTERN = 'extern'
    SEMICOLON = ';'
    for line in open(header):
        line = line.strip() # remove whitespace before and after line
        if not line:
            continue # skip blank lines
        if EXTERN in line and SEMICOLON in line:
            signatures.append(line)
        elif EXTERN in line:
            current_signature = line
            accumulating = True
        elif SEMICOLON in line and accumulating:
            current_signature += line
            signatures.append(current_signature)
            accumulating = False
        elif accumulating:
            current_signature += line
    return signatures

format_from_type = {
        "char":         "%c", 
        "int":          "%d",
        "long":         "%ld",
        "long long":    "%lld",
        "float":        "%f",
        "double":       "%lf",
        "Integer":      "%ld",
        "logical":      "%ld",
        }

intp_names = [
        'lo',
        'alo',
        'blo',
        'clo',
        'hi',
        'ahi',
        'bhi',
        'chi',
        'subscript',
        'dims',
        'ld',
        'width',
        'chunk',
        'map',
        'block',
        ]

class FunctionArgument(object):

    def __init__(self, signature):
        self.pointer = signature.count('*')
        self.array = '[' in signature
        signature = signature.replace('*','').strip()
        signature = signature.replace('[','').strip()
        signature = signature.replace(']','').strip()
        self.type,self.name = signature.split()
        self.intp = ((self.pointer == 1 or self.array)
                and self.type == 'Integer'
                and self.name in intp_names)

    def __str__(self):
        ret = self.type[:]
        ret += ' '
        for p in range(self.pointer):
            ret += '*'
        ret += self.name
        if self.array:
            ret += '[]'
        return ret

class Function(object):

    def __init__(self, signature):
        signature = signature.replace('extern','').strip()
        self.return_type,signature = signature.split(None,1)
        self.return_type = self.return_type.strip()
        signature = signature.strip()
        self.name,signature = signature.split('(',1)
        self.name = self.name.strip()
        signature = signature.replace(')','').strip()
        signature = signature.replace(';','').strip()
        self.args = []
        if signature:
            for arg in signature.split(','):
                self.args.append(FunctionArgument(arg.strip()))

    def get_call(self, name=None):
        sig = ''
        if not name:
            sig += self.name
        else:
            sig += name
        sig += '('
        if self.args:
            for arg in self.args:
                sig += arg.name
                sig += ', '
            sig = sig[:-2] # remove last ', '
        sig += ')'
        return sig

    def get_signature(self, name=None):
        sig = self.return_type[:]
        sig += ' '
        if not name:
            sig += self.name
        else:
            sig += name
        sig += '('
        if self.args:
            for arg in self.args:
                sig += str(arg)
                sig += ', '
            sig = sig[:-2] # remove last ', '
        sig += ')'
        return sig

    def args_have_nv(self):
        for arg in self.args:
            if arg.name == 'nv' and not arg.pointer:
                return True
        return False

    def args_have_ga(self):
        for arg in self.args:
            if arg.name == 'g_a' and not arg.pointer:
                return True
        return False

    def args_have_alpha(self):
        for arg in self.args:
            if arg.name == 'alpha':
                return True
        return False

    def args_have_beta(self):
        for arg in self.args:
            if arg.name == 'beta':
                return True
        return False

    def args_have_ndim(self):
        for arg in self.args:
            if 'ndim' in arg.name and not arg.pointer:
                return True
        return False

    def get_ndim_name(self):
        if 'set_irreg' in self.name:
            return 'pnga_set_data_last_ndim'
        else:
            for arg in self.args:
                if 'ndim' in arg.name and not arg.pointer:
                    return arg.name
            return '_ndim'

    def args_have_intp(self):
        for arg in self.args:
            if arg.intp:
                return True
        return False

    def get_tracer_body(self):
        slot = 0
        tracer = ''
        if self.args_have_intp() and 'set_irreg' not in self.name:
            if not self.args_have_ndim() and self.args_have_ga():
                tracer += '    int _ndim = pnga_ndim(g_a);\n'
        if self.args_have_alpha():
            tracer += '    Integer _atype;\n'
        if self.args_have_beta():
            tracer += '    Integer _btype;\n'
        if self.args_have_alpha():
            tracer += '    pnga_inquire_type(g_a, &_atype);\n'
        if self.args_have_beta():
            tracer += '    pnga_inquire_type(g_b, &_btype);\n'
        if 'void' not in func.return_type:
            tracer += '    %s retval;\n' % self.return_type
        pf,slot = self.get_tracer_printf(False,slot)
        tracer += '    %s;\n' % pf
        tracer += '    '
        if 'void' not in func.return_type:
            tracer += 'retval = '
        tracer += '%s;\n' % self.get_call()
        pf,slot = self.get_tracer_printf(True,slot)
        tracer += '    %s;\n' % pf
        if 'void' not in func.return_type:
            tracer += '    return retval;\n'
        return tracer

    def get_tracer_printf(self, end, slot):
        tracer = 'fprintf(fptrace, "%lf,'
        if end: tracer += '/'
        tracer += self.name + ','
        if self.args:
            tracer += '('
            for arg in self.args:
                if arg.pointer == 1 and 'char' in arg.type:
                    tracer += '%s;'
                elif arg.name == 'alpha' or arg.name == 'beta':
                    tracer += '%s;'
                elif arg.intp:
                    tracer += '%s;'
                elif arg.pointer or arg.array:
                    tracer += '%p;'
                else:
                    tracer += '%s;' % format_from_type[arg.type]
            tracer = tracer[:-1]
            tracer += ')'
        tracer += '\\n",MPI_Wtime()-first_wtime'
        if self.args:
            for arg in self.args:
                if arg.intp:
                    #if self.args_have_nv():
                    #    tracer += ',expand_intp(%s*nv,%s,%d)' % (
                    #            self.get_ndim_name(), arg.name, slot)
                    if 'map' in arg.name and 'locate_region' not in self.name:
                        tracer += ',expand_intp(sum_intp(%s,block),%s,%d)' % (
                                self.get_ndim_name(), arg.name, slot)
                    elif 'ld' in arg.name:
                        tracer += ',expand_intp(%s-1,%s,%d)' % (
                                self.get_ndim_name(), arg.name, slot)
                    else:
                        tracer += ',expand_intp(%s,%s,%d)' % (
                                self.get_ndim_name(), arg.name, slot)
                    slot += 1
                elif 'alpha' in arg.name:
                    tracer += ',expand_voidp(_atype,alpha,%d)' % slot
                    slot += 1
                elif 'beta' in arg.name:
                    tracer += ',expand_voidp(_btype,beta,%d)' % slot
                    slot += 1
                else:
                    tracer += ',%s' % arg.name
        tracer += ')'
        return tracer,slot

    def __str__(self):
        return self.get_signature()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'incorrect number of arguments'
        print 'usage: wapigen_trace.py <ga-papi.h> > <wapi_trace.c>'
        sys.exit(len(sys.argv))

    # print headers and other static stuff
    print '''
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ga-papi.h"
#include "typesf2c.h"

FILE *fptrace=NULL;
double first_wtime;
int me, nproc;

#if HAVE_PROGNAME
extern const char * PROGNAME;
#endif

#define LOHI_BUFSIZE 200
#define LOHI_SLOTS 20

static char **lohi_bufs;
static int pnga_set_data_last_ndim=-1;

static Integer sum_intp(int ndim, Integer *lohi) {
    Integer sum=0;

    assert(ndim>=0);
    if (NULL != lohi) {
        int i;
        for (i=0; i<ndim; i++) {
            sum += lohi[i];
        }
    }

    return sum;
}

static char* expand_voidp(int type, void *value, int slot) {
    char *str=NULL, *current_str=NULL;
    int size_written;
    int total_written = 0;

    assert(slot >= 0 && slot < LOHI_SLOTS);
    str = lohi_bufs[slot];
    current_str = str;
    if (NULL == value) {
        size_written = sprintf(current_str, "NULL");
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
    }
    else {
        switch (type) {
            case C_INT:     size_written = sprintf(current_str, "%d",
                            *((int*)value));
                            break;
            case C_LONG:    size_written = sprintf(current_str, "%ld",
                            *((long*)value));
                            break;
            case C_LONGLONG:size_written = sprintf(current_str, "%lld",
                            *((long long*)value));
                            break;
            case C_FLOAT:   size_written = sprintf(current_str, "%f",
                            *((float*)value));
                            break;
            case C_DBL:     size_written = sprintf(current_str, "%lf",
                            *((double*)value));
                            break;
            case C_SCPL:    size_written = sprintf(current_str, "%f#%f",
                            ((SingleComplex*)value)->real,
                            ((SingleComplex*)value)->imag);
                            break;
            case C_DCPL:    size_written = sprintf(current_str, "%lf#%lf",
                            ((DoubleComplex*)value)->real,
                            ((DoubleComplex*)value)->imag);
                            break;
        }
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
    }

    return str;
}

static char* expand_intp(int ndim, Integer *lohi, int slot) {
    char *str=NULL, *current_str=NULL;
    int i;
    int size_written;
    int total_written = 0;

    assert(ndim>=0);
    assert(slot >= 0 && slot < LOHI_SLOTS);
    str = lohi_bufs[slot];
    current_str = str;
    if (NULL == lohi) {
        size_written = sprintf(current_str, "{NULL}");
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
    }
    else if (0 == ndim) {
        size_written = sprintf(current_str, "{}");
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
    }
    else {
        size_written = sprintf(current_str, "{%ld", lohi[0]);
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
        for (i=1; i<ndim; i++) {
            size_written = sprintf(current_str, ":%ld", lohi[i]);
            total_written += size_written;
            assert(size_written > 0);
            assert(total_written < LOHI_BUFSIZE);
            current_str += size_written;
        }
        size_written = sprintf(current_str, "}\\0");
        total_written += size_written;
        assert(size_written > 0);
        assert(total_written < LOHI_BUFSIZE);
        current_str += size_written;
    }

    return str;
}  

static void init_lohi_bufs() {
    int i;
    lohi_bufs = (char**)malloc(sizeof(char*)*LOHI_SLOTS);
    for (i=0; i<LOHI_SLOTS; i++) {
        lohi_bufs[i] = (char*)malloc(LOHI_BUFSIZE);
    }
}

static void free_lohi_bufs() {
    int i;
    for (i=0; i<LOHI_SLOTS; i++) {
        free(lohi_bufs[i]);
    }
    free(lohi_bufs);
}

static void reset_lohi_bufs() {
    int i;
    for (i=0; i<LOHI_SLOTS; i++) {
        memset(&(lohi_bufs[i][0]), 0, LOHI_BUFSIZE);
    }
}

static void trace_finalize() {
    fclose(fptrace);
    free_lohi_bufs();
}

static void trace_initialize() {
    /* create files to write trace data */
    char *profile_dir=NULL;
    const char *program_name=NULL;
    char *file_name=NULL;
    struct stat f_stat;

    PMPI_Barrier(MPI_COMM_WORLD);
    PMPI_Comm_rank(MPI_COMM_WORLD, &me);
    PMPI_Comm_size(MPI_COMM_WORLD, &nproc);

    first_wtime = MPI_Wtime();
    init_lohi_bufs();

    profile_dir = getenv("PNGA_PROFILE_DIR");
#if HAVE_PROGNAME
    program_name = PROGNAME;
#else
    program_name = "unknown";
#endif
    if (0 == me) {
        int ret;

        if (!profile_dir) {
            fprintf(stderr, "You need to set PNGA_PROFILE_DIR env var\\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fprintf(stderr, "PNGA_PROFILE_DIR=%s\\n", profile_dir);
        if (-1 == stat(profile_dir, &f_stat)) {
            perror("stat");
            fprintf(stderr, "Cannot successfully stat to PNGA_PROFILE_DIR.\\n");
            fprintf(stderr, "Check %s profile dir\\n", profile_dir);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file_name = (char *)malloc(strlen(profile_dir)
                + 1 /* / */
                + strlen(program_name)
                + 2 /* NULL termination */);
        assert(file_name);
        sprintf(file_name, "%s/%s%c\\n", profile_dir, program_name, '\\0');
        ret = mkdir(file_name, 0755);
        if (ret) {
            perror("mkdir");
            fprintf(stderr, "%d: profile sub-directory creation failed: pathname=%s: exiting\\n", me, file_name);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        free(file_name);
    }
    PMPI_Barrier(MPI_COMM_WORLD);
    file_name = (char *)malloc(strlen(profile_dir)
            + 1 /* / */
            + strlen(program_name)
            + 1 /* / */
            + 7 /* mpi id */
            + 6 /* .trace */
            + 2 /* NULL termination */);
    assert(file_name);
    sprintf(file_name,"%s/%s/%07d.trace%c",profile_dir,program_name,me,'\\0');
    fptrace = fopen(file_name,"w");
    if(!fptrace) {
        perror("fopen");
        printf("%d: trace file creation failed: file_name=%s: exiting\\n", me, file_name);
        exit(0);
    }
    free(file_name);
}

'''

    functions = {}
    # parse signatures into the Function class
    for sig in get_signatures(sys.argv[1]):
        function = Function(sig)
        functions[function.name] = function

    # now process the functions
    for name in sorted(functions):
        func = functions[name]
        if name in ['pnga_initialize','pnga_initialize_ltd','pnga_terminate',
                'pnga_set_data']:
            continue
        func = functions[name]
        wnga_name = name.replace('pnga_','wnga_')
        print '''
%s
{
%s
}
''' % (func.get_signature(wnga_name), func.get_tracer_body())

    # output the pnga_set_data function
    name = 'pnga_set_data'
    wnga_name = name.replace('pnga_','wnga_')
    func = functions[name]
    print '''
%s
{
    pnga_set_data_last_ndim = ndim;
%s
}
''' % (func.get_signature(wnga_name), func.get_tracer_body())

    # output the initialize function
    name = 'pnga_initialize'
    wnga_name = name.replace('pnga_','wnga_')
    func = functions[name]
    print '''%s
{
    static int count_pnga_initialize=0;

    ++count_pnga_initialize;
    if (1 == count_pnga_initialize) {
        trace_initialize();
    }
%s
}
''' % (func.get_signature(wnga_name), func.get_tracer_body())

    # output the initialize_ltd function
    name = 'pnga_initialize_ltd'
    wnga_name = name.replace('pnga_','wnga_')
    func = functions[name]
    print '''%s
{
    static int count_pnga_initialize_ltd=0;

    ++count_pnga_initialize_ltd;
    if (1 == count_pnga_initialize_ltd) {
        trace_initialize();
    }
%s
}
''' % (func.get_signature(wnga_name), func.get_tracer_body())

    # prepare to output the terminate function
    name = 'pnga_terminate'
    wnga_name = name.replace('pnga_','wnga_')
    func = functions[name]
    # output the terminate function
    print '''%s
{
    static int count_pnga_terminate=0;

    ++count_pnga_terminate;
%s
    /* don't dump info if terminate more than once */
    if (1 == count_pnga_terminate) {
        trace_finalize();
    }
}
''' % (func.get_signature(wnga_name), func.get_tracer_body()) 
