#!/usr/bin/env python

'''Generate the armci_prof.c source from the parmci.h header.'''

import sys

def get_signatures(header):
    # first, gather all function signatures from parmci.h aka argv[1]
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

class FunctionArgument(object):
    def __init__(self, signature):
        self.pointer = signature.count('*')
        self.array = '[' in signature
        signature = signature.replace('*','').strip()
        signature = signature.replace('[','').strip()
        signature = signature.replace(']','').strip()
        self.type,self.name = signature.split()

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
        if '*' not in self.return_type and signature[0] == '*':
            # return type is void* not void
            self.return_type += '*'
            signature = signature[1:].strip()
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

    def __str__(self):
        return self.get_signature()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'incorrect number of arguments'
        print 'usage: prof_gen.py <parmci.h> > <armci_prof.c>'
        sys.exit(len(sys.argv))

    # print headers
    print '''
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include <mpi.h>

#include "armci.h"
#include "parmci.h"

static MPI_Comm comm = MPI_COMM_NULL;
static int me = -1;
static int nproc = -1;
static double global_start = 0.0;
static double global_stop = 0.0;
'''

    functions = {}
    # parse signatures into the Function class
    for sig in get_signatures(sys.argv[1]):
        function = Function(sig)
        functions[function.name] = function

    # for each function, generate a static count
    for name in sorted(functions):
        print 'static long count_%s = 0;' % name
    print ''

    # for each function, generate a static time
    for name in sorted(functions):
        print 'static double time_%s = 0;' % name
    print ''

    # now process the functions
    for name in sorted(functions):
        func = functions[name]
        maybe_declare = ''
        maybe_assign = ''
        maybe_return = ''
        if '*' in func.return_type or 'void' not in func.return_type:
            maybe_declare = '%s ret;' % func.return_type
            maybe_assign = 'ret = '
            maybe_return = 'return ret;'
        new_name = None
        if 'PARMCI_' in name:
            new_name = name.replace('PARMCI_','ARMCI_')
        elif 'parmci_' in name:
            new_name = name.replace('parmci_','armci_')
        if new_name in ['ARMCI_Finalize']:
            continue
        elif new_name in ['ARMCI_Init','ARMCI_Init_args']:
            print '''
%s
{
    double local_start;
    double local_stop;
    %s
    if (comm == MPI_COMM_NULL) {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
        MPI_Comm_rank(comm, &me);
        MPI_Comm_size(comm, &nproc);
    }
    ++count_%s;
    global_start = MPI_Wtime();
    local_start = MPI_Wtime();
    %s%s;
    local_stop = MPI_Wtime();
    time_%s += local_stop - local_start;
    %s
}
''' % (func.get_signature(new_name),
        maybe_declare,
        name,
        maybe_assign, func.get_call(),
        name,
        maybe_return)
        else:
            print '''
%s
{
    double local_start;
    double local_stop;
    %s
    ++count_%s;
    local_start = MPI_Wtime();
    %s%s;
    local_stop = MPI_Wtime();
    time_%s += local_stop - local_start;
    %s
}
''' % (func.get_signature(new_name),
        maybe_declare,
        name,
        maybe_assign, func.get_call(),
        name,
        maybe_return)

    # prepare to output the terminate function
    name = 'PARMCI_Finalize'
    new_name = 'ARMCI_Finalize'
    func = functions[name]
    the_code = ''
    # establish 'the_code' to use in the body of terminate
    # it's printing the count of each function if it was called at least once
    the_code += '''
        double recvbuf = 0.0;
'''
    for fname in sorted(functions):
        the_code += '''
        MPI_Reduce(&time_%s, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("%s,%%ld,%%lf\\n", count_%s, recvbuf);
        }
''' % (fname, fname, fname)
    # output the terminate function
    print '''%s
{
    int ret;
    ++count_%s;
    %s;
    global_stop = MPI_Wtime();
    /* don't dump info if terminate more than once */
    if (1 == count_%s) {
%s
        if (me == 0) {
            printf("global_stop-global_start,0,%%lf\\n",
                    global_stop-global_start);
        }
    }
    MPI_Comm_free(&comm);
}
''' % (func.get_signature(new_name),
        name,
        func.get_call(),
        name,
        the_code)

