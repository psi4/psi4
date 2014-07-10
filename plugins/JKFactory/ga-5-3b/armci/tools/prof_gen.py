#!/usr/bin/env python

'''Generate the armci_prof.c source from the parmci.h header.'''

import sys

name_to_event = {
"ARMCI_Acc"               : "ARMCI_PROF_ACC",
"ARMCI_AccS"              : "ARMCI_PROF_ACCS",
"ARMCI_AccV"              : "ARMCI_PROF_ACCV",
"ARMCI_AllFence"          : "ARMCI_PROF_ALLFENCE",
"ARMCI_Barrier"           : "ARMCI_PROF_BARRIER",
"ARMCI_Fence"             : "ARMCI_PROF_FENCE",
"ARMCI_Get"               : "ARMCI_PROF_GET",
"ARMCI_GetS"              : "ARMCI_PROF_GETS",
"ARMCI_GetV"              : "ARMCI_PROF_GETV",
"ARMCI_NbAcc"             : "ARMCI_PROF_NBACC",
"ARMCI_NbAccS"            : "ARMCI_PROF_NBACCS",
"ARMCI_NbAccV"            : "ARMCI_PROF_NBACCV",
"ARMCI_NbGet"             : "ARMCI_PROF_NBGET",
"ARMCI_NbGetS"            : "ARMCI_PROF_NBGETS",
"ARMCI_NbGetV"            : "ARMCI_PROF_NBGETV",
"ARMCI_NbPut"             : "ARMCI_PROF_NBPUT",
"ARMCI_NbPutS"            : "ARMCI_PROF_NBPUTS",
"ARMCI_NbPutV"            : "ARMCI_PROF_NBPUTV",
"ARMCI_Put"               : "ARMCI_PROF_PUT",
"ARMCI_PutS"              : "ARMCI_PROF_PUTS",
"ARMCI_PutV"              : "ARMCI_PROF_PUTV",
"ARMCI_Rmw"               : "ARMCI_PROF_RMW",
"ARMCI_Wait"              : "ARMCI_PROF_WAIT",
"armci_msg_barrier"       : "ARMCI_PROF_BARRIER",
"armci_msg_group_barrier" : "ARMCI_PROF_BARRIER",
"armci_notify_wait"       : "ARMCI_PROF_NOTIFY",
}

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

#include "armci.h"
#include "parmci.h"
#include "armci_profile.h"
#include "armci_profile.c"
'''

    functions = {}
    # parse signatures into the Function class
    for sig in get_signatures(sys.argv[1]):
        function = Function(sig)
        functions[function.name] = function

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
        func = functions[name]
        new_name = None
        if 'PARMCI_' in name:
            new_name = name.replace('PARMCI_','ARMCI_')
        elif 'parmci_' in name:
            new_name = name.replace('parmci_','armci_')
        if new_name in ['ARMCI_Init','ARMCI_Init_args']:
            print '''
%s
{
    int ret;
    ret = %s;
    armci_profile_init();
    return ret;
}
''' % (func.get_signature(new_name), func.get_call())
        elif new_name in ['ARMCI_Finalize']:
            print '''
%s
{
    armci_profile_terminate();
    %s;
}
''' % (func.get_signature(new_name), func.get_call())
        elif new_name in ['ARMCI_GetV',   'ARMCI_PutV',   'ARMCI_AccV',
                          'ARMCI_NbGetV', 'ARMCI_NbPutV', 'ARMCI_NbAccV']:
            print '''
%s
{
    %s
    armci_profile_start_vector(darr, len, proc, %s);
    %s%s;
    armci_profile_stop_vector(%s);
    %s
}
''' % (func.get_signature(new_name),
                maybe_declare,
                name_to_event[new_name],
                maybe_assign, func.get_call(),
                name_to_event[new_name],
                maybe_return)
        elif new_name in ['ARMCI_GetS',   'ARMCI_PutS',   'ARMCI_AccS',
                          'ARMCI_NbGetS', 'ARMCI_NbPutS', 'ARMCI_NbAccS']:
            print '''
%s
{
    %s
    armci_profile_start_strided(count, stride_levels, proc, %s);
    %s%s;
    armci_profile_stop_strided(%s);
    %s
}
''' % (func.get_signature(new_name),
                maybe_declare,
                name_to_event[new_name],
                maybe_assign, func.get_call(),
                name_to_event[new_name],
                maybe_return)
        elif new_name in ['ARMCI_Get',   'ARMCI_Put',   'ARMCI_Acc',
                          'ARMCI_NbGet', 'ARMCI_NbPut', 'ARMCI_NbAcc']:
            print '''
%s
{
    %s
    armci_profile_start_strided(&bytes, 0, proc, %s);
    %s%s;
    armci_profile_stop_strided(%s);
    %s
}
''' % (func.get_signature(new_name),
                maybe_declare,
                name_to_event[new_name],
                maybe_assign, func.get_call(),
                name_to_event[new_name],
                maybe_return)
        elif "ARMCI_Fence" in new_name:
            print '''
%s
{
    if (!SAMECLUSNODE(proc))
    armci_profile_start(ARMCI_PROF_FENCE);
    %s%s;
    if (!SAMECLUSNODE(proc))
    armci_profile_stop(ARMCI_PROF_FENCE);
}
''' % (func.get_signature(new_name), maybe_assign, func.get_call())
        elif new_name in name_to_event:
            print '''
%s
{
    %s
    armci_profile_start(%s);
    %s%s;
    armci_profile_stop(%s);
    %s
}
''' % (func.get_signature(new_name),
                maybe_declare,
                name_to_event[new_name],
                maybe_assign, func.get_call(),
                name_to_event[new_name],
                maybe_return)
        else:
            print '''
%s
{
    %s
    %s%s;
    %s
}
''' % (func.get_signature(new_name),
        maybe_declare,
        maybe_assign, func.get_call(),
        maybe_return)
