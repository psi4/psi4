#!/usr/bin/env python

'''Generate armci's capi.c source from the parmci.h header.'''

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
        print 'usage: capi_gen.py <parmci.h> > <capi.c>'
        sys.exit(len(sys.argv))

    # print headers
    print '''
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include "armci.h"
#include "parmci.h"
'''

    functions = {}
    # parse signatures into the Function class
    for sig in get_signatures(sys.argv[1]):
        function = Function(sig)
        functions[function.name] = function

    # now process the functions
    for name in sorted(functions):
        func = functions[name]
        maybe_return = ''
        if '*' in func.return_type or 'void' not in func.return_type:
            maybe_return = 'return '
        func = functions[name]
        new_name = None
        if 'PARMCI_' in name:
            new_name = name.replace('PARMCI_','ARMCI_')
        elif 'parmci_' in name:
            new_name = name.replace('parmci_','armci_')
        print '''
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak %s
#endif
%s
{
    %s%s;
}
''' % (new_name,
        func.get_signature(new_name),
        maybe_return, func.get_call())
