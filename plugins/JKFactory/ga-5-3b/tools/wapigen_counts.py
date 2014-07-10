#!/usr/bin/env python

'''Generate the wapi_counts.c source from the ga-papi.h header.'''

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
        print 'usage: wapigen_counts.py <ga-papi.h> > <wapi_counts.c>'
        sys.exit(len(sys.argv))

    # print headers
    print '''
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga-papi.h"
#include "typesf2c.h"
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

    # now process the functions
    for name in sorted(functions):
        func = functions[name]
        if 'terminate' in name:
            continue
        maybe_return = ''
        if 'void' not in func.return_type:
            maybe_return = 'return '
        func = functions[name]
        wnga_name = name.replace('pnga_','wnga_')
        print '''
%s
{
    ++count_%s;
    %s%s;
}
''' % (func.get_signature(wnga_name), name, maybe_return, func.get_call())

    # prepare to output the terminate function
    name = 'pnga_terminate'
    wnga_name = name.replace('pnga_','wnga_')
    func = functions[name]
    the_code = ''
    # establish 'the_code' to use in the body of terminate
    # it's printing the count of each function if it was called at least once
    for name in sorted(functions):
        the_code += '''
        if (count_%s) {
            printf("%s %%ld\\n", count_%s);
        }
''' % (name, name, name)
    # output the terminate function
    print '''%s
{
    ++count_pnga_terminate;
    %s;
    /* don't dump info if terminate more than once */
    if (1 == count_pnga_terminate) {
%s
    }
}
''' % (func.get_signature(wnga_name), func.get_call(), the_code) 
