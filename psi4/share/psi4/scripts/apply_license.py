# Checks all psi4 relevant files for proper boilerplate GNU license.
# This is sold as is with no warrenty-- probably should double check everything
# after running. I am not responsible if you break Psi4.
#
# Do not forget to do share/plugins by hand!

import os

# File type we know how to handle
ftypes = ['cc', 'h', 'py']

c_header ="""/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */"""

py_header = c_header.replace(' */', '#')
py_header = py_header.replace('/*', '#')
py_header = py_header.replace(' *', '#')

c_header =  c_header.splitlines()
py_header = py_header.splitlines()


def check_header(infile):
    f = open(infile, 'r+')
    data = f.read().splitlines()

    # Find the header location
    max_lines = 30
    try:
        symbol = None
        if filename.split('.')[-1] in ['py']:
            start = data.index("# @BEGIN LICENSE") - 1
            end = data.index("# @END LICENSE") + 1
            if data[start] != '#' or data[end] != '#':
                f.close()
                print('Did not find "wings" of license block in file %s' % infile)
                return
        else:
            start = data.index(" * @BEGIN LICENSE") - 1
            end = data.index(" * @END LICENSE") + 1
            if data[start] != '/*' or data[end] != ' */':
                f.close()
                print('Did not find "wings" of license block in file %s' % infile)
                return
    except:
        print('Could not find license block in file %s' % infile)
        f.close()
        return

    # Make sure the block actually looks like a license
    license = data[start:end+1]
    top = any("PSI4:" in x.upper() for x in license[:5])
    bot = any("51 Franklin Street" in x for x in license[5:])
    if not (top and bot):
        print('Did not understand infile %s' % infile)
        f.close()
        return

    # Replace license
    if filename.split('.')[-1] in ['cc', 'h']:
        data[start:end + 1] = c_header 
    elif filename.split('.')[-1] in ['py']:
        data[start:end + 1] = py_header 
    else:
        print('Did not understand infile end: %s' % infile)
        f.close()
        return
   
    # Write it out 
    f.seek(0)
    f.write("\n".join(data))
    f.truncate()
    f.close()

avoid_strings = ['qcdb', 'libJKFactory']
    
walk = list(os.walk('../../src/'))
walk += list(os.walk('../python'))

for root, dirnames, filenames in walk:
    if any(x in root for x in avoid_strings):
        continue

    for filename in filenames:

        if filename.split('.')[-1] not in ftypes:
            continue

        check_header(root + '/' + filename)
