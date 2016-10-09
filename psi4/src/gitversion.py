#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import os
import sys
import subprocess

top_srcdir = './'
if len(sys.argv) == 2:
    top_srcdir = sys.argv[1] + '/'


def write_version(branch, mmp, ghash, status):
    version_str = '#define GIT_VERSION "{%s} %s %s"\n' % \
                  (branch, ghash, status)

    mmp_str = '#define PSI_VERSION "%s"\n' % \
              (mmp if mmp else '(no tag)')

    with open('gitversion.h.tmp', 'w') as handle:
        handle.write(version_str)
        handle.write(mmp_str)

    with open('../psi4-config.tmp', 'r') as handle:
        f = handle.read()
    with open('../psi4-config', 'w') as handle:
        handle.write(f)
        handle.write('    psiver = "%s"\n' % (mmp))
        handle.write('    githash = "{%s} %s %s"\n' % (branch, ghash, status))
        handle.write('    sys.exit(main(sys.argv))\n\n')
    os.chmod('../psi4-config', 0o755)


# Use Git to find current branch name
#   Returns "refs/heads/BRANCHNAME"
#   Use [11:] to skip refs/heads/
try:
    command = "git symbolic-ref -q HEAD"
    process = subprocess.Popen(command.split(), stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE, cwd=top_srcdir,
                               universal_newlines=True)
    (out, err) = process.communicate()
    branch = str(out).rstrip()[11:]
    if process.returncode:
        branch = "detached?"
except:
    branch = "detached?"

# Use Git to find a sortable latest version number
try:
    command = "git describe --long --dirty --always"
    process = subprocess.Popen(command.split(), stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE, cwd=top_srcdir,
                               universal_newlines=True)
    (out, err) = process.communicate()
    fields = str(out).rstrip().split('-')

    #         a68d223        # tags not pulled, clean git directory
    #         a68d223-dirty  # tags not pulled, changes to git-controlled files
    # 0.1-62-ga68d223        # tags pulled, clean git directory
    # 0.1-62-ga68d223-dirty  # tags pulled, changes to git-controlled files

    if fields[-1] == 'dirty':
        status = fields.pop()
    else:
        status = ''

    if len(fields[-1]) == 7:
        ghash = fields.pop()
    elif len(fields[-1]) == 8:
        ghash = fields.pop()[1:]
    else:
        ghash = ''

    # major-minor-patch
    if len(fields) == 2:
        if fields[0].endswith('rc'):
            mmp = ''.join(fields)
        elif fields[0] == '1.0rc2':
            # special case. choose offset to allow monotonic versioning
            mmp = '1.0rc' + str(200 + int(fields[1]))
        else:
            mmp = '.'.join(fields)
    else:
        mmp = ''

    if process.returncode:
        try:
            # try to get some minimal version info from tarball not under git control
            zipname = top_srcdir.split('/')[-2]  # ending slash guaranteed
            command = "unzip -z ../" + zipname
            process = subprocess.Popen(command.split(), stderr=subprocess.PIPE,
                                       stdout=subprocess.PIPE, cwd=top_srcdir,
                                       universal_newlines=True)
            (out, err) = process.communicate()
            fields = str(out).rstrip().split()

            if zipname.endswith('master'):
                branch = 'master'
            status = ''
            ghash = fields.pop()[:7]
            mmp = ''

        except:
            status = ''
            ghash = ''
            mmp = ''

except:
    status = ''
    ghash = ''
    mmp = ''

write_version(branch, mmp, ghash, status)
