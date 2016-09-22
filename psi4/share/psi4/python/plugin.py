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

from __future__ import print_function
from __future__ import absolute_import
import os
import re
import sys
import glob


def sanitize_plugin_name(name):
    """Function to return *name* in coded form, stripped of
    characters that confuse filenames, characters into lowercase,
    ``+`` into ``p``, ``*`` into ``s``, and ``(``, ``)``, ``-``,
    & ``,`` into ``_``.

    """
    if name[0].isalpha():
        temp = name.lower()
        temp = temp.replace('+', 'p')
        temp = temp.replace('*', 's')
        temp = temp.replace('(', '_')
        temp = temp.replace(')', '_')
        temp = temp.replace(',', '_')
        temp = temp.replace('-', '_')
        return temp
    else:
        print("""Error: Plugin name must begin with a letter.""")
        sys.exit(1)


def get_template_path():
    """Return the PSIDATADIR/plugin directory where the plugin template files
    reside. Actually used is one directory up from location of this file.

    """
    thisdir = os.path.dirname(os.path.realpath(__file__))
    template_path = os.path.abspath(os.path.join(thisdir, '..', 'plugin'))
    return template_path


def create_new_plugin_makefile():
    """Generate here (.) a plugin Makefile with current build settings."""

    print("""Creating new plugin Makefile in the current directory.""")
    file_manager(name='.', files={'Makefile': 'Makefile.template'})


def create_new_plugin(name, template='plugin'):
    """Generate plugin *name* in sanitized directory of same name based
    upon *template*.

    """
    name = sanitize_plugin_name(name)
    ltemplate = template.lower()
    template_path = get_template_path()

    # select files for template
    template_files = {}
    ext = '.template'

    #   primary C++ file
    custom_fl = ltemplate + '.cc' + ext
    if os.path.isfile(os.path.join(template_path, custom_fl)):
        template_files[name + '.cc'] = custom_fl
    else:
        print("""Error: Template, {}, not found""".format(template))
        sys.exit(1)

    #   helper files
    for fl in ['__init__.py', 'doc.rst', 'input.dat', 'Makefile', 'pymodule.py']:
        custom_fl = ltemplate + '.' + fl + ext
        if os.path.isfile(os.path.join(template_path, custom_fl)):
            template_files[fl] = custom_fl
        else:
            template_files[fl] = fl + ext

    #   any extra files
    starting = os.path.join(template_path, ltemplate)
    for fl in glob.glob(starting + '.*.*' + ext):
        written_fl = fl[len(starting)+1:-len(ext)]
        template_files[written_fl] = fl

    # create, but not overwrite, plugin directory
    if os.path.exists(name):
        print("""Error: Plugin directory {} already exists.""".format(name))
        sys.exit(1)
    os.mkdir(name)
    print("""Created new plugin directory, {}, using '{}' template.""".format(
          name, ltemplate))
    file_manager(name=name, files=template_files)


def file_manager(name, files):
    """Process pairs of target files (to be written to path *name*) and
    source files (to be read from template library) in dictionary *files*
    through various string replacement macros.

    """
    template_path = get_template_path()

    for tgtfl, srcfl in files.items():
        source_name = os.path.join(template_path, srcfl)
        target_name = os.path.join(name, tgtfl)

        try:
            with open(source_name, 'r') as srchandle:
                srcflstr = srchandle.read()
        except IOError as err: 
            print("""Error: create_new_plugin: Unable to open {} template:""".format(source_name), err)
            sys.exit(1)

        # search and replace placeholders in the string
        srcflstr = srcflstr.replace('@plugin@', name)
        srcflstr = srcflstr.replace('@Plugin@', name.capitalize())
        srcflstr = srcflstr.replace('@PLUGIN@', name.upper())
        #srcflstr = srcflstr.replace('@PLUGIN_CXX@', format_cxx)
        #srcflstr = srcflstr.replace('@PLUGIN_DEFINES@', format_defines)
        #srcflstr = srcflstr.replace('@PLUGIN_FLAGS@', format_flags)
        #srcflstr = srcflstr.replace('@PLUGIN_INCLUDES@', format_includes)
        #srcflstr = srcflstr.replace('@PLUGIN_OBJDIR@', format_objdir)
        #srcflstr = srcflstr.replace('@PLUGIN_LDFLAGS@', format_ldflags)

        try:
            with open(target_name, 'w') as tgthandle:
                tgthandle.write(srcflstr)
                print("""    Created: {}""".format(tgtfl))
        except IOError as err:
            print("""Error: Unable to create {}:""".format(target_name), err)
            sys.exit(1)


if __name__ == "__main__":

    create_new_plugin_makefile()
    create_new_plugin('asdf')
    create_new_plugin('asdfscf', template='scf')
    create_new_plugin('asdf-fake', template='faKE')

