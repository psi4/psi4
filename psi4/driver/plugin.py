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

# flake8: noqa
import os
import sys

from psi4.driver.util.filesystem import *
from psi4.driver.util import tty
import psi4.config as config

def sanitize_name(name):
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
        tty.die("Plugin name must begin with a letter.")

# Determine the available plugins
available_plugins = []
plugin_path = join_path(config.share_dir, "plugin")
for dir in os.listdir(plugin_path):
    if os.path.isdir(join_path(plugin_path, dir)):
        available_plugins.append(dir)


# def create_new_plugin_makefile():
#     """Generate here (.) a plugin Makefile with current build settings."""
#
#     print("""Creating new plugin Makefile in the current directory.""")
#     file_manager(name='.', files={'Makefile': 'Makefile.template'})


def create_plugin(args):
    """Generate plugin in sanitized directory of same name based upon *type*"""

    name = sanitize_name(args['new_plugin'])
    type = args['new_plugin_template']
    template_path = join_path(plugin_path, type)

    # Create, but do not overwrite, plugin directory
    if os.path.exists(name):
        tty.error("""Plugin directory "{}" already exists.""".format(name))

    # Do a first pass to determine the template files
    template_files = os.listdir(template_path)
    source_files = []
    for file in template_files:
        target_file = file

        if file.endswith('.template'):
            target_file = file[0:-9]

        if file.endswith('.cc.template'):
            source_files.append(target_file)

    tty.hline("""Creating "{}" with "{}" template.""".format(name, type))

    os.mkdir(name)
    created_files = []
    for source_file in template_files:
        target_file = file

        if source_file.endswith('.template'):
            target_file = source_file[0:-9]

        try:
            with open(join_path(template_path, source_file), 'r') as file:
                contents = file.read()
        except IOError as err:
            tty.error("""Unable to open {} template.""".format(source_file))
            tty.error(err)
            sys.exit(1)

        contents = contents.replace('@plugin@', name)
        contents = contents.replace('@Plugin@', name.capitalize())
        contents = contents.replace('@PLUGIN@', name.upper())
        contents = contents.replace('@sources@', ' '.join(source_files))
        contents = contents.replace('@C@', config.c_compiler)
        contents = contents.replace('@CXX@', config.cxx_compiler)
        contents = contents.replace('@Fortran@', config.fortran_compiler)

        try:
            with open(join_path(name, target_file), 'w') as file:
                file.write(contents)
                created_files.append(target_file)
        except IOError as err:
            tty.error("""Unable to create {}""".format(target_file))
            tty.error(err)
            sys.exit(1)

    tty.info("Created plugin files: ", ", ".join(created_files))

    sys.exit(0)

