#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
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

from psi4.driver.util.filesystem import *
from psi4.driver.util import tty


def sanitize_name(name):
    """Function to return *name* in coded form, stripped of
    characters that confuse filenames, characters into lowercase,
    ``+`` into ``p``, ``*`` into ``s``, and ``(``, ``)``, ``-``,
    & ``,`` into ``_``.

    Also checks the sanitized name against a list of restricted C++ keywords.
    """
    if name[0].isalpha():
        temp = name.lower()
        temp = temp.replace('+', 'p')
        temp = temp.replace('*', 's')
        temp = temp.replace('(', '_')
        temp = temp.replace(')', '_')
        temp = temp.replace(',', '_')
        temp = temp.replace('-', '_')

        # Taken from http://en.cppreference.com/w/cpp/keyword
        cpp_keywords = [
            "alignas", "alignof", "and", "and_eq", "asm", "atomic_cancel",
            "atomic_commit", "atomic_noexcept", "auto", "bitand", "bitor",
            "bool", "break", "case", "catch", "char", "char16_t", "char32_t",
            "class", "compl", "concept", "const", "constexpr", "const_cast",
            "continue", "decltype", "default", "delete", "do", "double",
            "dynamic_cast", "else", "enum", "explicit", "export", "extern",
            "false", "float", "for", "friend", "goto", "if", "import", "inline",
            "int", "long", "module", "mutable", "namespace", "new", "noexcept",
            "not", "not_eq", "nullptr", "operator", "or", "or_eq", "private",
            "protected", "public", "register", "reinterpret_cast", "requires",
            "return", "short", "signed", "sizeof", "static", "static_assert",
            "static_cast", "struct", "switch", "synchronized", "template",
            "this", "thread_local", "throw", "true", "try", "typedef", "typeid",
            "typename", "union", "unsigned", "using", "virtual", "void",
            "volatile", "wchar_t", "while", "xor", "xor_eq",

            # Identifiers with special meanings"
            "override", "final", "transaction_safe", "transaction_safe_dynamic",

            # Preprocessor tokens
            "if", "elif", "else", "endif", "defined", "ifdef", "ifndef",
            "define", "undef", "include", "line", "error", "pragma",
            "_pragma"
        ]

        if temp in cpp_keywords:
            tty.die("The plugin name you provided is a C++ reserved keyword.  Please provide a different name.")

        return temp
    else:
        tty.die("Plugin name must begin with a letter.")


# Determine the available plugins
available_plugins = []
psidatadir = os.environ.get('PSIDATADIR', None)
plugin_path = join_path(psidatadir, "plugin")
for dir in os.listdir(plugin_path):
    if os.path.isdir(join_path(plugin_path, dir)):
        available_plugins.append(dir)


def create_plugin(name, template):
    """Generate plugin in directory with sanitized *name* based upon *template*."""

    name = sanitize_name(name)
    template_path = join_path(plugin_path, template)

    # Create, but do not overwrite, plugin directory
    if os.path.exists(name):
        tty.error("""Plugin directory "{}" already exists.""".format(name))

    # Do a first pass to determine the template temp_files
    template_files = os.listdir(template_path)
    source_files = []
    for temp_file in template_files:
        target_file = temp_file

        if temp_file.endswith('.template'):
            target_file = temp_file[0:-9]

        if temp_file.endswith('.cc.template'):
            source_files.append(target_file)

    tty.hline("""Creating "{}" with "{}" template.""".format(name, template))

    os.mkdir(name)
    created_files = []
    for source_file in template_files:
        target_file = source_file

        if source_file.endswith('.template'):
            target_file = source_file[0:-9]

        try:
            with open(join_path(template_path, source_file), 'r') as temp_file:
                contents = temp_file.read()
        except IOError as err:
            tty.error("""Unable to open {} template.""".format(source_file))
            tty.error(err)
            sys.exit(1)

        contents = contents.replace('@plugin@', name)
        contents = contents.replace('@Plugin@', name.capitalize())
        contents = contents.replace('@PLUGIN@', name.upper())
        contents = contents.replace('@sources@', ' '.join(source_files))

        try:
            with open(join_path(name, target_file), 'w') as temp_file:
                temp_file.write(contents)
                created_files.append(target_file)
        except IOError as err:
            tty.error("""Unable to create {}""".format(target_file))
            tty.error(err)
            sys.exit(1)

    tty.info("Created plugin files (in {} as {}): ".format(name, template), ", ".join(created_files))

    sys.exit(0)
