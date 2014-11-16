#!/usr/bin/python
# -*- python -*-
# -*- coding: utf-8 -*-
# vim:filetype=python:
# Create CMakeLists.txt template for leaf directories 
# (c) Roberto Di Remigio  <roberto.d.remigio@uit.no>
# licensed under the GNU Lesser General Public License

import sys
import glob

sys.path.append('../../cmake')
import argparse

parser = argparse.ArgumentParser(description='Create CMakeLists.txt template')
parser.add_argument('libname',
                     action='store',
                     metavar='LIBNAME',
                     nargs='?',
                     help='Name of the library to be created WITHOUT the lib prefix')
parser.add_argument('--lang',
        nargs='?',
        action='store',
        choices=('CXX', 'C', 'F'),
        default='CXX',
        const='CXX',
        help='Source file language')

args = parser.parse_args() 
libname = args.libname
lang    = args.lang

def glob_sources_cxx():
    """Create a list of C++ headers and sources to be used in a CMakeLists.txt file."""
    list_of_headers = glob.glob('*.h')
    headers = 'list(APPEND headers_list '
    for header in list_of_headers:
        headers += header + ' '
    headers += ')\n\n'
    list_of_sources = glob.glob('*.cc')
    sources = 'list(APPEND sources_list '
    for source in list_of_sources:
        sources += source + ' '
    sources += ')\n\n'

    message = 'set(headers_list "")\n'                                     \
    + '# List of headers\n' + headers                                      \
    + '# If you want to remove some headers specify them explictly here\n' \
    + 'if(DEVELOPMENT_CODE)\n'                                             \
    + '   list(REMOVE_ITEM headers_list "")\n'                             \
    + 'else()\n'                                                           \
    + '   list(REMOVE_ITEM headers_list "")\n'                             \
    + 'endif()\n'                                                          \
    + '# Sort alphabetically\n'                                            \
    + 'list(SORT headers_list)\n\n'                                        \
    + 'set(sources_list "")\n'                                             \
    + '# List of sources\n' + sources                                      \
    + '# If you want to remove some sources specify them explictly here\n' \
    + 'if(DEVELOPMENT_CODE)\n'                                             \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'else()\n'                                                           \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'endif()\n\n'                                                                 
    return message            

def glob_sources_c():
    """Create a list of C headers and sources to be used in a CMakeLists.txt file."""
    list_of_headers = glob.glob('*.h')
    headers = 'list(APPEND headers_list '
    for header in list_of_headers:
        headers += header + ' '
    headers += ')\n\n'
    list_of_sources = glob.glob('*.c')
    sources = 'list(APPEND sources_list '
    for source in list_of_sources:
        sources += source + ' '
    sources += ')\n\n'

    message = 'set(headers_list "")\n'                                     \
    + '# List of headers\n' + headers                                      \
    + '# If you want to remove some headers specify them explictly here\n' \
    + 'if(DEVELOPMENT_CODE)\n'                                             \
    + '   list(REMOVE_ITEM headers_list "")\n'                             \
    + 'else()\n'                                                           \
    + '   list(REMOVE_ITEM headers_list "")\n'                             \
    + 'endif()\n'                                                          \
    + '# Sort alphabetically\n'                                            \
    + 'list(SORT headers_list)\n\n'                                        \
    + 'set(sources_list "")\n'                                             \
    + '# List of sources\n' + sources                                      \
    + '# If you want to remove some sources specify them explictly here\n' \
    + 'if(DEVELOPMENT_CODE)\n'                                             \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'else()\n'                                                           \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'endif()\n\n'                                                                 
    return message            

def glob_sources_fortran():
    """Create a list of Fortran sources to be used in a CMakeLists.txt file."""
    list_of_sources = glob.glob('*.f')
    list_of_sources.extend(glob.glob('*.F'))
    list_of_sources.extend(glob.glob('*.f77'))
    list_of_sources.extend(glob.glob('*.F77'))
    list_of_sources.extend(glob.glob('*.f90'))
    list_of_sources.extend(glob.glob('*.F90'))

    sources = 'list(APPEND sources_list '
    for source in list_of_sources:
        sources += source + ' '
    sources += ')\n\n'

    message = '# List of headers (needed to avoid an "empty list error" from CMake)\n' \
    + 'set(headers_list "")\n\n'                                           \
    + 'set(sources_list "")\n'                                             \
    + '# List of sources\n' + sources                                      \
    + '# If you want to remove some sources specify them explictly here\n' \
    + 'if(DEVELOPMENT_CODE)\n'                                             \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'else()\n'                                                           \
    + '   list(REMOVE_ITEM sources_list "")\n'                             \
    + 'endif()\n\n'                                                                 
    return message            

f = open('CMakeLists.txt.try', 'w')
if (lang == 'CXX'):
    f.write(glob_sources_cxx())
elif (lang == 'C'):
    f.write(glob_sources_c())
else:
    f.write(glob_sources_fortran())

f.write('# Write list of files to be passed to cloc for counting lines of code.\n')
f.write('# Only files that are actually compiled are counted.\n')
f.write('set(to_count "${sources_list}" "${headers_list}")\n')
f.write('write_to_cloc_list("${to_count}")\n\n')

f.write('# Build static library\n')
f.write('add_library('+ libname + ' STATIC ${sources_list})\n')
f.write('# Specify dependencies for the library (if any)\n')
f.write('#add_dependencies('+ libname + ' )\n')
f.write('set(libs_to_merge '+ libname + ' ${libs_to_merge} PARENT_SCOPE)\n')
if (not (lang == 'C' or lang == 'F')):
    f.write('if(BUILD_CUSTOM_BOOST)\n')
    f.write('   add_dependencies('+ libname + ' custom_boost)\n')
    f.write('endif()\n')
f.write('install(TARGETS ' + libname + ' ARCHIVE DESTINATION lib)\n\n')

if (not lang == 'F'):
    f.write('# Sets install directory for all the headers in the list\n')
    f.write('install_list_FILES("${headers_list}" include/lib' + libname + ')\n')

print('Template for {} created'.format(args.libname))
print('Don\'t forget to fix excluded files and dependencies!!!')

# vim:et:ts=4:sw=4
