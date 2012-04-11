#!/usr/bin/python

import sys
import os
import glob
import re


DriverPath = ''
InsertPath = '/../../../'
if (len(sys.argv) == 2):
    DriverPath = sys.argv[1] + '/'
    sys.path.insert(0, os.path.abspath(os.getcwd()))
    import apply_relpath
    IncludePath = apply_relpath.get_topsrcdir_asrelativepathto_objdirsfnxsource()[1]


def pts(category, pyfile):
    print 'Auto-documenting %s file %s' % (category, pyfile)


# helper fn
def sphinxify_comment(text):

    text = text.replace('@@', '_')
    text = text.replace(' $', ' :math:`')
    text = text.replace('($', '(\ :math:`')
    text = text.replace('$ ', '` ')
    text = text.replace('$.', '`.')
    text = text.replace('$,', '`,')
    text = text.replace('$)', '`\ )')

    return text


# helper fn
# including the options abbr substitutions file in every SSSOUT option file slows
#   compilation by a factor of ten. so, back-translate |%s__%s| into :term:`%s`
def substitute_comment(cmnt):
    subst = re.compile(r'^(.*?)[\s\(]\|(\w+)__(\w+)\|[\s\).,](.*?)$')

    while True:
        if subst.match(cmnt):
            m = subst.match(cmnt)
            cmnt = m.group(1) + ' :term:`' + m.group(3).upper() + ' <' + m.group(3).upper() + ' (' + m.group(2).upper() + ')>` ' + m.group(4)
            continue
        else:
            break

    return cmnt


# helper fn
def determine_options(cfilename):

    module = re.compile(r'^(.*)name\s*==\s*"(.*)"(.*?)$', re.IGNORECASE)
    modulecomment = re.compile(r'^(\s*?)\/\*-\s*MODULEDESCRIPTION\s*(.*?)-\*\/(\s*?)$', re.IGNORECASE)
    modulecommentstart = re.compile(r'^(\s*?)\/\*-\s*MODULEDESCRIPTION\s*(.*?)(\s*?)$', re.IGNORECASE)
    subsection = re.compile(r'^(\s*?)\/\*-\s*SUBSECTION\s*(.*?)\s*-\*\/(\s*?)$', re.IGNORECASE)
    comment = re.compile(r'^(\s*?)\/\*-\s*(.*?)-\*\/(\s*?)$', re.IGNORECASE)
    commentend = re.compile(r'^(\s*)(.*?)-\*\/(\s*?)$', re.IGNORECASE)
    commentstart = re.compile(r'^(\s*?)\/\*-\s*(.*)(\s*?)$', re.IGNORECASE)
    kw_string_def_opt = re.compile(r'add_str\(\s*"(.*)"\s*,\s*"(.*)"\s*,\s*"(.*)"\s*\)')
    kw_string_def_opt_2 = re.compile(r'add_str_i\(\s*"(.*)"\s*,\s*"(.*)"\s*,\s*"(.*)"\s*\)')
    kw_string_def = re.compile(r'add_str\(\s*"(.*)"\s*,\s*"(.*)"\s*\)')
    kw_string_def_2 = re.compile(r'add_str_i\(\s*"(.*)"\s*,\s*"(.*)"\s*\)')
    kw_bool_def = re.compile(r'add_bool\(\s*"(.*)"\s*,\s*("?)([-\w]+)("?)\s*\)')
    kw_double_def = re.compile(r'add_double\(\s*"(.*)"\s*,\s*("?)([-/\.\w]+)("?)\s*\)')
    kw_generic_def = re.compile(r'add_(\w+)\(\s*"(\w+)"\s*,\s*("?)([-\w]+)("?)\s*\)')  # untested
    kw_complicated = re.compile(r'add\(\s*"(\w*)"\s*,\s*new\s+(\w+)\(\)\s*\)')  # untested

    fcfile = open(cfilename)
    contents = fcfile.readlines()
    fcfile.close()

    ii = 0
    while (ii < len(contents)):
        line = contents[ii]

        if module.match(line):
            currentmodule = module.match(line).group(2).upper()
            fmodule.write('.. toctree::\n   :hidden:\n   :glob:\n\n   %s__*\n\n' % (currentmodule.lower()))

        elif modulecommentstart.match(line):
            tag = ''
            while 1:
                if (not commentend.match(line)):
                    if modulecommentstart.match(line):
                        tag += modulecommentstart.match(line).group(2)
                    else:
                        tag += ' ' + line.strip()
                    ii += 1
                    line = contents[ii]
                    continue
                else:
                    if modulecomment.match(line):
                        tag += modulecomment.match(line).group(2)
                        break
                    else:
                        tag += ' ' + commentend.match(line).group(2)
                        break
            fglossary.write('**%s**: %s\n\n' % (currentmodule, tag))

        elif subsection.match(line):
            currentsubsection = subsection.match(line).group(2)
            fglossary.write('\n%s\n%s\n\n' % (currentsubsection, '^' * len(currentsubsection)))
            fglossary.write('.. glossary::\n   :sorted:\n\n')

        elif commentstart.match(line):
            tag = ''
            while 1:
                if (not commentend.match(line)):
                    if commentstart.match(line):
                        tag += commentstart.match(line).group(2)
                    else:
                        tag += ' ' + line.strip()
                    ii += 1
                    line = contents[ii]
                    continue
                else:
                    if comment.match(line):
                        tag += comment.match(line).group(2)
                        break
                    else:
                        tag += ' ' + commentend.match(line).group(2)
                        break
            tag = sphinxify_comment(tag)
            # capture option immediately after comment
            kw_name = ''
            kw_default = 'No Default'
            kw_type = ''
            kw_possible = ''
            ii += 1
            line = contents[ii]
            if (not line or line.isspace()):
                ii += 1
                line = contents[ii]

            if kw_string_def_opt.search(line):
                m = kw_string_def_opt.search(line)
                kw_name = m.group(1)
                kw_type = 'str'
                if not (not m.group(2) or m.group(2).isspace()):
                    kw_default = m.group(2)
                kw_possible = m.group(3)
            elif kw_string_def_opt_2.search(line):
                m = kw_string_def_opt_2.search(line)
                kw_name = m.group(1)
                kw_type = 'str'
                if not (not m.group(2) or m.group(2).isspace()):
                    kw_default = m.group(2)
                kw_possible = m.group(3)
            elif kw_string_def.search(line):
                m = kw_string_def.search(line)
                kw_name = m.group(1)
                kw_type = 'str'
                if not (not m.group(2) or m.group(2).isspace()):
                    kw_default = m.group(2)
            elif kw_string_def_2.search(line):
                m = kw_string_def_2.search(line)
                kw_name = m.group(1)
                kw_type = 'str'
                if not (not m.group(2) or m.group(2).isspace()):
                    kw_default = m.group(2)
            elif kw_bool_def.search(line):
                m = kw_bool_def.search(line)
                kw_name = m.group(1)
                kw_type = 'bool'
                if not (not m.group(3) or m.group(3).isspace()):
                    kw_default = m.group(3).lower()
                    if kw_default == '1':
                        kw_default = 'true'
                    if kw_default == '0':
                        kw_default = 'false'
            elif kw_double_def.search(line):
                m = kw_double_def.search(line)
                kw_name = m.group(1)
                kw_type = 'double'
                if not (not m.group(3) or m.group(3).isspace()):
                    kw_default = m.group(3).lower()
            elif kw_generic_def.search(line):
                m = kw_generic_def.search(line)
                kw_name = m.group(2)
                kw_type = m.group(1)
                if not (not m.group(4) or m.group(4).isspace()):
                    kw_default = m.group(4).lower()
            elif kw_complicated.search(line):
                m = kw_complicated.search(line)
                kw_name = m.group(1)
                kw_type = m.group(2)
                if kw_type == 'ArrayType':
                    kw_type = 'array'
                elif kw_type == 'MapType':
                    kw_type = 'map'
                elif kw_type == 'PythonDataType':
                    kw_type = 'python'
                else:
                    print 'ERROR: unrecognized type %s for %s' % (kw_type, kw_name)
                    sys.exit()

            if   kw_type == 'str':    kw_type = 'string'
            elif kw_type == 'int':    kw_type = 'integer'
            elif kw_type == 'bool':   kw_type = 'boolean'
            elif kw_type == 'double': pass
            elif kw_type == 'array':  pass
            elif kw_type == 'map':    pass
            elif kw_type == 'python': pass
            else:
                print 'ERROR: unrecognized type2 %s for %s' % (kw_type, kw_name)
                sys.exit()

            #print 'kw_name = \t', kw_name
            #print 'kw_type = \t', kw_type
            #print 'kw_dflt = \t', kw_default
            #print 'kw_poss = \t', kw_possible
            #print 'kw_tagl = \t', tag
            #print '\n'

            # substitution list file
            fabbr.write('.. |%s__%s| replace:: :term:`%s <%s (%s)>`\n' % 
              (currentmodule.lower(), kw_name.lower(), kw_name.upper(), kw_name.upper(), currentmodule.upper()))

            # individual option file for plugin options. rather pointless but consistent w/regular module options
            fsssdoc = open('source/autodir_plugins/'+currentmodule.lower()+'__'+kw_name.lower()+'.rst', 'w')
            div = '"' * (14 + len(currentmodule) + 2 * len(kw_name))
            fsssdoc.write(':term:`%s <%s (%s)>`\n%s\n\n' % (kw_name.upper(), kw_name.upper(), currentmodule.upper(), div))
            fsssdoc.write('      %s\n\n' % (substitute_comment(tag)))

            fglossary.write('   %s (%s)\n      %s\n\n' % (kw_name.upper(), currentmodule.upper(), tag))

            if kw_type == 'boolean':
                fglossary.write('      * **Type**: :ref:`boolean <op_c_boolean>`\n')
                fsssdoc.write('      * **Type**: :ref:`boolean <op_c_boolean>`\n')

            elif (kw_type == 'double') and ((kw_name.lower().find('conv') > -1) or (kw_name.lower().find('tol') > -1)):
                fglossary.write('      * **Type**: :ref:`conv double <op_c_conv>`\n')
                fsssdoc.write('      * **Type**: :ref:`conv double <op_c_conv>`\n')

            elif (kw_type == 'string') and ((kw_name.lower() == 'basis') or (kw_name.lower().startswith('df_basis'))):
                fglossary.write('      * **Type**: %s\n' % kw_type)
                fsssdoc.write('      * **Type**: %s\n' % kw_type)
                fglossary.write('      * **Possible Values**: :ref:`basis string <apdx:basisElement>`\n')
                fsssdoc.write('      * **Possible Values**: :ref:`basis string <apdx:basisElement>`\n')

            else:
                fglossary.write('      * **Type**: %s\n' % kw_type)
                fsssdoc.write('      * **Type**: %s\n' % kw_type)

            if not (not kw_possible or kw_possible.isspace()):
                sline = kw_possible.split()
                fglossary.write('      * **Possible Values**: %s\n' % (', '.join(sline)))
                fsssdoc.write('      * **Possible Values**: %s\n' % (', '.join(sline)))

            fglossary.write('      * **Default**: %s\n\n' % kw_default)
            fsssdoc.write('      * **Default**: %s\n\n' % kw_default)
            fsssdoc.close()

        if (line.find('extern "C" PsiReturnType') > -1):
            break

        ii += 1


# Objective #3
# Plugin directories in psi4/tests/plugin_
fdriver = open('source/autodoc_available_plugins.rst', 'w')
fdriver.write('\n.. index:: plugins; available\n')
fdriver.write('.. _`sec:availablePlugins`:\n\n')
fdriver.write('====================================================\n')
fdriver.write('Emerging Theoretical Methods: Plugins DFADC to RQCHF\n')
fdriver.write('====================================================\n\n')
fdriver.write('.. toctree::\n   :maxdepth: 1\n\n')

fabbr = open('source/autodoc_abbr_options_plugins.rst', 'w')

# from each plugin directory ...
for pydir in glob.glob(DriverPath + '../../tests/plugin_*'):
    dirname = os.path.split(pydir)[1]
    div = '=' * len(dirname)

    if dirname not in []:

        pts('plugin', dirname)
        fdriver.write('   autodir_plugins/module__%s' % (dirname))

        fmodule = open('source/autodir_plugins/module__'+dirname+'.rst', 'w')
        fmodule.write('\n.. _`sec:%s`:\n' % (dirname.lower()))
        fmodule.write('.. index:: plugin; %s\n\n' % (dirname.lower()))
        fmodule.write(':srcplugin:`' + dirname.lower() + '`\n')
        fmodule.write(div + '=============' + '\n\n')
        #fmodule.write(dirname.lower() + '\n')
        #fmodule.write(div + '\n\n')
        #fmodule.write('.. toctree::\n   :hidden:\n   :glob:\n\n   %s__*\n\n' % (dirname.lower()))
        fmodule.write('.. toctree::\n   :hidden:\n\n   /autodir_plugins/glossary__%s\n\n' % (dirname.lower()))
        fglossary = open('source/autodir_plugins/glossary__'+dirname+'.rst', 'w')
        fglossary.write('\n.. include:: /autodoc_abbr_options_c.rst\n')
        fglossary.write('.. include:: /autodoc_abbr_options_plugins.rst\n\n')
        fglossary.write('.. glossary::\n   :sorted:\n\n')

        # ... include doc.rst file
        docfile = '%s/doc.rst' % (pydir)
        if os.path.isfile(docfile):
            fmodule.write('.. include:: %stests/%s/doc.rst\n\n' % (IncludePath, dirname))

        # ... include docstrings from any *.py files
        pyfiles = glob.glob(pydir + '/*.py')
        if len(pyfiles) > 0:
            fmodule.write('Py-side Documentation\n')
            fmodule.write('---------------------\n\n')
            for pyfile in pyfiles:
                filename = os.path.split(pyfile)[1]
                basename = os.path.splitext(filename)[0]
                fmodule.write('.. automodule:: %s.%s\n' % (dirname, basename))
                fmodule.write('   :members:\n')
                fmodule.write('   :undoc-members:\n\n')

        # ... include keywords section from any *.cc files
              # todo: turn this into a fn and store in a dictionary
        cfiles = glob.glob(pydir + '/*.cc') + glob.glob(pydir + '/*.cc.in')
        if len(cfiles) > 0:
            fmodule.write('C-side Documentation\n')
            fmodule.write('--------------------\n\n')


            for cfile in cfiles:
                determine_options(cfile)

            fmodule.write('.. include:: /autodir_plugins/glossary__%s.rst' % (dirname))
        fmodule.write('\n\n')
        fmodule.close()
        fglossary.write('\n\n')
        fglossary.close()
    fdriver.write('\n')
fdriver.write('\n')
fdriver.close()
fabbr.write('\n')
fabbr.close()
