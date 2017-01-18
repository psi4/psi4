#! /usr/bin/env python

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
from __future__ import print_function

import os
import re
import sys
import argparse
import subprocess


#    with open('../psi4-config.tmp', 'r') as handle:
#        f = handle.read()
#    with open('../psi4-config', 'w') as handle:
#        handle.write(f)
#        handle.write('    psiver = "%s"\n' % (mmp))
#        handle.write('    githash = "{%s} %s %s"\n' % (branch, ghash, status))
#        handle.write('    sys.exit(main(sys.argv))\n\n')
#    os.chmod('../psi4-config', 0o755)


def collect_version_input_from_fallback(meta_file='metadata.py'):
    """From *meta_file*, collect lines matching ``_version_{key} = {value}``
    and return as dictionary.

    """
    cwd = os.path.dirname(os.path.abspath(__file__))
    res = dict(re.findall("__version_([a-z_]+)\s*=\s*'([^']+)'", open(cwd + '/' + meta_file).read()))
    res.pop('_')
    return res


def is_git_repo(cwd='./', dot_git_qualifies=False, no_git_cmd_result=False):
    """Returns boolean as to whether *cwd* is under git control. When no ``git``
    command available in environment, *no_git_cmd_result* returned. If within
    the .git directory of a git repository, *dot_git_qualifies* returned.

    """
    command = 'git rev-parse --is-inside-work-tree'
    try:
        process = subprocess.Popen(command.split(),
                                   stderr=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   cwd=cwd,
                                   universal_newlines=True)
    except EnvironmentError as e:
        # most likely, git command not available
        return no_git_cmd_result

    (out, err) = process.communicate()

    if process.returncode != 0:
        # fatal: Not a git repository (or any of the parent directories): .git
        return False

    if out.strip() == 'true':
        # in a git repo and not within .git dir
        return True

    if out.strip() == 'false':
        # in a git repo in .git dir
        return dot_git_qualifies


def collect_version_input_from_git():
    """Returns a dictionary filled with ``git describe`` results, clean/dirty
    flag, and branch status. *cwd* should already be confirmed as a git
    repository; this doesn't catch returncodes or EnvironmentErrors because the
    raised errors are preferred to incomplete return dictionary.

    """
    cwd = os.path.dirname(os.path.abspath(__file__))
    res = {}

    # * only want annotated tags, so not --all
    # * in case *no* tags (impossible in Psi4), --always gets at least hash
    # * get commits & hash info even if on tag using --long
    command = 'git describe --abbrev=7 --long --always HEAD'
    process = subprocess.Popen(command.split(),
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               cwd=cwd,
                               universal_newlines=True)
    (out, err) = process.communicate()

    fields = str(out).rstrip().split('-')
    if len(fields) == 3:
        # normal: 0.1-62-ga68d223
        res['latest_annotated_v_tag'] = fields[0][1:]  # drop the "v"; tag mismatch caught later
        res['commits_since_tag'] = fields[1]
        res['seven_char_hash'] = fields[2][1:]  # drop the "g" git identifier
    else:
        # no tag present: a68d223
        res['latest_annotated_v_tag'] = ''
        res['commits_since_tag'] = ''
        res['seven_char_hash'] = fields[0]  # no prepended "g"

    command = 'git diff-index --name-only HEAD'
    process = subprocess.Popen(command.split(),
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               cwd=cwd,
                               universal_newlines=True)
    (out, err) = process.communicate()

    res['is_clean'] = False if str(out).rstrip() else True

    command = 'git rev-parse --abbrev-ref HEAD'  # returns HEAD when detached
    process = subprocess.Popen(command.split(),
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               cwd=cwd,
                               universal_newlines=True)
    (out, err) = process.communicate()

    res['branch_name'] = str(out).rstrip()

    return res


def reconcile_and_compute_version_output(quiet=False):
    res = collect_version_input_from_fallback(meta_file='metadata.py')
    meta_latest_annotated_v_tag, _, meta_seven_char_hash = res['long'].partition('+')

    # this is the tag format (PEP440 compliant) that our machinery is expecting.
    #   let's catch any deviations with Travis before it can corrupt versioning.
    sane_tag = re.compile("""^(?P<tag>(?P<forwardseries>\d+\.\d+(?P<patch>\.[1-9]+)?)(?(patch)|(?P<prere>((a)|(b)|(rc))\d+)?))$""")

    mobj = sane_tag.match(meta_latest_annotated_v_tag)
    if mobj:
        # some versioning machinery (looking at you, CMake) does strictly
        #   numerical comparisons such as M.m.p.t and thus can't handle
        #   prereleases and dev snapshots. we compute a Most Rescent Ancestral
        #   Release tag (e.g., 1.0 or 1.12.1) for a backward release series.
        backwardseries = mobj.group('forwardseries')
        if mobj.group('prere'):
            tmp = backwardseries.split('.')
            bumpdown = str(int(tmp[-1]) - 1)
            if bumpdown == '-1':
                print("""Unavoidable snag. Probably "2.0". Can't predict backward series from present prerelease.""")
                sys.exit()
            else:
                tmp[-1] = bumpdown
                backwardseries = '.'.join(tmp)
    else:
        print("""Tag in {} is malformed: {}""".format(
            'metadata.py', meta_latest_annotated_v_tag))
        sys.exit()

    cwd = os.path.dirname(os.path.abspath(__file__))
    if is_git_repo(cwd=cwd):
        res.update(collect_version_input_from_git())

        # establish the default response
        project_release = False
        project_prerelease = False
        project_version = 'undefined'
        project_version_long = 'undefined+' + res['seven_char_hash']

        if res['latest_annotated_v_tag'] == meta_latest_annotated_v_tag:

            trial_version_long_release = res['latest_annotated_v_tag'] + '+' + res['seven_char_hash']
            trial_version_devel = res['upcoming_annotated_v_tag'] + '.dev' + res['commits_since_tag']
            trial_version_long_devel = trial_version_devel + '+' + res['seven_char_hash']

            if int(res['commits_since_tag']) == 0:

                if trial_version_long_release == res['long']:
                    print("""Amazing, this can't actually happen that git hash stored at git commit.""")
                    sys.exit()
                else:
                    if meta_seven_char_hash == 'zzzzzzz':
                        if not quiet:
                            print("""Defining {} version: {} (recorded and computed)""".format(
                                'prerelease' if mobj.group('prere') else 'release', trial_version_long_release))
                        project_release = res['is_clean'] and not mobj.group('prere')
                        project_prerelease = res['is_clean'] and mobj.group('prere')
                        project_version = meta_latest_annotated_v_tag
                        project_version_long = trial_version_long_release

                    else:
                        print("""Undefining version for irreconcilable hashes: {} (computed) vs {} (recorded)""".format(
                            trial_version_long_release, res['long']))

            else:
                if res['branch_name'].endswith('.x'):
                    print("""Undefining version as development snapshots not allowed on maintenance branch: {} (rejected computed)""".format(
                        trial_version_long_devel))

                # TODO prob should be undef unless on master
                else:
                    if not quiet:
                        print("""Defining development snapshot version: {} (computed)""".format(
                            trial_version_long_devel))
                    project_version = trial_version_devel
                    project_version_long = trial_version_long_devel

        else:
            print("""Undefining version for irreconcilable tags: {} (computed) vs {} (recorded)""".format(
                res['latest_annotated_v_tag'], meta_latest_annotated_v_tag))

    else:
        print("""Blindly (no git) accepting release version: {} (recorded)""".format(
            res['long']))
        # assumes that zip only comes from [pre]release. GitHub hides others, but they're there.
        project_release = not bool(mobj.group('prere'))
        project_prerelease = bool(mobj.group('prere'))
        project_version = meta_latest_annotated_v_tag
        project_version_long = res['long']
        res['is_clean'] = True
        res['branch_name'] = ''

    def mapped_cmake_version(last_release, is_release):
        """CMake expects MAJOR.MINOR.PATCH.TWEAK. The ancestral *last_release*
        is padded into the first three roles. If not *is_release*, the tweak role
        collects all postrelease states (prereleases and devel snapshots) into
        dummy 999 that at least gets them sorted correctly between releases and
        allows EXACT CMake version comparisons. Returns, for example, 1.1.0.0 for
        release 1.1, 1.3.4.0 for maintenance release 1.3.4, and 1.0.0.999 for
        prerelease 1.1a1 or snapshot 1.1.dev600

        """
        cm = last_release.split('.')
        cm += ['0'] * (4 - len(cm))
        if not is_release:
            cm[-1] = '999'
        cm = '.'.join(cm)
        return cm

    return {'__version__': project_version,
            '__version_long': project_version_long,
            '__version_is_clean': res['is_clean'],
            '__version_branch_name': res['branch_name'],
            '__version_last_release': backwardseries,
            '__version_cmake': mapped_cmake_version(backwardseries, project_release),
            '__version_release': project_release,
            '__version_prerelease': project_prerelease}


def write_new_metafile(versdata, outfile='metadata.out.py'):
    formatter_fn = """
def version_formatter(formatstring='{version}'):
    if formatstring == 'all':
        formatstring = '{version} {{{branch}}} {githash} {cmake} {clean} {release} {lastrel} <-- {versionlong}'

    release = 'release' if (__version_release == 'True') else ('prerelease' if (__version_prerelease == 'True') else '')

    ans = formatstring.format(version=__version__,
                              versionlong=__version_long,
                              githash=__version_long[len(__version__)+1:],
                              clean='' if __version_is_clean == 'True' else 'dirty',
                              branch=__version_branch_name,
                              lastrel=__version_last_release,
                              cmake=__version_cmake,
                              release=release)
    return ans
"""
    main_fn = """
if __name__ == '__main__':
    print(version_formatter(formatstring='all'))
"""
    with open(os.path.abspath(outfile), 'w') as handle:
        for k in sorted(versdata):
            handle.write("""{} = '{}'\n""".format(k, versdata[k]))
        handle.write(formatter_fn)
        handle.write(main_fn)


def version_formatter(versdata, formatstring="""{version}"""):
    """Return version information string with data from *versdata* when
    supplied with *formatstring* suitable for ``formatstring.format()``.
    Use plaintext and any placeholders among: version, versionlong, githash,
    branch, clean, release, lastrel, cmake. For example, '{branch}@{githash}'
    returns something like 'fix200@1234567'.

    """
    if formatstring == 'all':
        formatstring = '{version} {{{branch}}} {githash} {cmake} {clean} {release} {lastrel} <-- {versionlong}'

    release = 'release' if versdata['__version_release'] else ('prerelease' if versdata['__version_prerelease'] else '')

    ans = formatstring.format(version=versdata['__version__'],
                              versionlong=versdata['__version_long'],
                              githash=versdata['__version_long'][len(versdata['__version__']) + 1:],
                              clean='' if versdata['__version_is_clean'] else 'dirty',
                              branch=versdata['__version_branch_name'],
                              lastrel=versdata['__version_last_release'],
                              cmake=versdata['__version_cmake'],
                              release=release)
    return ans


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to extract Psi4 version from source. Use psi4.version_formatter(fmt_string) after build.')
    parser.add_argument('--metaout', default='metadata.out.py', help='file to which the computed version info written')
    parser.add_argument('--format', default='all', help='string like "{version} {githash}" to be filled in and returned')
    parser.add_argument('--formatonly', action='store_true', help='print only the format string, not the detection info')
    args = parser.parse_args()

    ans = reconcile_and_compute_version_output(quiet=args.formatonly)
    write_new_metafile(ans, args.metaout)
    ans2 = version_formatter(ans, formatstring=args.format)
    print(ans2)
