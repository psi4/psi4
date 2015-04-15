import sys
import subprocess

top_srcdir = './'
if len(sys.argv) == 2:
    top_srcdir = sys.argv[1] + '/'


def write_version(branch, mmp, ghash, status):
    if ghash:
        version_str = "#define GIT_VERSION \"{%s} %s %s\"\n" % \
                      (branch, ghash, status)
    else:
        version_str = "#undef GIT_VERSION"

    if mmp:
        mmp_str = "#define PSI_VERSION \"%s\"\n" % (mmp)
    else:
        mmp_str = "#undef PSI_VERSION"

    with open('gitversion.h', 'w') as handle:
        handle.write(version_str)
        handle.write(mmp_str)


# Use Git to find current branch name
#   Returns "refs/heads/BRANCHNAME"
#   Use [11:] to skip refs/heads/
try:
    command = "git symbolic-ref -q HEAD"
    process = subprocess.Popen(command.split(), stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE, cwd=top_srcdir)
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
                               stdout=subprocess.PIPE, cwd=top_srcdir)
    (out, err) = process.communicate()
    fields = str(out).rstrip().split('-')

    #         a68d223        # tags not pulled, clean git directory
    #         a68d223-dirty  # tags not pulled, changes to git-controlled files
    # 0.1-62-ga68d223        #  tags pulled, clean git directory
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

    if len(fields) == 2:
        mmp = '.'.join(fields)
    else:
        mmp = ''

    if process.returncode:
        status = ''
        ghash = ''
        mmp = ''
except:
    status = ''
    ghash = ''
    mmp = ''

write_version(branch, mmp, ghash, status)
