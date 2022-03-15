__version__ = '1.6a1.dev12'
__version_branch_name = 'new_tamps'
__version_cmake = '1.5.0.999'
__version_is_clean = 'False'
__version_last_release = '1.5'
__version_long = '1.6a1.dev12+5249bc3'
__version_prerelease = 'False'
__version_release = 'False'

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

if __name__ == '__main__':
    print(version_formatter(formatstring='all'))
