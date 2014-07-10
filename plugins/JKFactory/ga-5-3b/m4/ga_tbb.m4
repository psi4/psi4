# GA_TBB()
# --------
# Locate Intel Threading Building Blocks.
# 
# They may have been installed to a typical root/include and root/lib location,
# or the user may have downloaded a source or binary tarball.
#
# The source tarballs contain a Makefile to be used with gmake while the binary
# tarball contains a number of libraries.  In either case, the TBB README
# recommends setting the $TBBROOT environment variable to point to the base
# directory of either the source or binary untarred directories. This isn't
# sufficient for locating libtbb.so since the source release builds both
# 'debug' and 'release' versions of the library into a 'work_dir' location
# based on the user's system information (uname, etc) and the binary release
# contains multiple versions of the library based on many different systems and
# CPU architectures. Ugh.
#
# If $TBBROOT is set (or $tbb_root), we first search for the
# $TBBROOT/build/common.inc Makefile fragment since this will help us find the
# appropriate libtbb.so.
#
# If $TBBROOT is not set, we try a standard lookup of the package which will
# most likely fail unless the user helps.
#
AC_DEFUN([GA_TBB], [
tbb_ok=no
# has the user set TBBROOT or tbb_root?
AS_IF([test "x$TBBROOT" != x || test "x$tbb_root" != x], [
rm -f the_makefile
rm -f result.txt
rm -f out.txt
cat >the_makefile <<"EOF"
ifdef TBBROOT
tbb_root=$(TBBROOT)
endif
include $(tbb_root)/build/common.inc
result.txt:
	@echo "$(work_dir)" > result.txt
EOF
    AS_IF([gmake -f the_makefile &> out.txt],
        [tbb_work_dir=`cat result.txt`; tbb_ok=yes],
        [cat out.txt >&AS_MESSAGE_LOG_FD])
    AS_IF([test "x$tbb_work_dir" != x],
        [TBB_LDFLAGS="-L${tbb_work_dir}_release"])
    AS_IF([test "x$TBBROOT" != x],
        [TBB_CPPFLAGS="-I$TBBROOT/include"],
        [test "x$tbb_root" != x],
        [TBB_CPPFLAGS="-I$tbb_root/include"])
    TBB_LIBS="-ltbb"
    rm -f the_makefile
    rm -f result.txt
    rm -f out.txt
])
AS_IF([test "x$tbb_ok" = xno], [
GA_CHECK_PACKAGE([tbb], [tbb/tbb.h],
    [tbb], [TBB_runtime_interface_version],
    [], [tbb_ok=yes],
    [AC_MSG_ERROR([could not locate Intel Thread Building Blocks])])
])
AC_SUBST([TBB_CPPFLAGS])
AC_SUBST([TBB_LDFLAGS])
AC_SUBST([TBB_LIBS])
])
