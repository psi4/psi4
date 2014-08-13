# GA_F77_LD_OVERRIDE
# ------------------
# Replaced $(F77LINK) in our Makefile.am's with $(LINK) if needed.
#
# There is a ton of garbage here.  We must replicate what automake will do in
# the case where we enable F77 code (and therefore require the F77LINK).  We
# use autoconf quadrigraphs to avoid problems with $(foo) and $@.  This is yet
# another hack to disable F77 since automake selects the linker based on the
# static list of source files.
AC_DEFUN([GA_F77_LD_OVERRIDE], [
m4_pattern_allow([AM_V_lt])
AS_IF([test "x$enable_f77" = xyes],
    [F77LINK='@S|@@{:@LIBTOOL@:}@ @S|@@{:@AM_V_lt@:}@ --tag=F77 @S|@@{:@AM_LIBTOOLFLAGS@:}@ @S|@@{:@LIBTOOLFLAGS@:}@ --mode=link @S|@@{:@F77LD@:}@ @S|@@{:@AM_FFLAGS@:}@ @S|@@{:@FFLAGS@:}@ @S|@@{:@AM_LDFLAGS@:}@ @S|@@{:@LDFLAGS@:}@ -o @S|@@'
     am__v_F77LD_0='@echo "  F77LD " @S|@@;'],
    [F77LINK='@S|@@{:@LINK@:}@'
     am__v_F77LD_0='@S|@@{:@am__v_CCLD_0@:}@'])
AC_SUBST([F77LINK])
AC_SUBST([am__v_F77LD_0])
])dnl
