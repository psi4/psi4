
.PHONY:	default all install depend clean dclean targetclean

default: all

all:

install_man:
depend:
clean:
dclean:
targetclean:

top_srcdir = /Users/masaaki/psi4
prefix = /usr/local/psi4
srcdir = /Users/masaaki/psi4/include
incdir = ${prefix}/include
VPATH = /Users/masaaki/psi4/include
top_objdir = ..
MKDIRS = $(top_srcdir)/bin/mkdirs.sh
INSTALL = /usr/bin/install -c
INSTALL_INCLUDE = $(top_srcdir)/bin/installinc.sh

HEADERS = psiconfig.h psitypes.h ccfiles.h masses.h physconst.h \
psifiles.h symmetry.h rgb.h vdw_radii.h chkpt_params.h psi4-dec.h psi4-def.h \
element_to_Z.h Z_to_element.h molecular_system.h fragment.h

install: install_inc

install_inc: $(HEADERS)
	$(MKDIRS) $(DESTDIR)$(incdir)
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)

$(top_srcdir)/configure: $(top_srcdir)/configure.ac $(top_srcdir)/aclocal.m4
	cd $(top_srcdir) && autoconf

$(top_objdir)/config.status: $(top_srcdir)/configure
	cd $(top_objdir) && ./config.status --recheck

Makefile: $(srcdir)/Makefile.in $(top_objdir)/config.status
	cd $(top_objdir) && CONFIG_FILES=include/Makefile ./config.status
