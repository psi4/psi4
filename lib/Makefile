
.PHONY:	default all install depend clean dclean targetclean

#subdirs = ruby

default: all

all:

install_inc:
depend:
clean:
dclean:
targetclean:

top_srcdir = /Users/masaaki/psi4
prefix = /usr/local/psi4
srcdir = /Users/masaaki/psi4/lib
datadir = ${datarootdir}
datarootdir = ${prefix}/share
pkgdatadir = $(datadir)/psi
VPATH = /Users/masaaki/psi4/lib
top_objdir = ..
MKDIRS = $(top_srcdir)/bin/mkdirs.sh
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644

install: psi.dat
	$(MKDIRS) $(DESTDIR)$(pkgdatadir)
	$(INSTALL_DATA) $(srcdir)/psi.dat $(DESTDIR)$(pkgdatadir)
	for dir in $(subdirs); \
          do \
            (cd $${dir}; echo Making in $${dir}; $(MAKE) install) || exit 1; \
          done

install_man: macro.psi
	$(MKDIRS) $(DESTDIR)$(pkgdatadir)
	$(INSTALL_DATA) $< $(DESTDIR)$(pkgdatadir)

$(top_srcdir)/configure: $(top_srcdir)/configure.ac $(top_srcdir)/aclocal.m4
	cd $(top_srcdir) && autoconf

$(top_objdir)/config.status: $(top_srcdir)/configure
	cd $(top_objdir) && ./config.status --recheck

Makefile: $(srcdir)/Makefile.in $(top_objdir)/config.status
	cd $(top_objdir) && CONFIG_FILES=lib/Makefile ./config.status
