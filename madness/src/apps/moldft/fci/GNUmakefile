#
# $Id: GNUmakefile,v 1.2 1999/07/28 00:23:43 d3e129 Exp $
#

   ARFLAGS = urs
    RANLIB = ranlib

# SUN SOLARIS
#       CC = cc
#   CFLAGS = -O
#       FC = f77
#   FFLAGS = -fast -O4 -stackvar -fsimple=2 -depend -dalign
#     LIBS = 
# INCLUDES = -I.
#  DEFINES = -DSUN
# CPPFLAGS = $(INCLUDES) $(DEFINES)
#    BLAS := dgemm.o idamax.o crayio.o dscal.o daxpy.o dcopy.o ddot.o 
#  LDFLAGS =

# LINUX with g77
       CC = gcc
   CFLAGS = -g -O
       FC = gfortran
   FFLAGS = -O2 -malign-double -ffast-math -fno-second-underscore
## -funroll-loops ... introduces bugs in mfind_weight
##g -Wimplicit -Wunused 
     LIBS = -L.
 INCLUDES = -I.
  DEFINES = -DLINUX
 CPPFLAGS = $(INCLUDES) $(DEFINES)
 BLAS := linux_cputime.o dgemm.o idamax.o crayio.o dscal.o daxpy.o dcopy.o ddot.o
  LDFLAGS =

# IBM AIX
#EXPLICITF = TRUE
#       FC = xlf
#   FFLAGS = -O3 -qnosave -qalign=4k 
#     LIBS = -lblas
# INCLUDES = -I.
#  DEFINES = -DIBM -DAIX
# CPPFLAGS = $(INCLUDES) $(DEFINES)
# FCONVERT = $(CPP) $(CPPFLAGS)  $< | sed '/^\#/D'  > $*.f
#     BLAS = 


############################################
# Below here should be machine independent #
############################################


 OBJ :=  fci_mfndwt.o  fci_ofndoc.o  fci_num.o \
         fci_ifill.o   fci_mgraph.o  fci_ofndwt.o  fci_setup.o \
         fci_mfndoc.o  fci_ograph.o  \
         fci_err.o     fci_tinit.o   fci_cntone.o  fci_spnad.o  \
         fci_prspnad.o \
         fci_dfill.o   fci_spnfrc.o  screen.o      renorm.o     \
         fci_onelcc.o  fci_sigma.o  fci_diags.o \
         fci_load.o \
         fci_srt1cc.o  yacobi.o      fci_hmod.o   fci_prntri.o \
         fci_make_eaj.o fci_prteop.o fci_hsort.o   fci_detwalk.o \
         fci_owalk_info.o fci_couple1.o fci_1pdm.o fci_output.o \
         fci_2pdm.o    fci_couple2.o  fci_3pdm.o fci_couple3.o \
         fci_hami.o    fci_onepdm.o  fci_intchk.o \
         fci_hv.o      fci_axb.o fci_solve.o fci_davids.o rsg.o $(BLAS)

fci:	fci.o libfci.a
	$(FC) $(LDFLAGS) $(FFLAGS) -o $@ $^ $(LIBS)

timer:	timer.o libfci.a
	$(FC) $(LDFLAGS) $(FFLAGS) -o $@ timer.o libblas.a libfci.a $(LIBS)

libfci.a:	$(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

clean:
	-$(RM) -f core fci *.o libfci.a

realclean:	clean
	-$(RM) -f *~ \#*\# *mon.out *.prof a.out TAGS

$(OBJ):	fciinfo.fh fciparam.fh

fci_tinit.o:	fcitime.fh

ifdef EXPLICITF
#
# Needed on when FC does not preprocess .F files with CPP to get .f files
#
.SUFFIXES:	
.SUFFIXES:	.o .s .F .f .c

 FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
		$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1

.F.o:	
	@echo Converting $*.F '->' $*.f
	@$(FCONVERT)
	$(FC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.F.f:
	@echo Converting $*.F '->' $*.f
	@$(FCONVERT)
endif
# 

