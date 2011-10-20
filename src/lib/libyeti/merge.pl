#!/usr/bin/perl

# This script copies the standalone Yeti and smartptr codes into
# the Psi distribution
#
# Andy Simmonett Jan '11

use strict;
use warnings;
use File::Copy;

die "Usage:\n\t$0 yeti_topdir\n\n\tN.B., this should be run from PSI4/src/lib/libyeti"
     unless @ARGV == 1;

my $Yeti = shift;

# Remove existing Yeti files
opendir(DIR, ".") or die "\nI can't open this directory\n";
foreach my $File (readdir DIR){
    next unless $File =~ /\.(cc|hpp|h)$/;
    next if $File =~ /psi4_interface/; # Don't remove the interface!
    print "Removing $File...\n";
    unlink($File) or die "\nFailed!\n";
}
    
# Copy over the Yeti files
my $YetiSrc = $Yeti."/src/";
opendir(YETIDIR, "$YetiSrc") or die "I can't read $YetiSrc\n";
foreach my $YetiFile (readdir YETIDIR){
    next unless $YetiFile =~ /\.(cc|hpp|h)$/;
    next if $YetiFile =~ /mpqc/; # Don't add the MPQC interface files!
    print "Copying $YetiSrc$YetiFile to $YetiFile...\n";
    copy($YetiSrc.$YetiFile, $YetiFile) or die "Failed!\n";
}

# Remove existing SmartPtr files
opendir(DIR, "../libsmartptr") or die "\nI can't open ../libsmartptr\n";
foreach my $File (readdir DIR){
    next unless $File =~ /\.(cc|hpp|h)$/;
    print "Removing ../libsmartptr/$File...\n";
    unlink("../libsmartptr/".$File) or die "\nFailed!\n";
}

# Copy over the smartptr files
my $SmartPtrSrc = $Yeti."/src/smartptr/src/";
opendir(SMARTPTRDIR, "$SmartPtrSrc") or die "I can't read $SmartPtrSrc\n";
foreach my $SmartPtrFile (readdir SMARTPTRDIR){
    next unless $SmartPtrFile =~ /\.(cc|hpp|h)$/;
    print "Copying $SmartPtrSrc$SmartPtrFile to $SmartPtrFile...\n";
    copy($SmartPtrSrc.$SmartPtrFile, "../libsmartPtr/".$SmartPtrFile) or die "Failed!\n";
}

# Tweak the paths in the headers
opendir(DIR, ".") or die "I can't read this directory\n";
foreach my $File (readdir DIR){
    next unless $File =~ /\.(h|cc)$/;
    print "Processing $File...\n";
    open(FILE, "<$File") or die "I can't read $File\n";
    my $Contents;
    while(<FILE>){
        # config.h
        next if /config\.h/;
        # src/smartptr/src/... -> libsmartptr/...
        s/src\/smartptr\/src/libsmartptr/g;
        # src/... -> ...
        s/src\//libyeti\//g;
        $Contents .= $_;
    }
    # Put the correct psi MPI definitions in there
    $Contents =~ s/\#if HAVE_MPI_H/\#include <psiconfig.h>\n\#if HAVE_MPI/g;
    # Strip the BLAS/LAPACK declarations
    $Contents =~ s/(extern.*\/\/EndExternC)/#include "blas.h"\n$1/smg;
    # Change all the YETI_DGEMM like calls to C_DGEMM.  This might need to change in the future
    $Contents =~ s/YETI_D/F_D/g;
    close FILE;
    open(FILE, ">$File") or die "I can't write to $File\n";
    print FILE $Contents;
    close FILE;
}

open(BLAS,">blas.h");
print BLAS
"#if FC_SYMBOL==2
#define F_DGBMV dgbmv_
#define F_DGEMM dgemm_
#define F_DGEMV dgemv_
#define F_DGER dger_
#define F_DSBMV dsbmv_
#define F_DSPMV dspmv_
#define F_DSPR dspr_
#define F_DSPR2 dspr2_
#define F_DSYMM dsymm_
#define F_DSYMV dsymv_
#define F_DSYR dsyr_
#define F_DSYR2 dsyr2_
#define F_DSYR2K dsyr2k_
#define F_DSYRK dsyrk_
#define F_DTBMV dtbmv_
#define F_DTBSV dtbsv_
#define F_DTPMV dtpmv_
#define F_DTPSV dtpsv_
#define F_DTRMM dtrmm_
#define F_DTRMV dtrmv_
#define F_DTRSM dtrsm_
#define F_DTRSV dtrsv_
#elif FC_SYMBOL==1
#define F_DGBMV dgbmv
#define F_DGEMM dgemm
#define F_DGEMV dgemv
#define F_DGER dger
#define F_DSBMV dsbmv
#define F_DSPMV dspmv
#define F_DSPR dspr
#define F_DSPR2 dspr2
#define F_DSYMM dsymm
#define F_DSYMV dsymv
#define F_DSYR dsyr
#define F_DSYR2 dsyr2
#define F_DSYR2K dsyr2k
#define F_DSYRK dsyrk
#define F_DTBMV dtbmv
#define F_DTBSV dtbsv
#define F_DTPMV dtpmv
#define F_DTPSV dtpsv
#define F_DTRMM dtrmm
#define F_DTRMV dtrmv
#define F_DTRSM dtrsm
#define F_DTRSV dtrsv
#elif FC_SYMBOL==3
#define F_DGBMV DGBMV
#define F_DGEMM DGEMM
#define F_DGEMV DGEMV
#define F_DGER DGER
#define F_DSBMV DSBMV
#define F_DSPMV DSPMV
#define F_DSPR DSPR
#define F_DSPR2 DSPR2
#define F_DSYMM DSYMM
#define F_DSYMV DSYMV
#define F_DSYR DSYR
#define F_DSYR2 DSYR2
#define F_DSYR2K DSYR2K
#define F_DSYRK DSYRK
#define F_DTBMV DTBMV
#define F_DTBSV DTBSV
#define F_DTPMV DTPMV
#define F_DTPSV DTPSV
#define F_DTRMM DTRMM
#define F_DTRMV DTRMV
#define F_DTRSM DTRSM
#define F_DTRSV DTRSV
#elif FC_SYMBOL==4
#define F_DGBMV DGBMV_
#define F_DGEMM DGEMM_
#define F_DGEMV DGEMV_
#define F_DGER DGER_
#define F_DSBMV DSBMV_
#define F_DSPMV DSPMV_
#define F_DSPR DSPR_
#define F_DSPR2 DSPR2_
#define F_DSYMM DSYMM_
#define F_DSYMV DSYMV_
#define F_DSYR DSYR_
#define F_DSYR2 DSYR2_
#define F_DSYR2K DSYR2K_
#define F_DSYRK DSYRK_
#define F_DTBMV DTBMV_
#define F_DTBSV DTBSV_
#define F_DTPMV DTPMV_
#define F_DTPSV DTPSV_
#define F_DTRMM DTRMM_
#define F_DTRMV DTRMV_
#define F_DTRSM DTRSM_
#define F_DTRSV DTRSV_
#endif
#if FC_SYMBOL==2
#define F_DGEEV dgeev_
#define F_DGESV dgesv_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DPOTRF dpotrf_
#define F_DPOTRI dpotri_
#define F_DPOTRS dpotrs_
#define F_DGESVD dgesvd_
#define F_DSYEV dsyev_
#define F_DSPSVX dspsvx_
#elif FC_SYMBOL==1
#define F_DGEEV dgeev
#define F_DGESV dgesv
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DPOTRF dpotrf
#define F_DPOTRI dpotri
#define F_DPOTRS dpotrs
#define F_DGESVD dgesvd
#define F_DSYEV dsyev
#define F_DSPSVX dspsvx
#elif FC_SYMBOL==3
#define F_DGEEV DGEEV
#define F_DGESV DGESV
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DPOTRF DPOTRF
#define F_DPOTRI DPOTRI
#define F_DPOTRS DPOTRS
#define F_DGESVD DGESVD
#define F_DSYEV DSYEV
#define F_DSPSVX DSPSVX
#elif FC_SYMBOL==4
#define F_DGEEV DGEEV_
#define F_DGESV DGESV_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DPOTRF DPOTRF_
#define F_DPOTRI DPOTRI_
#define F_DPOTRS DPOTRS_
#define F_DGESVD DGESVD_
#define F_DSYEV DSYEV_
#define F_DSPSVX DSPSVX_
#endif
";
close BLAS;
