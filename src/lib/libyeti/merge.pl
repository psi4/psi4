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
    print "Removing $File...\n";
    unlink($File) or die "\nFailed!\n";
}
    
# Copy over the Yeti files
my $YetiSrc = $Yeti."/src/";
opendir(YETIDIR, "$YetiSrc") or die "I can't read $YetiSrc\n";
foreach my $YetiFile (readdir YETIDIR){
    next unless $YetiFile =~ /\.(cc|hpp|h)$/;
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
        # src/smartptr/src/... -> libsmartptr/...
        s/src\/smartptr\/src/libsmartptr/g;
        # src/... -> ...
        s/src\//libyeti\//g;
        $Contents .= $_;
    }
    close FILE;
    open(FILE, ">$File") or die "I can't write to $File\n";
    print FILE $Contents;
    close FILE;
}
