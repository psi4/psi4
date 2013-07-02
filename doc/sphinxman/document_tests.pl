#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(remove_tree);

# This script looks for comments in each test case to identify
# the type of calculation.  The results of this parsing is put into TeX files, which are
# inlined into the manual.
#
# Then, we look at the test cases, using the specially formatted comments to add a listing
# of samples in the users' manual; the tests are copied (sans validation stuff) to the appropriate
# location in the samples folder.


my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }


#
# Raid the test cases looking for tags
#

my $SamplesFolder = $DriverPath . "../../samples";
opendir(SAMPLES, $SamplesFolder) or die "I can't read $SamplesFolder\n";
foreach my $File(readdir SAMPLES){
    next unless $File =~ /\w+/; # Make sure we don't nuke .. or . !
    next if $File =~ /example_psi4rc_file/; # Keep the example psi4rc file
    remove_tree("$SamplesFolder/$File");
}

my $TestsFolder = $DriverPath . "../../tests";
opendir(TESTS, $TestsFolder) or die "I can't read $TestsFolder\n";
my $TexSummary = "tests_descriptions.tex";
my $RstSummary = "source/autodoc_testsuite.rst";
# Create a plain-text summary in the samples directory
my $Summary = $SamplesFolder."/SUMMARY";
open(SUMMARY,">$Summary") or die "I can't write to $Summary\n";
# Make a LaTeX version for the manual, too
open(TEXSUMMARY,">$TexSummary") or die "I can't write to $TexSummary\n";
open(RSTSUMMARY,">$RstSummary") or die "I can't write to $RstSummary\n";
print "Auto-documenting samples directory inputs\n";
print RSTSUMMARY "\n=============================================   ============\n";
print RSTSUMMARY   "Input File                                      Description \n";
print RSTSUMMARY   "=============================================   ============\n";
foreach my $Dir(readdir TESTS){
    my $Input = $TestsFolder."/".$Dir."/input.dat";
    # Look for an input file in each subdirectory, or move on
    open(INPUT, "<$Input") or next;
    #
    # Remember, we want to copy the test suite over to psi4/samples, omitting the test functions
    #
    my $SamplesDirectory = $SamplesFolder."/".$Dir;
    unless(-d $SamplesDirectory){
        # This directory doesn't exist in psi4/samples, make it now
        mkdir $SamplesDirectory or die "\nI can't create $SamplesDirectory\n";
    }
    my $SampleInput = $SamplesDirectory."/input.dat";
    open(SAMPLE, ">$SampleInput") or die "\nI can't write to $SampleInput\n";
    my $Description;
    my $TestOnlyNoSample = 0;
    while(<INPUT>){
        # If this line isn't associated with testing, put it in the sample
        print SAMPLE unless /\#TEST/;
        # Now we only want to grab the comments.  Move on if this is not a comment.
        next unless s/\#\!//;
        $Description .= $_;
        chomp $Description;
        if ($Description =~ /\!nosample/) {
            $TestOnlyNoSample = 1;
        }
    }
    close INPUT;
    close SAMPLE;

    if ($TestOnlyNoSample) {
        unlink $SampleInput or warn "Could not unlink $SampleInput: $!";
        rmdir $SamplesDirectory or die "\nI can't remove $SamplesDirectory\n";
    }

    # Process the comment that we grabbed from the input.
    if($Description){
        # make directory name tex safe
        my $Dir_tex = $Dir;
        $Dir_tex =~ s/_/\\_/g;
        my $Description_tex = $Description;
        $Description_tex =~ s/_/\\_/g;
        $Description_tex =~ s/@@/_/g;
        my $Description_rst = $Description_tex;
        $Description_rst =~ s/ \$/ :math:`/g;
        $Description_rst =~ s/\(\$/(\\ :math:`/g;
        $Description_rst =~ s/\$ /` /g;
        $Description_rst =~ s/\$\./`./g;
        $Description_rst =~ s/\$,/`,/g;
        $Description_rst =~ s/\$\)/`\\ \)/g;
        $Description_rst =~ s/\\_/_/g;
        if ($TestOnlyNoSample == 0) {
            print TEXSUMMARY "\\begin{tabular*}{\\textwidth}[tb]{p{0.2\\textwidth}p{0.8\\textwidth}}\n";
            print TEXSUMMARY "{\\bf $Dir_tex} & $Description_tex \\\\\n\\\\\n";
            print TEXSUMMARY "\\end{tabular*}\n";
            my $srcfilename = ":srcsample:`" . $Dir_tex . "`";
            printf RSTSUMMARY "%-45s  %s\n", $srcfilename, $Description_rst;
            printf SUMMARY "%-12s %s\n\n\n", $Dir.":", $Description;
        }
    }else{
        warn "Warning!!! Undocumented input: $Input\n";
    }
}
print RSTSUMMARY "=============================================   ============\n\n";
close TEXSUMMARY ;
close RSTSUMMARY;
close SUMMARY;
closedir TESTS;

