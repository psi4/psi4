#! /usr/bin/perl

# This script extracts trajectory data from EFPMD output

use 5.008;
use strict;
use warnings;

die "usage: trajectory.pl <input>\n" if scalar(@ARGV) != 1;

open(FH, "<", $ARGV[0]) || die "$!";

while (<FH>) {
	next unless /GEOMETRY/;
	<FH>;

	my @lines;

	while (<FH>) {
		last if /^$/;
		push @lines, $_;
	}

	print scalar(@lines), "\n";
	print "xyz", "\n";

	foreach (@lines) {
		$_ =~ s/A[0-9]+([a-zA-Z]+)[0-9]*/$1/;
		print $_;
	}
}
