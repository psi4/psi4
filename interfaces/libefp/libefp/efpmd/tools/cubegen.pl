#! /usr/bin/perl

# This script generates a rectangular box of fragments for EFPMD input

use 5.008;
use strict;
use warnings;

use Math::Trig;

if (scalar(@ARGV) != 6) {
	die <<EOF;
usage: cubegen.pl <n1:n2:...> <r1:r2:...> <space> <nx> <ny> <nz>

  <n1:n2:...>  a colon separated list of fragment names
  <r1:r2:...>  a colon separated list of fragment ratios
      <space>  distance between the fragments
         <nx>  number of fragments in x direction
         <ny>  number of fragments in y direction
         <nz>  number of fragments in z direction

example: cubegen.pl c2h5oh_l:h2o_l 40:60 4.0 20 20 20 > vodka.in
EOF
}

my @name = split ':', $ARGV[0];
my @ratio = split ':', $ARGV[1];
my $space = $ARGV[2];
my $nx = $ARGV[3];
my $ny = $ARGV[4];
my $nz = $ARGV[5];

die "ratios do not match names" unless scalar(@name) == scalar(@ratio);
die "positive number expected" unless $space > 0.0 && $nx > 0 && $ny > 0 && $nz > 0;

foreach my $i (1 .. scalar(@ratio) - 1) {
	$ratio[$i] += $ratio[$i - 1];
}

foreach my $x (0 .. $nx - 1) {
	foreach my $y (0 .. $ny - 1) {
		foreach my $z (0 .. $nz - 1) {
			print_fragment($x, $y, $z);
		}
	}
}

sub select_fragment {
	my $n = int(rand($ratio[-1]));

	foreach my $i (0 .. scalar(@ratio) - 1) {
		if ($ratio[$i] > $n) {
			return $i;
		}
	}

	die;
}

sub print_fragment {
	my ($x, $y, $z) = @_;
	my $idx = select_fragment();
	my @xyzabc = ($space * $x,
		      $space * $y,
		      $space * $z,
		      2.0 * pi * rand,
		      pi * rand,
		      2.0 * pi * rand);

	print "\n", "fragment ", $name[$idx], "\n";

	foreach (@xyzabc) {
		printf " %8.3f", $_;
	}

	print "\n";
}
