#! /usr/bin/perl

# Compute Fragment Center of Mass Radial Distribution Function

# This script is 30x times slower than C version

use 5.008;
use strict;
use warnings;

use Getopt::Long;
use Math::Trig;

my $help = 0;
my $max_dist = 10.0;
my $n_bins = 100;

my $res = GetOptions('r=f' => \$max_dist, 'n=i' => \$n_bins, 'help' => \$help);
if ($help or !$res) {
	usage();
}

my @hist = (0) x $n_bins;
while (<>) {
	next unless /RESTART/;
	readline(STDIN);

	my @box = (0.0, 0.0, 0.0);
	my @points;

	while (<>) {
		last unless /fragment/;
		my @xyz = split / +/, readline(STDIN);
		push @points, [ $xyz[1], $xyz[2], $xyz[3] ];

		readline(STDIN);
		readline(STDIN);
		readline(STDIN);
	}

	while (<>) {
		last if /STATE AFTER/;

		if (/PERIODIC BOX SIZE/) {
			$_ =~ s/\s*PERIODIC BOX SIZE\s*//;
			@box = split / +/;
			last;
		}
	}

	for (my $i1 = 0; $i1 < scalar @points; $i1++) {
		for (my $i2 = $i1 + 1; $i2 < scalar @points; $i2++) {
			my $dist = get_dist($points[$i1], $points[$i2], \@box);

			if ($dist < $max_dist) {
				$hist[int($dist / $max_dist * $n_bins)] += 2;
			}
		}
	}
}

for (my $i = 1; $i < $n_bins; $i++) {
	my $r = $max_dist * $i / $n_bins;
	my $volume = 4.0 * pi * $r * $r * $max_dist / $n_bins;
	my $rdf = $hist[$i] / $volume;

	printf "%10.3f %10.3f\n", $r, $rdf;
}

sub get_dist {
	my ($p1, $p2, $box) = @_;

	my $dx = abs($p2->[0] - $p1->[0]);
	my $dy = abs($p2->[1] - $p1->[1]);
	my $dz = abs($p2->[2] - $p1->[2]);

	if ($box->[0] > 1.0e-8 && $box->[1] > 1.0e-8 && $box->[2] > 1.0e-8) {
		$dx -= $box->[0] * int(0.5 + $dx / $box->[0]);
		$dy -= $box->[1] * int(0.5 + $dy / $box->[1]);
		$dz -= $box->[2] * int(0.5 + $dz / $box->[2]);
	}

	return sqrt($dx * $dx + $dy * $dy + $dz * $dz);
}

sub usage {
	die <<EOF;
Usage:
	rdf.pl [-r <max_dist>] [-n <n_bins>]
	rdf.pl [-h | --help]

Defaults:
	max_dist 10.0
	  n_bins 100

Input is read from stdin
Output is written to stdout
EOF
}
