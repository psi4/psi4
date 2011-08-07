use strict;
use warnings;

die "Usage\n\t\t$0 input_file target\n" unless @ARGV == 2;

my $psi = "../../bin/psi4";
my $input = shift;
my $target = shift;

sub backtick(@)
{
    my $pid = open(KID, '-|');
    die "fork: $!" unless defined($pid);
    if ($pid) {
        my $output;
        while (<KID>) {
            print STDOUT $_;
            $output .= $_; # could be improved...
        }
        close(KID);
        return $output;
    } else {
        exec @_;
    }
}

my @cmd = ($psi, $input);
my $output = backtick(@cmd);
my $return = $? ? 1 : 0;
open(TARGET,">$target") or die "I can't write to $target\n";
print TARGET $output;
close TARGET;

exit $return;


