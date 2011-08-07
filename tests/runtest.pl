use strict;
use warnings;

die "Usage\n\t\t$0 input_file logfile\n" unless @ARGV == 2;

my $psi = "../../bin/psi4";
my $input = shift;
my $logfile = shift;

open(LOGFILE,">$logfile") or die "I can't write to $logfile\n";

sub backtick(@)
{
    my $pid = open(KID, '-|');
    die "fork: $!" unless defined($pid);
    if ($pid) {
        my $output;
        while (<KID>) {
            print STDOUT;
            print LOGFILE;
        }
        close(KID);
        return $output;
    } else {
        exec @_;
    }
}

my @cmd = ($psi, $input);
backtick(@cmd);

my $return = $? ? 1 : 0;
close LOGFILE;

exit $return;


