# Author: Shujia Huang
# Modify Date: 2017-03-05

#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

main();
exit;

sub main {
  &usage if (@ARGV < 1);
  my $command = shift(@ARGV);
  my %func = (subsam=>\&subsam);
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}};
}

sub subsam {
  die(qq/Usage: extract.pl subsam <in.vcf> [samples]\n/) if (@ARGV == 0);
  my ($fh, %h);
  my $fn = shift(@ARGV);
  my @col;
  open($fh, ($fn =~ /\.gz$/)? "gzip -dc $fn |" : $fn) || die;
  $h{$_} = 1 for (@ARGV);
  while (<$fh>) {
    if (/^##/) {
      print;
    } elsif (/^#/) {
      my @t = split;
      my @s = @t[0..8]; # all fixed fields + FORMAT
      for (9 .. $#t) {
        if ($h{$t[$_]}) {
          push(@s, $t[$_]);
          push(@col, $_);
        }
      }
      pop(@s) if (@s == 9); # no sample selected; remove the FORMAT field
      print join("\t", @s), "\n";
    } else {
      my @t = split;
      if (@col == 0) {
        print join("\t", @t[0..7]), "\n";
      } else {
        my ($ac,$af) = calcu_af(map {$t[$_]} @col);
        next if $af eq 0 or $af eq "NA" or $ac eq 0 or $ac eq "NA";
        $t[7] =~ s/AC=[^;]+/AC=$ac/;
        $t[7] =~ s/;AF=[^;]+/;AF=$af/;
        print join("\t", @t[0..8], map {$t[$_]} @col), "\n";
      }
    }
  }
  close($fh);
}

sub calcu_af {
  my @gnt = @_;
  my ($a, $f) = (0, 0);
  for my $gnt (@gnt) {
    $gnt = (split /:/, $gnt)[0];
    next if $gnt =~ /\./;
    my @g = split /[\/\|]/, $gnt;
    @g = ($g[0], 0) if (@g < 2);

    $g[0] = 1 if $g[0] > 1;
    $g[1] = 1 if $g[1] > 1;

    $a += ($g[0] + $g[1]);
    $f += 2;
  }
  if($f == 0){
    return ("NA","NA");
    next;
  }
  return ($a, sprintf ("%.4f", $a/$f));
} 

sub usage {
  die(qq/
Usage: perl $0 <command> [<arguments>]\n
Command: subsam       get a subset of samples
\n/);
}
