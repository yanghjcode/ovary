#!/usr/bin/perl


use strict;
use List::Util qw/min max/;

scalar @ARGV ==3 or die "\nUsage: perl $0 <RCTD.results.xls> <outperfix.pos> <outfile>\n\n";
my $rctd = shift;
my $pos = shift;
my $out = shift;

my %rctd;
my %name;
open IN,"<$rctd" or die;
my %idx;
<IN>;
while(<IN>){
     chomp;
     $_=~s/"//g;
     my @ln=split /,/,$_;
     my ($x,$y)=(split /[:_]/,$ln[0])[-2,-1];
     my $id="$x-$y"; 
     $rctd{$id}=$ln[2];
}
close IN;

open IN,"<$pos" or die;
open OUT,">$out" or die;
my $head2=<IN>;chomp($head2);print OUT "$head2\tType\n";
while(<IN>){
     chomp;
     my @ln=split;
     my $id="$ln[1]-$ln[2]";
     if(exists $rctd{$id}){
           print OUT "$ln[0]\t$ln[1]\t$ln[2]\t$ln[3]\t$ln[4]\t$rctd{$id}\n";
     }
}
close IN;
close OUT;


