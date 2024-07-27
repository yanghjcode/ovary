#!/usr/bin/perl

use strict;
use List::Util qw/min max/;

scalar @ARGV ==2 or die "\nUsage: perl $0 <gem> <outperfix>\n\n";
my $gem = shift;
my $out = shift;


my %sum;
my %pos;
if($gem=~/.gz/){open IN,"gzip -dc $gem|" or die;}else{open IN,"<$gem" or die;}
while(<IN>){
      chomp;
      next if(/^#/); 
      next if(/^geneID/);       
      my @ln=split;
      next if($ln[-1]==0);
      $sum{$ln[-1]}{mid}{$ln[0]}+=$ln[3];
      $sum{$ln[-1]}{pos}{"$ln[1]_$ln[2]"}=1;
      
}
close IN;
close OUT;

open POS,">$out.pos" or die;
open GEM,">$out.gem" or die;
print POS "CellID\tx_mean\ty_mean\tx\ty\n";
print GEM "geneID\tx\ty\tMIDCount\n";
foreach my $label(sort keys %sum){
      #pos
      my @x;
      my @y;
      my @pos;
      foreach my $pos(sort keys %{$sum{$label}{pos}}){            
            my ($x,$y)=(split /_/,$pos)[0,1];
            push @x,$x;
            push @y,$y;
      }
      my $avg_x = aver(\@x); $avg_x=int($avg_x);my $x=join(",",@x);
      my $avg_y =aver(\@y); $avg_y=int($avg_y); my $y=join(",",@y);
      foreach my $num(0..$#x){
             print POS "$label\t$avg_x\t$avg_y\t$x[$num]\t$y[$num]\n";
      }
      #mid
      foreach my $gene(sort keys %{$sum{$label}{mid}}){
             print GEM "$gene\t$avg_x\t$avg_y\t$sum{$label}{mid}{$gene}\n";       

      }
}
close OUT;

sub aver{
    my $arr = shift;
    my $s = 0;
    grep {$s += $_}@$arr;
    return $s/@$arr;
}
sub var {
    my $arr = shift;
    my $v = aver($arr);
    my $d = 0;
    grep {$d += ($_-$v)**2;}@$arr;
    return $d/(@$arr-1);
}

sub sd {
    my $arr = shift;
    my $d = 0;
    my $var=var($arr);
    return sqrt($var);
}

