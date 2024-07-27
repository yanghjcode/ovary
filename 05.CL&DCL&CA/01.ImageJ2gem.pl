#!/usr/bin/perl

use strict;
use File::Basename qw(basename dirname);
use List::Util qw/max min/;
use FindBin '$Bin';
use Getopt::Long;

=head1
       Description
           generate gem from ImageJ
       usage
           perl -image .txt -gem gem -out outgem


=cut


my($image,$gem,$out);
GetOptions("image:s"=>\$image,"gem:s"=>\$gem,"out:s"=>\$out);
die `pod2text $0` unless($image or $gem or $out);

my %image;
open IN,"<$image" or die;
while(<IN>){
     chomp;
          my @ln=split;
          next if($ln[2]<20);
          my $x=int($ln[0]/20)*20;
          my $y=int($ln[1]/20)*20;
               $image{"$x\t$y"}=1;
          }
close IN;

my %gem;
if($gem=~/gz/){
    open IN,"gzip -dc $gem|" or die;
}else{
    open IN,"<$gem" or die;
}
<IN>;
while(<IN>){
      chomp;
      my @ln=split;    
      my $x=int($ln[1]/20)*20;
      my $y=int($ln[2]/20)*20;       
      $gem{$ln[0]}{"$x\t$y"}+=$ln[3];
}
close IN;

open OUT,"|gzip >$out.gem.gz" or die;
print OUT "geneID\tx\ty\tMIDCount\n";
foreach my $gene(sort keys %gem){
    foreach my $id(sort keys %image){
         if(exists $gem{$gene}{$id}){
               print OUT "$gene\t$id\t$gem{$gene}{$id}\n";
         }
    }
}
close OUT;

sub sum{
    my $arr = shift;
    my $s = 0;
    grep {$s += $_}@$arr;
    return $s;
}


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











