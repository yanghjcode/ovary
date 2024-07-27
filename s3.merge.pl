#!/usr/bin/perl


use strict;
use List::Util qw/min max/;

scalar @ARGV ==3 or die "\nUsage: perl $0 <meta.data.xls> <outperfix.pos> <outfile>\n\n";
my $meta = shift;
my $pos = shift;
my $out = shift;

my %meta;
my %name;
open IN,"<$meta" or die;
my %idx;
my $head=<IN>;chomp($head);my @head=split /\t/,$head;foreach my $num(0..$#head){$idx{$head[$num]}=$num;}
while(<IN>){
     chomp;
     my @ln=split /\t/,$_;
     my $id="$ln[$idx{x}]-$ln[$idx{y}]"; 
     $meta{$id}{"OO.genes"}=$ln[$idx{"OO.genes1"}];
     $meta{$id}{"GC.genes"}=$ln[$idx{"GC.genes1"}];
     $meta{$id}{"GC.C.genes"}=$ln[$idx{"GC.C.genes1"}];
     $meta{$id}{"GC.M.genes"}=$ln[$idx{"GC.M.genes1"}];
     $meta{$id}{"GC.3.genes"}=$ln[$idx{"GC.3.genes1"}];
     $meta{$id}{"TC.genes"}=$ln[$idx{"TC.genes1"}];
     $meta{$id}{"smc.genes"}=$ln[$idx{"smc.genes1"}];
     $meta{$id}{"endo.genes"}=$ln[$idx{"endo.genes1"}];
     $meta{$id}{"Macrophage.genes"}=$ln[$idx{"Macrophage.genes1"}];
     $name{"OO.genes"}=1;
     $name{"GC.genes"}=1;
     $name{"GC.C.genes"}=1;
     $name{"GC.M.genes"}=1;
     $name{"GC.3.genes"}=1;
     $name{"TC.genes"}=1;
     $name{"smc.genes"}=1;
     $name{"endo.genes"}=1;
     $name{"Macrophage.genes"}=1;
}
close IN;

my @names=sort keys %name;
my $names=join("\t",@names);
open IN,"<$pos" or die;
open OUT,">$out" or die;
my $head2=<IN>;chomp($head2);print OUT "$head2\t$names\n";
while(<IN>){
     chomp;
     my @ln=split;
     my $id="$ln[1]-$ln[2]";
     if(exists $meta{$id}){
           print OUT "$ln[0]\t$ln[1]\t$ln[2]\t$ln[3]\t$ln[4]";
           foreach my $n(@names){
                 print OUT "\t$meta{$id}{$n}";
           }  
           print OUT "\n";    
     }
}
close IN;
close OUT;


