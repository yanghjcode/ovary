#!/usr/bin/perl
use strict;
use File::Basename qw(basename dirname);
use List::Util qw/max min/;
use List::MoreUtils ':all';
use FindBin '$Bin';
use Getopt::Long;

=head1
       Description
           generate oo gem
       usage
           perl  -rctd RCTD.results.csv -neighbor num  -outfile 


=cut


my($rctd,$neighbor_num,$out);
GetOptions("rctd:s"=>\$rctd,"neighbor:s"=>\$neighbor_num, "outfile:s"=>\$out);
die `pod2text $0` unless($rctd or $neighbor_num or $out);

######################################################

my $bin=20;

my %score;

open IN,"$rctd" or die;
my %index;
my $head=<IN>;chomp($head);$head=~s/"//g; my @head=split /,/,$head;
foreach my $num(0..$#head){
     $index{$head[$num]}=$num; print "$num\t$head[$num]\t$index{$head[$num]}\n";
}
while(<IN>){
     chomp;
     $_=~s/"//g;
     my @ln=split /,/,$_;
     next if( $ln[$index{"first_type"}] ne "7_oo" && $ln[$index{"first_type"}] ne "3_gran");

     $score{$ln[0]}{"OO.genes1"}=$ln[$index{"OO.genes1"}];
     my $idoo="OO.genes1";
     $score{$ln[0]}{"GC.genes1"}=$ln[$index{"GC.genes1"}];
     my $idgc="GC.genes1";
     $score{$ln[0]}{"TC.genes1"}=$ln[$index{"TC.genes1"}];
     my $idtc="TC.genes1";
     $score{$ln[0]}{"celltype"}=$ln[$index{"first_type"}]; 
}
close IN;


my %core;
foreach my $spot(sort keys %score){           
      my $neighbor=0;
      my ($core_x,$core_y)=(split /[:_]/,$spot)[-2,-1];     
      my $sample=(split /:/,$spot)[0];
      my @x=($core_x-20,$core_x,$core_x+20);
      my @y=($core_y-20,$core_y,$core_y+20);
      foreach my $x(@x){
            foreach my $y(@y) {
                 next if ($x==$core_x && $y==$core_y);
                        my $id="$sample:$x\_$y";
                        if(exists  $score{$id}){
                             $neighbor++;
                        }
           }
      }
      $score{$spot}{"neighbor"}=$neighbor;
      if($neighbor<=$neighbor_num){
             delete $score{$spot};
      } 
}

my $num=0;
my %pos;
foreach my $c(sort keys %score){
      my ($x,$y)=(split /[:_]/,$c)[-2,-1];
      $pos{"$x\_$y"}=1;
}



my %score_raw=%score;
my $merge_num=1;
      LABLE: foreach my $spots(sort keys %score){
            my @spots=split /,/,$spots;
            foreach my $spot(@spots){             
                    my ($core_x,$core_y)=(split /[:_]/,$spot)[-2,-1];
                    my $sample=(split /:/,$spot)[0];
                    foreach my  $spots2(sort keys %score){  #: print "$spots\t$spots2\n";
                          my @spots2=split /,/,$spots2;
                          foreach my $spot2(@spots2){
                                if($spot ne $spot2){
                                     my ($core_x2,$core_y2)=(split /[:_]/,$spot2)[-2,-1];
                                     my @x=($core_x-20,$core_x,$core_x+20);
                                     my @y=($core_y-20,$core_y,$core_y+20);
                                     foreach my $x(@x){
                                        foreach my $y(@y) {
                                              if ($x==$core_x2 && $y==$core_y2){
                                                     my $s="$spots,$spots2";
                                                     my @s=split /,/,$s;
                                                     my @s_u=uniq(@s);
                                                     my $id=join(",",@s_u);
                                                     $score{$id}=1;   
                                                     delete $score{$spots};   
                                                     delete $score{$spots2};
                                                     next LABLE;
                                              }                                               
                                        }
                                     }  
                                }
                          }
                   }
           }
      }

open OUT,">$out" or die;
foreach my $s(sort keys %score){
       print "$s\n";
       my %spot;
       my @spot=split /,/,$s;my $num=$#spot+1;
       my @x;my @y;
       foreach my $spot(@spot){
              my ($x,$y)=(split /[:_]/,$spot)[-2,-1];
              push @x,$x;
              push @y,$y;              
       } 
       my $x_mean=int(aver(\@x));
       my $y_mean=int(aver(\@y));
       foreach my $spot2(@spot){
             my ($x,$y)=(split /[:_]/,$spot2)[-2,-1];
             my $distance=sqrt( ($y-$y_mean)*($y-$y_mean)+($x-$x_mean)*($x-$x_mean) );
             $spot{"distance"}{$spot2}=$distance;
             $spot{"OO.genes1"}{$spot2}=$score_raw{$spot2}{"OO.genes1"};
             $spot{"GC.genes1"}{$spot2}=$score_raw{$spot2}{"GC.genes1"};
             $spot{"TC.genes1"}{$spot2}=$score_raw{$spot2}{"TC.genes1"};
    
       }
       
       my @spot_sort=sort {$spot{"distance"}{$b} <=> $spot{"distance"}{$a}} keys %{$spot{"distance"}};
       my @score_oo;my @score_gc;my @score_tc;my @dist;
       foreach my $spot_sort(@spot_sort){
              push @dist, $spot{"distance"}{$spot_sort};
              push @score_oo,$spot{"OO.genes1"}{$spot_sort};
              push @score_gc,$spot{"GC.genes1"}{$spot_sort};
              push @score_tc,$spot{"TC.genes1"}{$spot_sort};
       }
       
         
       if($num >8){
              foreach my $num2(1..4){
                     if($score_gc[-$num2] < $score_tc[-$num2] && $score_oo[-$num2] < $score_tc[-$num2]){
                              $dist[-$num2]=(5-$num2)*200;
                     }
               }

               my @outlier_num;
               foreach my $num3(1..($#spot_sort)){
                      if(abs($dist[$num3]-$dist[$num3-1]) >20){
                             push @outlier_num,$num3;
                      }
               }
               my $outlier_num=min(@outlier_num);
               if($outlier_num[0]){
                    @dist=@dist[0..($outlier_num-1)];#@dist=@dist_n;
                    @spot_sort=@spot_sort[0..($outlier_num-1)];#@id=@id_n;
               }
      }
      foreach my $s(@spot_sort){
           my ($sam,$x,$y)=(split /[:_]/,$s)[0,-2,-1];
      }
      my $spot_sort=join(",",@spot_sort);
      my $spot_sort_num=$#spot_sort+1;
      print OUT "$x_mean\t$y_mean\t$spot_sort_num\t$spot_sort\n";
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

