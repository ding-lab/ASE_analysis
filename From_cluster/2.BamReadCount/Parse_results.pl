use strict;

open FILE, "$ARGV[0]";
open OUT, ">$ARGV[0].parsed.tsv";
while (<FILE>){
    chomp;
    my $str=$_;
    my @arr=split /\t/, $str;
    for (0..2){
	print OUT "$arr[$_]\t";
    }
    print OUT "$arr[3]";
    for (5..$#arr){
	my $ind=$_;
	my @arr_2=split /\:/, $arr[$ind];
	print OUT "\t$arr_2[0]\t$arr_2[1]";
    }
    print OUT "\n";
}
