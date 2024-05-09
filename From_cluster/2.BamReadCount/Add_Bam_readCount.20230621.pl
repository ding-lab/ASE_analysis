use strict;
use List::Util qw(max);

my $dir=$ARGV[0];
my $type=$ARGV[1];
my $chr_prefix="chr";

open OUT, ">BamReadCount_added/$type/$dir.BamReadCount";
open (INPUT, "Site_lists/All_sites.tsv") or die "Can't open file: $!";
my $c_str=0;

while (<INPUT>){
    chomp;
    $c_str++;
    my $str=$_;
    if ($c_str==1){
	next;
    }	
    my @arr=split /\t/, $str;
    my $sample=$arr[2];
    my $chr=$arr[0];
    my $pos=$arr[1];
    my $ref;
    
    my $ref_s=$arr[3];
    my $alt_s=$arr[4];
#   my $alt=$arr[4];
#    if (length($ref_ch)>1 | length($alt_ch)>1){
#	next;
#    }
    
#    my ($ref_c_t, $alt_c_t, $vaf_t);

    open READC_T, "$type/$dir/$sample.read_counts.tsv.parsed.tsv";
    while (<READC_T>){
	chomp;
	my $str_r=$_;
	my @arr_r=split /\t/, $str_r;
	if ($arr_r[0] eq $chr && $arr_r[1]==$pos){
	    $ref=$arr_r[2];
	    if ($ref ne $ref_s) {
		print "$sample\t$chr\t$ref_s\t$ref\tCHECK REF\n";
	    }

	    my ($ref_c, $alt_c);
	    if ($ref_s eq $arr_r[4]){
		$ref_c=$arr_r[5];
	    }
	    if ($ref_s eq $arr_r[6]){
		$ref_c=$arr_r[7];
	    }
	    if ($ref_s eq $arr_r[8]){
		$ref_c=$arr_r[9];
	    }
	    if ($ref_s eq $arr_r[10]){
		$ref_c=$arr_r[11];
	    }
	    if ($alt_s eq $arr_r[4]){
		$alt_c=$arr_r[5];
	    }
	    if ($alt_s eq $arr_r[6]){
		$alt_c=$arr_r[7];
	    }
	    if ($alt_s eq $arr_r[8]){
		$alt_c=$arr_r[9];
	    }
	    if ($alt_s eq $arr_r[10]){
		$alt_c=$arr_r[11];
	    }

	    
	    print OUT "$str\t$ref_c\t$alt_c\n";

	}
    }
    
    
}
	    

	    
	    
	    
