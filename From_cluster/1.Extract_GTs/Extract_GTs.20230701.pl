use strict;

my $file='/diskmnt/Projects/Users/nvterekhanova/Germline_analysis/DATA/MAF_files/CPTAC.PanCan1093.merged.AD5.noLongIndels.ROI.norm.sorted.vep100.gencode34.vcf.gz';
open (INPUT, "gunzip -c $file |");
open OUT, ">out/All_germline_var_sites.20230701.tsv";

my @arr_ids;
while (<INPUT>){
    chomp;
    my $str=$_;
    unless ($str=~/^#CHROM/ || $str=~/^chr/){
	next;
    }
    my (@all_ids, %gts);
    if ($str=~/#CHROM/){
	@arr_ids=split /\s+/, $str;
    }
    else{
    
	my @arr_str=split /\s+/, $str;
	for (9..$#arr_str){
	    if ($arr_str[$_] eq './.:.:.:.:.:.:.:.:.:.:.:.:.:.:.'){
		next;
	    }
	    else {
		my @arr_gts=split /\:/,  $arr_str[$_];

		
#		if ($arr_gts[0] eq '0/1' || $arr_gts[0] eq '1/1' || $arr_gts[0] eq '1/0'){ 
		if ($arr_gts[0] ne './.'){
		    print OUT "$arr_str[0]\t$arr_str[1]\t$arr_str[3]\t$arr_str[4]\t$_\t$arr_ids[$_]\t$arr_gts[0]\n";
		}
	    }
	}
    }
}
