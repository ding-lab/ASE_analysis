use strict;


my @all_cases;
my $c_case=0;
my $r_dir=$ARGV[0];

open FILE, "Site_lists/$r_dir/Sample_list.tsv";
while (<FILE>){
    chomp;
    my $str=$_;
    push @all_cases, $str;
    $c_case++;
}

for (@all_cases){
    my $case=$_;
    
    `bash run_bamreadcount.sh -I BAM_paths_NAT/out/$r_dir/$case.RNA-Seq.genomic.N.hg38.bam -l Site_lists/$r_dir/$case.sites.tsv -R /diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa -o $case -O NAT/$r_dir/`;    

}

print "\nNumber of processed samples: $c_case\n";
