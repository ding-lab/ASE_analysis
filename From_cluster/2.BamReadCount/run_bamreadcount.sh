#2022-09-29: Updated for katmai
#!/bin/bash
# ------------------------------------------------------------------
#title                  :run_bamreadcount.sh
#description            :This script will run bam-readcount on a given .bam file using the positions provided.
#author                 :Fernanda Martins Rodrigues (fernanda)
#date                   :20190327
#usage                  :bash run_bamreadcount.sh -I [input bam] -l [site list] -R [reference genome to use] -o [output basename] -O [output_directory] -L [logs_directory]
#notes                  :Install bam-readcount to use this script
#bash_version           :4.2.46(2)-release (x86_64-redhat-linux-gnu)
#bam-readcount_version	:0.8
# ------------------------------------------------------------------

USAGE="Usage: bash run_bamreadcount.sh -I [input bam] -l [site list] -R [reference genome to use] -o [output basename] -O [output_directory]"
# --- Options processing -------------------------------------------

while getopts "I:l:R:o:O:" opt; do
        case $opt in
                I)
                        INPUT_BAM=$OPTARG # input BAM file
                        ;;
		l)
			SITE_LIST=$OPTARG # ile containing a list of regions to report readcounts within
			;;
		R)
			REF=$OPTARG #reference genome .fa file
			;;
                o)
                        OUT_BASENAME=$OPTARG # basename for output file
                        ;;
                O)
                        OUT_DIR=$OPTARG # output directory
                        ;;
                \?)
                        echo "Invalid option -$OPTARG"
                        echo $USAGE
                        exit 1
                        ;;
                :)
                        echo "Option -$OPTARG requires an argument."
                        echo $USAGE
                        exit 1
                        ;;
        esac
done

if [ "$#" -ne 10 ] ; then
        echo "Invalid number of arguments."
        echo $USAGE
        exit 1;
fi

# --- Body --------------------------------------------------------

if [ ! -d ${OUT_DIR} ]; then
	mkdir ${OUT_DIR}
fi

/diskmnt/Software/bam-readcount-0.7.4/mybuild/bin/bam-readcount -i -q 10 -b 15 -l ${SITE_LIST} -f ${REF} ${INPUT_BAM} > ${OUT_DIR}/${OUT_BASENAME}.read_counts.tsv 
#bam-readcount -q 10 -b 15 -l ${SITE_LIST} -f ${REF} ${INPUT_BAM} > ${OUT_DIR}/${OUT_BASENAME}.readcounts.tsv 

# --- End ---------------------------------------------------------








