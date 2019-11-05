#!/bin/bash



while getopts I:O:D:v:R: option
do
 case "${option}"
 in
 I) INPUT_BAM=${OPTARG};;
 O) BNAME=${OPTARG};;
 D) OUT_DIR=${OPTARG};;
 v) NTHREADS=${OPTARG};;
 R) REFERENCE=${OPTARG};;
 *) echo "Unknown option";;
 esac
done



module load samtools/1.3.1
module load bowtie/1.2

mkdir -p ${OUT_DIR}


# convert unaligned bam into fastq, while adding UMI to readid for usage with umi_tools, and excluding reads with N in UMI

samtools view ${INPUT_BAM} |  awk '{if(!match(substr($14,6,10),/N/)){print "@"$1":"substr($14,6,10)"\n"substr($10,1,20)"\n+\n"substr($11,1,20)}}' > ${OUT_DIR}/${BNAME}.fastq



# use bowtie mapping and filtering to retain only perfect matches to guides, add header after filtering

#  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
#  -m <int>           suppress all alignments if > <int> exist (def: no limit)
#  --best             hits guaranteed best stratum; ties broken by quality
#  --strata           hits in sub-optimal strata aren't reported (requires --best)
# 0x0004 Query unmapped
bowtie -p ${NTHREADS} -m 1 --best --strata -S -v 0 ${REFERENCE} ${OUT_DIR}/${BNAME}.fastq | samtools view -S -F 0x0004 - |  samtools view -ht ${REFERENCE}.fai - | samtools sort --threads ${NTHREADS} --output-fmt BAM -o ${OUT_DIR}/${BNAME}.bowtie.bam -

samtools index  ${OUT_DIR}/${BNAME}.bowtie.bam


# umi_tools group flatfile
#    read_id
#    contig
#    position
#    umi = raw umi
#    umi_count = how many times was this umi observed at the same alignment coordinates
#    final_umi = the error corrected umi
#    final_umi_count = how many times was the umi observed at the same alignment coordinates, inc. error correction
#    unique_id = the unique identifier for this group

#This is because adjacency only groups barcodes together if they are all separated from the same hub barcode and each by a single edge. Hence, any barcode with 2 or more errors will be always be in a separate group from the true barcode from which it originates. In contrast, the directional method is not limited to single edge paths.
#source activate py2.7; umi_tools group  --method=adjacency --umi-separator=:  --in-sam  -I  ${OUT_DIR}/${BNAME}.bowtie.bam  --group-out=${OUT_DIR}/${BNAME}_UMI_group.txt
source activate py2.7; umi_tools group  --umi-separator=:  --in-sam  -I  ${OUT_DIR}/${BNAME}.bowtie.bam  --group-out=${OUT_DIR}/${BNAME}_UMI_group.txt

##
awk '{if ( ($4==$6) || ($5>10 && $7>100 && $7>($5+10)) ){ print } }' ${OUT_DIR}/${BNAME}_UMI_group.txt | awk '{print $2,$6,$5,$4}' | sort | uniq > ${OUT_DIR}/${BNAME}_UMI_group_highconf.txt

awk '{arr[$1" "$2]+=$3} END {for (i in arr) {print i" "arr[i]}}'  ${OUT_DIR}/${BNAME}_UMI_group_highconf.txt >  ${OUT_DIR}/${BNAME}_guide_UMI_counts.txt
