#!/bin/bash


while IFS=\  read -r "BAMpwd" "CONDITION" "REPLICATE"; do
    if [[ $BAMpwd =~ bam ]]; then
	qsub -pe smp 1 -q public.q -cwd -b y -shell y -N map "./bam2fastq_map_UMIgroup.sh -v 10 -I ${BAMpwd} -O ${CONDITION}_${REPLICATE} -D ${RESULTS_DIR} -R ${REF}"
    fi
done < "$TARGETS_FILE"

