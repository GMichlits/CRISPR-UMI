#!/usr/bin/env bash

module load samtools/1.3.1
samtools faidx ${REF}

module load bowtie/1.2
bowtie-build ${REF} ${REF}

