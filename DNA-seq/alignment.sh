#!usr/bin/env bash

####################################
# DNA-seq alignment 
# Author: Deepali L. Kundnani
# Institute: Georgia institute of Technology(GT)
# Description: Trimming, alignment, sorting and indexing to get SAM/BAM files
# Dependencies: TrimGalore, Bowtie, Samtools, Preexisting bowtie2 indexes for Reference
##################################

##################################
# Edit the following parameters as needed before running the script
##################################

raw_reads='~/DNA-seq/FS41/' #Raw reads folder
trimmed_reads='~/DNA-seq/FS41/trimmed' #Foldername for trimmed reads
aligned='~/DNA-seq/FS41/aligned' #Foldername for aligned data
ref=$hgref #.fa file for reference genome

sample_list="CD4T-1 CD4T-2 CD4T-3 CD4T-4 DLTB-8 DLTB-P TLTB-8 TLTB-P H9-1 H9-2 HEK293T-RNASEH2A-KO-T17 HEK293T-RNASEH2A-KO-T8 HEK293T-WT" #Make sure you are using the identifies of fq files, such as 'sample_R1.fq' and 'sample_R2.fq'

##################################
# Donot modify below this line to make sure of consistent runs
##################################
mkdir $trimmed_reads $aligned $aligned/log $aligned/fixmate $aligned/temp
run() {
    trim_galore -q 15 --length 50 --fastqc --paired ${raw_reads}/${sample}_R1.fq ${raw_reads}/${sample}_R2.fq -o $trimmed_reads 
    bowtie2 --threads 4 -x $(dirname $ref)/$(basename $ref .fa) -1 ${trimmed_reads}/${sample}_R1_val_1.fq -2 ${trimmed_reads}/${sample}_R2_val_2.fq -S ${aligned}/${sample}.sam 2> ${aligned}/log/${sample}.log
    samtools sort -O bam -T $aligned/temp/ ${aligned}/fixmate/fix_${sample}.bam | samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' - -o ${aligned}/${sample}.bam
    samtools index ${aligned}/${sample}.bam
}

for sample in ${sample_list}
do 
run &
done
wait

