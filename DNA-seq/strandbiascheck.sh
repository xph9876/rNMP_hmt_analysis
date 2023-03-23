#!usr/bin/env bash

####################################
# DNA-seq alignment 
# Author: Deepali L. Kundnani
# Institute: Georgia institute of Technology(GT)
# Description: Strand Bias check for mitochondrial DNA
# Dependencies: Samtools
##################################

##################################
# Edit the following parameters as needed before running the script
##################################


raw_reads='~/DNA-seq/FS41/' #Raw reads folder
trimmed_reads='~/DNA-seq/FS41/trimmed' #Foldername for trimmed reads
aligned='~/DNA-seq/FS41/aligned' #Foldername for aligned data
ref=$hgref #.fa file for reference genome

sample_list="CD4T-1 CD4T-2 CD4T-3 CD4T-4 DLTB-8 DLTB-P TLTB-8 TLTB-P H9-1 H9-2 HEK293T-RNASEH2A-KO-T17 HEK293T-RNASEH2A-KO-T8 HEK293T-WT" #Make sure you are using the identifies of fq files, such as 'sample_R1.fq' and 'sample_R2.fq'

chr='chrM'

##################################
# Donot modify below this line to make sure of consistent runs
##################################

strand_bias() {
    forward=$(samtools view $file ${chr} | gawk '(and(16,$2))' | wc -l)
    reverse=$(samtools view $file ${chr} | gawk '(! and(16,$2))' | wc -l)
    f=$(basename $file .bam)
    echo $f'\t'$forward'\t'$reverse
}

strand_bias_all() {
    forward=$(samtools view $file | gawk '(and(16,$2))' | wc -l)
    reverse=$(samtools view $file | gawk '(! and(16,$2))' | wc -l)
    f=$(basename $file .bam)
    echo $f $forward $reverse
}

echo "Strand bias in selected chromosome"

for sample in sample_list
do
file=${aligned/${sample}.bam}
strand_bias 
done

echo "Strand bias in all chromosomes of reference"
for sample in sample_list
do
file=${aligned/${sample}.bam}
strand_bias_all
done