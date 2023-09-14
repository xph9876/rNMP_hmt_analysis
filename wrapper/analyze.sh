#!/bin/bash

RibosePreferenceAnalysis='/storage/home/hcoda1/0/pxu64/bio-storici/scripts/RibosePreferenceAnalysis'
scripts='/storage/home/hcoda1/0/pxu64/bio-storici/scripts/rNMP_hmt_analysis'
bed_folder='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt/bed'
order='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt/order_human.tsv'
output='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt/results'
wrapper="$scripts/wrapper"

# make folders
if [ ! -d $output ]
then
    mkdir $output
fi

# Remove plots folder
if [ -d $output/plots ]
then
    rm -rf $output/plots
fi

for folder in heatmap_barplot logs enriched_zone control_region_figures gene_analysis plots bed_split strand_split
do
    if ! [ -d $output/$folder ]
    then 
        mkdir $output/$folder
    fi
done


# # Count rNMPs in each library
# curr=$(pwd)
# cd $bed_folder
# wc -l *.bed | grep -v total | sed 's/.bed//;s/^ *//;s/ /\t/' | awk 'BEGIN{OFS="\t"}{print $2,$1}' > $output/mito_count.tsv
# cd $scripts/gene_analysis/control_bed
# wc -l *.bed | grep -v total | sed 's/.bed//;s/^ *//;s/ /\t/' | awk 'BEGIN{OFS="\t"}{print $2,$1}' >> $output/mito_count.tsv
# cd $curr

# # Generate strand split bed files
# for file in $(ls $bed_folder)
# do
#     grep '+' $bed_folder/$file > $output/bed_split/${file::-4}_light.bed &
#     grep -v '+' $bed_folder/$file > $output/bed_split/${file::-4}_heavy.bed &
# done
# wait

# # heatmap
# eval $wrapper/heatmap_barplot.sh \
#     $RibosePreferenceAnalysis \
#     $scripts/heatmap_barplot \
#     $scripts/refseq/hg38_chrM.fa \
#     $bed_folder \
#     $order \
#     $output/heatmap_barplot \
#     2>&1 > $output/logs/heatmap_barplot.log &

# # Strand split figures
# eval $wrapper/draw_strand_split.sh \
#     $RibosePreferenceAnalysis \
#     $scripts/heatmap_barplot \
#     $scripts/refseq \
#     $output/bed_split \
#     $order $output \
#     2>&1 > $output/logs/strand_split.log &

# # enriched zones
# eval $scripts/enriched_zone/enriched_zone_analysis.py \
#     $bed_folder/*.bed \
#     $scripts/refseq/hg38_chrM.fa.fai \
#     $order \
#     -o $output/enriched_zone/chrM \
#     2>&1 > $output/logs/enriched_zone.log &                                                        

# # Draw control region figures
# eval $wrapper/control_region_figures.sh \
#     $RibosePreferenceAnalysis \
#     $scripts/control_region \
#     $output/control_region_figures \
#     $bed_folder \
#     $order \
#     $output/mito_count.tsv \
#     $scripts/refseq/hg38_chrM.fa \
#     2>&1 > $output/logs/control_region_figures.log &

# Perform gene analysis
eval $wrapper/gene_analysis.sh \
    $scripts/gene_analysis \
    $RibosePreferenceAnalysis \
    $output/gene_analysis \
    $bed_folder \
    $scripts/gene_analysis/order_human.tsv \
    $output/mito_count.tsv > $output/logs/gene_analysis.log &

wait

# for folder in heatmap_barplot control_region_figures gene_analysis strand_split
# do
#     mv $output/$folder/plots $output/plots/$folder
# done
# mv $output/enriched_zone $output/plots/enriched_zone

echo 'Done!'