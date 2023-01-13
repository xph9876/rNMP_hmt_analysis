#!/bin/bash

RibosePreferenceAnalysis='/storage/home/hcoda1/0/pxu64/bio-storici/scripts/RibosePreferenceAnalysis'
scripts='/storage/home/hcoda1/0/pxu64/bio-storici/scripts/rNMP_hmt_analysis'
genome='/storage/home/hcoda1/0/pxu64/bio-storici/genome/hg38_chrM.fa'
bed_folder='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt_new/raw_bed'
order='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt_new/order_human.tsv'
output='/storage/home/hcoda1/0/pxu64/bio-storici/human_mt_new/results_raw'
wrapper="$scripts/wrapper"

# make folders
if [ ! -d $output ]
then
    mkdir $output
fi
for folder in heatmap_barplot logs enriched_zone control_region_figures gene_analysis plots
do
    if ! [ -d $output/$folder ]
    then 
        mkdir $output/$folder
    fi
done


# Count rNMPs in each library
curr=$(pwd)
cd $bed_folder
wc -l *.bed | grep -v total | sed 's/.bed//;s/^ *//;s/ /\t/' | awk 'BEGIN{OFS="\t"}{print $2,$1}' > $output/mito_count.tsv
cd $curr

# heatmap
eval $wrapper/heatmap_barplot.sh $RibosePreferenceAnalysis $scripts/heatmap_barplot $genome $bed_folder $order $output/heatmap_barplot 2>&1 > $output/logs/heatmap_barplot.log &

# enriched zones
eval $scripts/enriched_zone/enriched_zone_analysis.py $bed_folder/*.bed ${genome}.fai $order -o $output/enriched_zone/chrM 2>&1 > $output/logs/enriched_zone.log &                                                        

# Draw control region figures
eval $wrapper/control_region_figures.sh $RibosePreferenceAnalysis $scripts/control_region $output/control_region_figures $bed_folder $order $output/mito_count.tsv $genome 2>&1 > $output/logs/control_region_figures.log &

# Perform gene analysis
eval $wrapper/gene_analysis.sh $scripts/gene_analysis $output/gene_analysis $bed_folder $order $output/mito_count.tsv > $output/logs/gene_analysis.log &

wait

for folder in heatmap_barplot control_region_figures gene_analysis
do
    mv $output/$folder/plots $output/plots/$folder
done
mv $output/enriched_zone $output/plots/enriched_zone

echo 'Done!'