#!/bin/bash

# Draw heatmaps and barplots
# Usage ./heatmap_barplot.sh <RibosePreferenceAnalysis location> <heatmap_barplot folder location> <chrM genome> <bed folder> <library order> <output folder> 
RibosePreferenceAnalysis=$1
heatmap_barplot=$2
genome=$3
bed_folder=$4
order=$5
output_folder=$6

# make folders
if ! [ -d $output_folder ]
then
    mkdir $output_folder
fi

for folder in individual raw normalized tsv plots bg
do
    if ! [ -d $output_folder/$folder ]
    then
        mkdir $output_folder/$folder
    fi
done

# background
eval $RibosePreferenceAnalysis/count_background.py $genome --mono -o $output_folder/bg/chrM_mono.raw &
for dis in $(seq 1 5) 100
do
    eval $RibosePreferenceAnalysis/count_background.py $genome -d ${dis} -o $output_folder/bg/chrM_dinuc_d${dis}.raw &
done
eval $RibosePreferenceAnalysis/count_background.py $genome --trinuc -o $output_folder/bg/chrM_trinuc.raw &
wait

# count individual
for file in $(ls $bed_folder)
do
    eval $RibosePreferenceAnalysis/count_rNMP.py $genome $bed_folder/$file -m -d --dist 1 2 3 4 5 100 -t -o $output_folder/individual/chrM &
done
wait
rename 's/chrM_//' $output_folder/individual/* -f

# Get all files
eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*mono -o $output_folder/raw/chrM_mono.raw &
for dis in 1 2 3 4 5 100
do
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*d${dis}_nr -o $output_folder/raw/chrM_dinuc_d${dis}_nr.raw &
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*d${dis}_rn -o  $output_folder/raw/chrM_dinuc_d${dis}_rn.raw &
done
for pat in nnr nrn rnn 
do
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*${pat} -o $output_folder/raw/chrM_trinuc_${pat}.raw &
done
wait

# normalize
eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/chrM_mono.raw $output_folder/bg/chrM_mono.raw --name chrM -o $output_folder/normalized/chrM_mono.norm &
for dis in 1 2 3 4 5 100
do
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/chrM_dinuc_d${dis}_nr.raw $output_folder/bg/chrM_dinuc_d${dis}.raw --name chrM --group_len 4 -o $output_folder/normalized/chrM_dinuc_d${dis}_nr.norm &
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/chrM_dinuc_d${dis}_rn.raw $output_folder/bg/chrM_dinuc_d${dis}.raw --name chrM --group_len 4 -o $output_folder/normalized/chrM_dinuc_d${dis}_rn.norm &
done
for pat in nnr nrn rnn
do
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/chrM_trinuc_${pat}.raw $output_folder/bg/chrM_trinuc.raw --name chrM --group_len 16 -o $output_folder/normalized/chrM_trinuc_${pat}.norm &
done
wait

# resort
ls $output_folder/normalized | xargs -I aa -P 8 $RibosePreferenceAnalysis/resort.py $output_folder/normalized/aa $order -c 2 -o $output_folder/tsv/aa.tsv
rename 's/norm.tsv/tsv/' $output_folder/tsv/* -f

# draw
eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/chrM_mono.tsv -b $output_folder/bg/chrM_mono.raw -o $output_folder/plots/chrM_mono.png --palette RdBu_r &
for dis in 1 2 3 4 5 100
do
    for ty in nr rn
    do
        eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/chrM_dinuc_d${dis}_${ty}.tsv -b $output_folder/bg/chrM_dinuc_d${dis}.raw -o $output_folder/plots/chrM_dinuc_d${dis}_${ty}.png --palette RdBu_r &
    done
done
for pat in nnr nrn rnn
do
    eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/chrM_trinuc_${pat}.tsv -b $output_folder/bg/chrM_trinuc.raw -o $output_folder/plots/chrM_trinuc_${pat}.png --no_annot --palette RdBu_r &
done
wait


# draw barplot
eval $heatmap_barplot/generate_bar_plot.py $output_folder/tsv/chrM_mono.tsv -o $output_folder/plots/chrM_barplot_normalized.png &
eval $RibosePreferenceAnalysis/sum1.py $output_folder/raw/chrM_mono.raw -o $output_folder/normalized/chrM_sum1.norm
eval $RibosePreferenceAnalysis/resort.py $output_folder/normalized/chrM_sum1.norm $order -c 2 -o $output_folder/tsv/chrM_mono_sum1.tsv
eval $heatmap_barplot/generate_bar_plot.py $output_folder/tsv/chrM_mono_sum1.tsv -o $output_folder/plots/chrM_barplot_raw.png &


for pat in nnr
do
    eval $heatmap_barplot/draw_heatmap.py $output_folder/tsv/chrM_trinuc_${pat}.tsv -b $output_folder/bg/chrM_trinuc.raw -o $output_folder/plots/chrM_trinuc_${pat}_for_crop.png --no_annot --palette RdBu_r &
done
wait