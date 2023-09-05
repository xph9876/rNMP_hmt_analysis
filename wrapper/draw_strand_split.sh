#!/bin/bash

# Draw heatmaps and barplots
# Usage ./draw_strand_split.sh <RibosePreferenceAnalysis location> <heatmap_barplot folder location> <genome folder> <bed split folder> <library order> <output folder> 
RibosePreferenceAnalysis=$1
heatmap_barplot=$2
genome_folder=$3
genome=$genome_folder/hg38_chrM.fa
genome_rc=$genome_folder/hg38_chrM_rc.fa
bed_split_folder=$4
order=$5
output=$6


# Make output folders
if ! [ -d $output/strand_split ]
then
    mkdir $output/strand_split
fi

for folder in individual raw normalized tsv plots bg
do
    if ! [ -d $output/strand_split/$folder ]
    then
        mkdir $output/strand_split/$folder
    fi
done


# background light
$RibosePreferenceAnalysis/count_background.py $genome --mono -s -o $output/strand_split/bg/chrM_mono_light.raw &
for dis in $(seq 1 5) 100
do
    $RibosePreferenceAnalysis/count_background.py $genome -d ${dis} -s -o $output/strand_split/bg/chrM_dinuc_d${dis}_light.raw &
done
$RibosePreferenceAnalysis/count_background.py $genome --trinuc -s -o $output/strand_split/bg/chrM_trinuc_light.raw &


# background heavy
$RibosePreferenceAnalysis/count_background.py $genome_rc --mono -s -o $output/strand_split/bg/chrM_mono_heavy.raw &

for dis in $(seq 1 5) 100
do
    $RibosePreferenceAnalysis/count_background.py $genome_rc -d ${dis} -s -o $output/strand_split/bg/chrM_dinuc_d${dis}_heavy.raw &
done
$RibosePreferenceAnalysis/count_background.py $genome_rc --trinuc -s -o $output/strand_split/bg/chrM_trinuc_heavy.raw &
wait

# count individual
for file in $(ls $bed_split_folder)
do
    $RibosePreferenceAnalysis/count_rNMP.py $genome $bed_split_folder/$file -m -d --dist 1 2 3 4 5 100 -t -o $output/strand_split/individual/chrM &
done
wait
rename 's/chrM_//' $output/strand_split/individual/* -f

# Get all files
for st in light heavy
do
    $RibosePreferenceAnalysis/get_chrom.py $output/strand_split/individual/*${st}*mono -o $output/strand_split/raw/chrM_mono_${st}.raw &
    for dis in 1 2 3 4 5 100
    do
        $RibosePreferenceAnalysis/get_chrom.py $output/strand_split/individual/*${st}*d${dis}_nr -o $output/strand_split/raw/chrM_dinuc_d${dis}_nr_${st}.raw &
        $RibosePreferenceAnalysis/get_chrom.py $output/strand_split/individual/*${st}*d${dis}_rn -o $output/strand_split/raw/chrM_dinuc_d${dis}_rn_${st}.raw &
    done
    for pat in nnr nrn rnn 
    do
        $RibosePreferenceAnalysis/get_chrom.py $output/strand_split/individual/*${st}*${pat} -o $output/strand_split/raw/chrM_trinuc_${pat}_${st}.raw &
    done
done
wait
sed -i 's/_light//;s/_heavy//' $output/strand_split/raw/*

# normalize
for st in light heavy
do
    $RibosePreferenceAnalysis/normalize.py $output/strand_split/raw/chrM_mono_${st}.raw $output/strand_split/bg/chrM_mono_${st}.raw --name chrM -o $output/strand_split/normalized/chrM_mono_${st}.norm &
    for dis in 1 2 3 4 5 100
    do
        $RibosePreferenceAnalysis/normalize.py $output/strand_split/raw/chrM_dinuc_d${dis}_nr_${st}.raw $output/strand_split/bg/chrM_dinuc_d${dis}_${st}.raw --name chrM --group_len 4 -o $output/strand_split/normalized/chrM_dinuc_d${dis}_nr_${st}.norm &
        $RibosePreferenceAnalysis/normalize.py $output/strand_split/raw/chrM_dinuc_d${dis}_rn_${st}.raw $output/strand_split/bg/chrM_dinuc_d${dis}_${st}.raw --name chrM --group_len 4 -o $output/strand_split/normalized/chrM_dinuc_d${dis}_rn_${st}.norm &
    done
    for pat in nnr nrn rnn
    do
        $RibosePreferenceAnalysis/normalize.py $output/strand_split/raw/chrM_trinuc_${pat}_${st}.raw $output/strand_split/bg/chrM_trinuc_${st}.raw --name chrM --group_len 16 -o $output/strand_split/normalized/chrM_trinuc_${pat}_${st}.norm &
    done
done
wait

# resort
ls $output/strand_split/normalized | xargs -I aa -P 8 $RibosePreferenceAnalysis/resort.py $output/strand_split/normalized/aa $order -c 2 -o $output/strand_split/tsv/aa.tsv
rename 's/norm.tsv/tsv/' $output/strand_split/tsv/* -f

# draw
for st in light heavy
do
    $RibosePreferenceAnalysis/draw_heatmap.py $output/strand_split/tsv/chrM_mono_${st}.tsv -b $output/strand_split/bg/chrM_mono_${st}.raw -o $output/strand_split/plots/chrM_mono_${st}.png --palette RdBu_r &
    for dis in 1 2 3 4 5 100
    do
        for ty in nr rn
        do
            $RibosePreferenceAnalysis/draw_heatmap.py $output/strand_split/tsv/chrM_dinuc_d${dis}_${ty}_${st}.tsv -b $output/strand_split/bg/chrM_dinuc_d${dis}_${st}.raw -o $output/strand_split/plots/chrM_dinuc_d${dis}_${ty}_${st}.png --palette RdBu_r &
        done
    done
    for pat in nnr nrn rnn
    do
        $RibosePreferenceAnalysis/draw_heatmap.py $output/strand_split/tsv/chrM_trinuc_${pat}_${st}.tsv -b $output/strand_split/bg/chrM_trinuc_${st}.raw -o $output/strand_split/plots/chrM_trinuc_${pat}_${st}.png --no_annot --palette RdBu_r &
    done
done
wait


# draw barplot
for st in light heavy
do
    $heatmap_barplot/generate_bar_plot.py $output/strand_split/tsv/chrM_mono_${st}.tsv -o $output/strand_split/plots/chrM_barplot_normalized_${st}.png &
    $RibosePreferenceAnalysis/sum1.py $output/strand_split/raw/chrM_mono_${st}.raw -o $output/strand_split/normalized/chrM_sum1_${st}.norm
    $RibosePreferenceAnalysis/resort.py $output/strand_split/normalized/chrM_sum1_${st}.norm $order -c 2 -o $output/strand_split/tsv/chrM_mono_sum1_${st}.tsv
    $heatmap_barplot/generate_bar_plot.py $output/strand_split/tsv/chrM_mono_sum1_${st}.tsv -o $output/strand_split/plots/chrM_barplot_raw_${st}.png &
done
wait

# draw barplot light/heavy comparison
ls $output/strand_split/raw | grep mono| xargs -I aa -P 2 $RibosePreferenceAnalysis/resort.py $output/strand_split/raw/aa $order -c 2 -o $output/strand_split/tsv/aa.tsv
rename 's/.raw.tsv/_raw.tsv/' $output/strand_split/tsv/* -f
$heatmap_barplot/strand_split_plot.py $output/strand_split/tsv/chrM_mono_light_raw.tsv $output/strand_split/tsv/chrM_mono_heavy_raw.tsv -o $output/strand_split/plots/chrM_count.png --no_annot
$heatmap_barplot/strand_split_plot.py $output/strand_split/tsv/chrM_mono_light_raw.tsv $output/strand_split/tsv/chrM_mono_heavy_raw.tsv -o $output/strand_split/plots/chrM_count_annot.png

# Draw contribution for strand split
$heatmap_barplot/strand_split_contribution.py \
    $output/strand_split/raw/chrM_mono_light.raw \
    $output/strand_split/raw/chrM_mono_heavy.raw \
    $output/strand_split/bg/chrM_mono_light.raw \
    $output/strand_split/bg/chrM_mono_heavy.raw \
    $order -o $output/strand_split/plots/strand_split

