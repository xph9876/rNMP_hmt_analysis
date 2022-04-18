#!/bin/bash

if ! [ -e strand_split ]
then
    mkdir strand_split
fi

for folder in individual raw normalized tsv plots bg
do
    if ! [ -e strand_split/$folder ]
    then
        mkdir strand_split/$folder
    fi
done

# background forward
~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa --mono -s -o strand_split/bg/chrM_mono_forward.raw &
for dis in $(seq 1 5) 100
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa -d ${dis} -s -o strand_split/bg/chrM_dinuc_d${dis}_forward.raw &
done
~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa --trinuc -s -o strand_split/bg/chrM_trinuc_forward.raw &


# background reverse
~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM_rc.fa --mono -s -o strand_split/bg/chrM_mono_reverse.raw &

for dis in $(seq 1 5) 100
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM_rc.fa -d ${dis} -s -o strand_split/bg/chrM_dinuc_d${dis}_reverse.raw &
done
~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM_rc.fa --trinuc -s -o strand_split/bg/chrM_trinuc_reverse.raw &
wait

# count individual
for file in $(ls ../bed_split)
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/count_rNMP.py ../../genome/hg38_chrM.fa ../bed_split/$file -m -d --dist 1 2 3 4 5 100 -t -o strand_split/individual/chrM &
done
wait
rename 's/chrM_//' strand_split/individual/* -f

# Get all files
for st in forward reverse
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py strand_split/individual/*${st}*mono -o strand_split/raw/chrM_mono_${st}.raw &
    for dis in 1 2 3 4 5 100
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py strand_split/individual/*${st}*d${dis}_nr -o strand_split/raw/chrM_dinuc_d${dis}_nr_${st}.raw &
        ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py strand_split/individual/*${st}*d${dis}_rn -o strand_split/raw/chrM_dinuc_d${dis}_rn_${st}.raw &
    done
    for pat in nnr nrn rnn 
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py strand_split/individual/*${st}*${pat} -o strand_split/raw/chrM_trinuc_${pat}_${st}.raw &
    done
done
wait
sed -i 's/_forward//;s/_reverse//' strand_split/raw/*

# normalize
for st in forward reverse
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py strand_split/raw/chrM_mono_${st}.raw strand_split/bg/chrM_mono_${st}.raw --name chrM -o strand_split/normalized/chrM_mono_${st}.norm &
    for dis in 1 2 3 4 5 100
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py strand_split/raw/chrM_dinuc_d${dis}_nr_${st}.raw strand_split/bg/chrM_dinuc_d${dis}_${st}.raw --name chrM --group_len 4 -o strand_split/normalized/chrM_dinuc_d${dis}_nr_${st}.norm &
        ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py strand_split/raw/chrM_dinuc_d${dis}_rn_${st}.raw strand_split/bg/chrM_dinuc_d${dis}_${st}.raw --name chrM --group_len 4 -o strand_split/normalized/chrM_dinuc_d${dis}_rn_${st}.norm &
    done
    for pat in nnr nrn rnn
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py strand_split/raw/chrM_trinuc_${pat}_${st}.raw strand_split/bg/chrM_trinuc_${st}.raw --name chrM --group_len 16 -o strand_split/normalized/chrM_trinuc_${pat}_${st}.norm &
    done
done
wait

# resort
ls strand_split/normalized | xargs -I aa -P 8 ~/bio-storici/scripts/RibosePreferenceAnalysis/resort.py strand_split/normalized/aa ../order_human.tsv -c 2 -o strand_split/tsv/aa.tsv
rename 's/norm.tsv/tsv/' strand_split/tsv/* -f

# draw
for st in forward reverse
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py strand_split/tsv/chrM_mono_${st}.tsv -b strand_split/bg/chrM_mono_${st}.raw -o strand_split/plots/chrM_mono_${st}.png --palette RdBu_r &
    for dis in 1 2 3 4 5 100
    do
        for ty in nr rn
        do
            ~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py strand_split/tsv/chrM_dinuc_d${dis}_${ty}_${st}.tsv -b strand_split/bg/chrM_dinuc_d${dis}_${st}.raw -o strand_split/plots/chrM_dinuc_d${dis}_${ty}_${st}.png --palette RdBu_r &
        done
    done
    for pat in nnr nrn rnn
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py strand_split/tsv/chrM_trinuc_${pat}_${st}.tsv -b strand_split/bg/chrM_trinuc_${st}.raw -o strand_split/plots/chrM_trinuc_${pat}_${st}.png --no_annot --palette RdBu_r &
    done
done
wait


# draw barplot
./generate_bar_plot.py strand_split/tsv/chrM_mono_${st}.tsv -o strand_split/plots/chrM_barplot_normalized_${st}.png &
~/bio-storici/scripts/RibosePreferenceAnalysis/sum1.py strand_split/raw/chrM_mono_${st}.raw -o strand_split/normalized/chrM_sum1_${st}.norm
~/bio-storici/scripts/RibosePreferenceAnalysis/resort.py strand_split/normalized/chrM_sum1_${st}.norm ../order_human.tsv -c 2 -o strand_split/tsv/chrM_mono_sum1_${st}.tsv
./generate_bar_plot.py strand_split/tsv/chrM_mono_sum1_${st}.tsv -o strand_split/plots/chrM_barplot_raw_${st}.png &

# draw barplot light/heavy comparison
ls strand_split/raw | grep mono| xargs -I aa -P 2 ~/bio-storici/scripts/RibosePreferenceAnalysis/resort.py strand_split/raw/aa ../order_human.tsv -c 2 -o strand_split/tsv/aa.tsv
rename 's/.raw.tsv/_raw.tsv/' strand_split/tsv/* -f
./strand_split_plot.py strand_split/tsv/chrM_mono_forward_raw.tsv strand_split/tsv/chrM_mono_reverse_raw.tsv -o strand_split/plots/chrM_count.png

