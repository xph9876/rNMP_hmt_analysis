#!/bin/bash

for folder in individual raw normalized tsv plots bg
do
    if ! [ -e $folder ]
    then
        mkdir $folder
    fi
done

# # background
# ~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa --mono -o bg/chrM_mono.raw &
# for dis in $(seq 1 5) 100
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa -d ${dis} -o bg/chrM_dinuc_d${dis}.raw &
# done
# ~/bio-storici/scripts/RibosePreferenceAnalysis/count_background.py ../../genome/hg38_chrM.fa --trinuc -o bg/chrM_trinuc.raw &
# wait

# # count individual
# for file in $(ls ../bed)
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/count_rNMP.py ../../genome/hg38_chrM.fa ../bed/$file -m -d --dist 1 2 3 4 5 100 -t -o individual/chrM &
# done
# wait
# rename 's/chrM_//' individual/* -f

# # Get all files
# ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py individual/*mono -o raw/chrM_mono.raw &
# for dis in 1 2 3 4 5 100
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py individual/*d${dis}_nr -o raw/chrM_dinuc_d${dis}_nr.raw &
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py individual/*d${dis}_rn -o raw/chrM_dinuc_d${dis}_rn.raw &
# done
# for pat in nnr nrn rnn 
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/get_chrom.py individual/*${pat} -o raw/chrM_trinuc_${pat}.raw &
# done
# wait

# # normalize
# ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py raw/chrM_mono.raw bg/chrM_mono.raw --name chrM -o normalized/chrM_mono.norm &
# for dis in 1 2 3 4 5 100
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py raw/chrM_dinuc_d${dis}_nr.raw bg/chrM_dinuc_d${dis}.raw --name chrM --group_len 4 -o normalized/chrM_dinuc_d${dis}_nr.norm &
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py raw/chrM_dinuc_d${dis}_rn.raw bg/chrM_dinuc_d${dis}.raw --name chrM --group_len 4 -o normalized/chrM_dinuc_d${dis}_rn.norm &
# done
# for pat in nnr nrn rnn
# do
#     ~/bio-storici/scripts/RibosePreferenceAnalysis/normalize.py raw/chrM_trinuc_${pat}.raw bg/chrM_trinuc.raw --name chrM --group_len 16 -o normalized/chrM_trinuc_${pat}.norm &
# done
# wait

# # resort
# ls normalized | xargs -I aa -P 8 ~/bio-storici/scripts/RibosePreferenceAnalysis/resort.py normalized/aa ../order_human.tsv -c 2 -o tsv/aa.tsv
# rename 's/norm.tsv/tsv/' tsv/* -f

# draw
~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py tsv/chrM_mono.tsv -b bg/chrM_mono.raw -o plots/chrM_mono.png --palette RdBu_r &
for dis in 1 2 3 4 5 100
do
    for ty in nr rn
    do
        ~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py tsv/chrM_dinuc_d${dis}_${ty}.tsv -b bg/chrM_dinuc_d${dis}.raw -o plots/chrM_dinuc_d${dis}_${ty}.png --palette RdBu_r &
    done
done
for pat in nnr nrn rnn
do
    ~/bio-storici/scripts/RibosePreferenceAnalysis/draw_heatmap.py tsv/chrM_trinuc_${pat}.tsv -b bg/chrM_trinuc.raw -o plots/chrM_trinuc_${pat}.png --no_annot --palette RdBu_r &
done
wait


# # draw barplot
# ./generate_bar_plot.py tsv/chrM_mono.tsv -o plots/chrM_barplot_normalized.png &
# ~/bio-storici/scripts/RibosePreferenceAnalysis/sum1.py raw/chrM_mono.raw -o normalized/chrM_sum1.norm
# ~/bio-storici/scripts/RibosePreferenceAnalysis/resort.py normalized/chrM_sum1.norm ../order_human.tsv -c 2 -o tsv/chrM_mono_sum1.tsv
# ./generate_bar_plot.py tsv/chrM_mono_sum1.tsv -o plots/chrM_barplot_raw.png &

