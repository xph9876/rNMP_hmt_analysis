#!/bin/bash

scripts="~/bio-storici/scripts/rNMP_hmt_analysis/replication_origin"
heatmap="~/bio-storici/scripts/RibosePreferenceAnalysis"

for aa in bed_same bed_oppo plots same oppo
do
    if ! [ -e $aa ]
    then
        mkdir $aa
    fi
done

for aa in distribution same oppo
do
    if ! [ -e plots/$aa ]
    then
        mkdir plots/$aa
    fi
done


# please copy bg!!!
if ! [ -e bg ]
then
    echo 'Please copy bg from github!!!'
    exit 1
fi


# generate intersect
for file in $(ls ../bed/ | grep FS)
do
    bedtools intersect -a bg/same.bed -b ../bed/${file} -s > bed_same/$file &
    bedtools intersect -a bg/oppo.bed -b ../bed/${file} -s > bed_oppo/$file &
done
wait

# analysis
for st in same oppo
do
    # seperate cr oh and ol
    for region in CR OH OL
    do
        # create folders
        for folder in bed mono
        do
            if ! [ -e ${st}/${folder}_${region} ]
            then
                mkdir ${st}/${folder}_${region}
            fi
        done

        # seperate
        for lib in $(ls bed_${st})
        do
            grep ${region} bed_${st}/$lib > ${st}/bed_${region}/$lib &
        done
        wait

        # count mono
        eval $heatmap/count_rNMP.py ../../genome/hg38_chrM.fa ${st}/bed_${region}/* -m -o ${st}/mono_${region}/test &
        wait
        rename s'/test_//' ${st}/mono_${region}/* -f

        # get chrom
        eval $heatmap/get_chrom.py ${st}/mono_${region}/* -o ${st}/${region}_mono.raw
        eval $heatmap/normalize.py ${st}/${region}_mono.raw bg/${st}.tsv --name ${region} -o ${st}/${region}_mono.norm
        eval $heatmap/resort.py ${st}/${region}_mono.norm ../order_human.tsv -c 2 -o ${st}/${region}_mono.tsv
        eval $heatmap/draw_heatmap.py ${st}/${region}_mono.tsv -b bg/${st}.tsv --background_chrom ${region} --palette RdBu_r -o plots/${st}/${region}

        # barplot
        eval $scripts/../heatmap_barplot/generate_bar_plot.py ${st}/${region}_mono.tsv -o plots/${st}/${region}_barplot &
    done
done

# control region distribution
eval $scripts/replication_distribution.py ../mito_count.tsv ../order_human.tsv --same same/bed_CR/* --oppo oppo/bed_CR/* -o plots/distribution/CR

for aa in OH OL CR
do
    eval $scripts/compare_strands.py same/${aa}_mono.raw oppo/${aa}_mono.raw ../order_human.tsv -o plots/${aa}_strands
    
# done

eval $scripts/replication_distribution_combined.py ../mito_count.tsv ../order_human.tsv --same same/bed_CR/* --oppo oppo/bed_CR/* -o plots/distribution/all_wt
eval $scripts/replication_distribution_combined.py ../mito_count.tsv ../order_human.tsv --same same/bed_CR/* --oppo oppo/bed_CR/* --selected "HEK293T" "RNH2A-KO\ T3-8" "RNH2A-KO\ T3-17" --palette Set2 -o plots/distribution/all_hek

wait