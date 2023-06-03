#!/bin/bash

# usage control_region_figures.sh <RibosePreferenceAnalysis> <control_region folder> <output folder> <bed folder> <order> <mito_count>
heatmap=$1
scripts=$2
output=$3
bg=$scripts/bg
chipseq=$scripts/chipseq
bed_folder=$4
order=$5
mito_count=$6
genome=$7

for aa in bed_same bed_oppo plots same oppo
do
    if ! [ -d $output/$aa ]
    then
        mkdir $output/$aa
    fi
done

for aa in distribution same oppo
do
    if ! [ -e $output/plots/$aa ]
    then
        mkdir $output/plots/$aa
    fi
done


# generate intersect
for file in $(ls $bed_folder | grep -v RD)
do
    bedtools intersect -a $bg/same.bed -b $bed_folder/${file} -s > $output/bed_same/$file &
    bedtools intersect -a $bg/oppo.bed -b $bed_folder/${file} -s > $output/bed_oppo/$file &
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
            if ! [ -d $output/${st}/${folder}_${region} ]
            then
                mkdir $output/${st}/${folder}_${region}
            fi
        done

        # seperate
        for lib in $(ls $output/bed_${st})
        do
            grep ${region} $output/bed_${st}/$lib > $output/${st}/bed_${region}/$lib &
        done
        wait

        # count mono
        eval $heatmap/count_rNMP.py $genome $output/${st}/bed_${region}/* -m -o $output/${st}/mono_${region}/test &
        wait
        rename s'/test_//' $output/${st}/mono_${region}/* -f

        # get chrom
        eval $heatmap/get_chrom.py $output/${st}/mono_${region}/* -o $output/${st}/${region}_mono.raw
        eval $heatmap/normalize.py $output/${st}/${region}_mono.raw $bg/${st}.tsv --name ${region} -o $output/${st}/${region}_mono.norm
        eval $heatmap/resort.py $output/${st}/${region}_mono.norm $order -c 2 -o $output/${st}/${region}_mono.tsv
        eval $heatmap/draw_heatmap.py $output/${st}/${region}_mono.tsv -b $bg/${st}.tsv --background_chrom ${region} --palette RdBu_r -o $output/plots/${st}/${region}

        # barplot
        eval $scripts/../heatmap_barplot/generate_bar_plot.py $output/${st}/${region}_mono.tsv -o $output/plots/${st}/${region}_barplot &
    done
done

# control region distribution
eval $scripts/replication_distribution.py $mito_count $order --same $output/same/bed_CR/* --oppo $output/oppo/bed_CR/* -o $output/plots/distribution/CR &

for aa in OH OL CR
do
    eval $scripts/compare_strands.py $output/same/${aa}_mono.raw $output/oppo/${aa}_mono.raw $order -o $output/plots/${aa}_strands --no_annot &
    eval $scripts/compare_strands.py $output/same/${aa}_mono.raw $output/oppo/${aa}_mono.raw $order -o $output/plots/${aa}_strands_annot &
done

eval $scripts/replication_distribution_combined.py $mito_count $order --same $output/same/bed_CR/* --oppo $output/oppo/bed_CR/* -o $output/plots/distribution/all_wt &
eval $scripts/replication_distribution_combined.py $mito_count $order --same $output/same/bed_CR/* --oppo $output/oppo/bed_CR/* --selected "HEK293T" "RNH2A-KO\ T3-8" "RNH2A-KO\ T3-17" --palette Set2 -o $output/plots/distribution/all_hek &
eval $scripts/replication_distribution_combined_vs_chipseq.py $mito_count $order $chipseq/*.bed --same $output/same/bed_CR/* --oppo $output/oppo/bed_CR/* -o $output/plots/distribution/with_chipseq &

wait