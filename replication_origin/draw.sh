#!/bin/bash

if ! [ -e plots ]
then
    mkdir plots
fi

if ! [ -e plots/distribution ]
then
    mkdir plots/distribution
fi

for st in same oppo
do
    if ! [ -e ${st} ]
    then
        mkdir ${st}
    fi

    if ! [ -e plots/${st} ]
    then
        mkdir plots/${st}
    fi

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
        ../../distributed_scripts/RibosePreferenceAnalysis/count_rNMP.py ../../../genome/saccer.fa ${st}/bed_${region}/* -m -o ${st}/mono_${region}/test &
        wait
        rename test_ '' ${st}/mono_${region}/*

        # get chrom
        ../../distributed_scripts/RibosePreferenceAnalysis/get_chrom.py ${st}/mono_${region}/* -o ${st}/${region}_mono.raw
        ../../distributed_scripts/RibosePreferenceAnalysis/normalize.py ${st}/${region}_mono.raw bg/${st}/mito_mono.raw --name ${region} -o ${st}/${region}_mono.norm
        ../../distributed_scripts/RibosePreferenceAnalysis/resort.py ${st}/${region}_mono.norm ../order_human.tsv -c 2 -o ${st}/${region}_mono.tsv
        ../../distributed_scripts/RibosePreferenceAnalysis/draw_heatmap.py ${st}/${region}_mono.tsv -b bg/${st}/mito_mono.raw --background_chrom ${region} --palette RdBu_r -o plots/${st}/${region}

        # barplot
        ../generate_bar_plot.py ${st}/${region}_mono.tsv -o plots/${st}/${region}_barplot &
    done
done

# control region distribution
./replication_distribution.py ../mito_count.tsv ../order_human.tsv --same same/bed_CR/* --oppo oppo/bed_CR/* -o plots/distribution/CR

for aa in OH OL CR
do
    ./compare_strands.py same/${aa}_mono.raw oppo/${aa}_mono.raw ../order_human.tsv -o plots/${aa}_strands
    
done