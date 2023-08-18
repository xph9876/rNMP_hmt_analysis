#!/bin/bash

# usage gene_analysis.sh <gene_analysis folder> <heatmap script> <output folder> <bed folder> <order> <mito_count>
scripts=$1
heatmap=$2
output=$3
bed_folder=$4
order=$5
mito_count=$6
ref_bed="$1/ref_bed"
ref="$1/ref"
control_bed="$1/control_bed"
genome="$1/../refseq/hg38_chrM.fa"

for aa in intersect raw info plots cg
do
    if ! [ -d $output/$aa ]
    then
        mkdir $output/$aa
    fi
done
 
for aa in intersect
do
    for bb in cds noncoding random
    do
        if ! [ -d $output/${aa}/${bb} ]
        then
            mkdir $output/${aa}/${bb}
        fi
    done
done

for aa in heatmap barplot regplot regplot_cg
do
    if ! [ -d $output/plots/$aa ]
    then
        mkdir $output/plots/$aa
    fi
done

cp $control_bed/RD*.bed $bed_folder
for aa in $(ls $bed_folder)
do
    for ty in cds noncoding random
    do
        for bb in $(ls $ref_bed/${ty})
        do
            bedtools intersect -a $bed_folder/$aa -b $ref_bed/${ty}/$bb -nonamecheck -s > $output/intersect/${ty}/${aa}_${bb}_nt &
            bedtools intersect -a $bed_folder/$aa -b $ref_bed/${ty}/$bb -nonamecheck -S > $output/intersect/${ty}/${aa}_${bb}_t &
        done
    done
done
wait

rename 's/.bed_nt/_nontemplate.bed/;s/.bed_t/_template.bed/;s/.bed_/_/' $output/intersect/*/* -f

# count intersections
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        wc -l $output/intersect/${ty}/*_${st}.bed | head -n -1 | sed "s/^ *//;s/ /\t/;s/\t.*\//\t/;s/_${st}.bed//;s/_/\t/g" | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > $output/raw/${ty}_${st}.raw &
    done
done
wait

# add cds information
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        eval $scripts/add_info.py $output/raw/${ty}_${st}.raw $order $scripts/ref/${ty}.gtf $mito_count -c 2 -o $output/info/${ty}_${st}.tsv &   
    done
done
wait

# draw heatmap
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        eval $scripts/draw_heatmap.py $output/info/${ty}_${st}.tsv --no_cbar -o $output/plots/heatmap/${ty}_${st} &   
    done
done
wait

# generate barplot
for st in template nontemplate
do
    cat $output/info/cds_${st}.tsv > $output/info/${st}.tsv
    tail -n +2 $output/info/noncoding_${st}.tsv >> $output/info/${st}.tsv
    sed -i 's/MT-//g' $output/info/${st}.tsv
    eval $scripts/draw_barplot.py $output/info/${st}.tsv $scripts/gene_names.tsv -o $output/plots/barplot/${st} &
done
wait

# generate regression plot
for st in template nontemplate
do
    for ty in cds
    do
        eval $scripts/draw_regplot.py $output/info/${ty}_${st}.tsv -o $output/plots/regplot/${ty}_${st} &
    done
done
wait

# remove random bed files
rm $bed_folder/RD1.bed $bed_folder/RD2.bed $bed_folder/RD3.bed -f

# Correlation between CG count to PPB
bedtools getfasta -fi $genome -bed $ref/hg38_mt_cds_ensembl.bed -fo $output/cg/cds.fa -s
eval $heatmap/count_background.py $output/cg/cds.fa --mono -s -o $output/cg/cds_mono.tsv 
eval $scripts/append_mono.py $output/info/cds_nontemplate.tsv $output/cg/cds_mono.tsv -o $output/cg/cds_nontemplate.tsv
eval $scripts/draw_regplot_cg.py $output/cg/cds_nontemplate.tsv -o $output/plots/regplot_cg/cds_nontemplate &


