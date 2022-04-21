#!/bin/bash

for aa in intersect ref_bed raw info plots
do
    if ! [ -e $aa ]
    then
        mkdir $aa
    fi
done
 
for aa in intersect ref_bed
do
    for bb in cds noncoding random
    do
        if ! [ -e ${aa}/${bb} ]
        then
            mkdir ${aa}/${bb}
        fi
    done
done

for aa in heatmap barplot regplot
do
    if ! [ -e plots/$aa ]
    then
        mkdir plots/$aa
    fi
done

# scripts location
scripts="~/bio-storici/scripts/rNMP_hmt_analysis"

for aa in $(ls ../bed)
do
    for ty in cds noncoding random
    do
        for bb in $(ls ref_bed/${ty})
        do
            bedtools intersect -a ../bed/$aa -b ref_bed/${ty}/$bb -nonamecheck -s > intersect/${ty}/${aa}_${bb}_nt &
            bedtools intersect -a ../bed/$aa -b ref_bed/${ty}/$bb -nonamecheck -S > intersect/${ty}/${aa}_${bb}_t &
        done
    done
done
wait

rename 's/.bed_nt/_nontemplate.bed/;s/.bed_t/_template.bed/;s/.bed_/_/' intersect/*/*

# count intersections
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        wc -l intersect/${ty}/*_${st}.bed | head -n -1 | sed 's/^ *//;s/ /\t/;s/_${st}.bed//;s/_/\t/g;s/inter.*FS/FS/;s/inter.*RD/RD/' | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > raw/${ty}_${st}.raw &
    done
done
wait

# add cds information
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        eval $scripts/gene_analysis/add_info.py raw/${ty}_${st}.raw ../order_human.tsv $scripts/gene_analysis/ref/${ty}.gtf ../mito_count.tsv -c 2 -o info/${ty}_${st}.tsv &   
    done
done
wait

# draw heatmap
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        eval $scripts/gene_analysis/add_info.py raw/${ty}_${st}.raw ../order_human.tsv $scripts/gene_analysis/ref/${ty}.gtf ../mito_count.tsv -c 2 -o info/${ty}_${st}.tsv &   
    done
done
wait

# draw heatmap
for st in template nontemplate
do
    for ty in cds noncoding random
    do
        eval $scripts/gene_analysis/draw_heatmap.py info/${ty}_${st}.tsv --no_cbar -o plots/heatmap/${ty}_${st} &   
    done
done
wait

# generate barplot
for st in template nontemplate
do
    cat info/cds_${st}.tsv > info/${st}.tsv
    tail -n +2 info/noncoding_${st}.tsv >> info/${st}.tsv
    sed -i 's/MT-//g' info/${st}.tsv
    eval $scripts/gene_analysis/draw_barplot.py info/${st}.tsv $scripts/gene_analysis/gene_names.tsv -o plots/barplot/${st} &
done
wait

# generate regression plot
for st in template nontemplate
do
    for ty in cds
    do
        eval $scripts/gene_analysis/draw_regplot.py info/${ty}_${st}.tsv -o plots/regplot/${ty}_${st} &
    done
done
wait