#!/bin/bash
set -x 
#########################################################################
# File Name: 3-summarize_res.sh
# Created on: 2021-03-28 15:27:39
# Author: JFF
# Last Modified: 2022-03-30 15:35:31
# Description: 去掉低质量的imputeSNP，筛maf，geno，hwe等。合并结果，转成plink格式，再转成vcf格式
# Usage: 
# Input: 
# Output: 
#########################################################################

folder="OA_217.TOP.processed.chr1_23.QC.split"

sample=`echo ${folder}|sed 's/split$/sample/g'`
if [ -f "${folder}.merge_list" ];then
    rm ${folder}.merge_list
fi

ls ${folder}/*.imputed.gz |while read imputed
do
    #modify chr
    chr=`basename $imputed|awk -F'_' '{print $1}'|awk -F'-' '{print $1}'`
    gzip -d -c $imputed |sed 's/^---/'$chr'/g' > temp.$imputed
    mv temp.$imputed $imputed
    #to bed
    info=`echo $imputed|sed 's/\.gz$/_info/g'`
    prefix=`echo $imputed|sed 's/\.gz$//g'`
    plink --gen $imputed --sample $sample --qual-scores $info 7 2 1 --qual-threshold 0.3  --hard-call-threshold 0.3 --make-bed --out $prefix && echo ${prefix}.bed ${prefix}.bim ${prefix}.fam >> ${folder}.merge_list
done

#merge & to vcf & filter geno + mind + maf + hwe
plink --merge-list ${folder}.merge_list --geno 0.05 --make-bed --out ${folder}.merge.geno_0.05
plink --bfile ${folder}.merge.geno_0.05 --mind 0.05 --make-bed --out ${folder}.merge.geno_0.05.ind_0.05
plink --bfile ${folder}.merge.geno_0.05.ind_0.05 --maf 0.05 --make-bed --out ${folder}.merge.geno_0.05.ind_0.05.maf_0.05
plink --bfile ${folder}.merge.geno_0.05.ind_0.05.maf_0.05 --hwe 1e-6 midp --make-bed --out ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6

#update rs id
awk '{if($2~/^rs/){split($2,a,":");print $2"\t"a[1]}else{print $2"\t"$2}}' ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.bim > ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.bim.updatename
plink --bfile ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6 --update-name ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.bim.updatename --make-bed --out ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.updatename

#add chr
sed 's/^/chr/g' ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.updatename.bim -i
#to vcf，23号染色体用 --split-par 处理,将男性X染色体作为二倍体处理
#plink2 --bfile ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.updatename  --recode vcf-4.2 id-paste=iid --split-par 0 1 --output-chr chr26 --out ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.hwe_1e-6.updatename

#remove temp file
#rm ${folder}.merge.geno_0.05.bed ${folder}.merge.geno_0.05.bim ${folder}.merge.geno_0.05.fam ${folder}.merge.geno_0.05.log
#rm ${folder}.merge.geno_0.05.ind_0.05.bed ${folder}.merge.geno_0.05.ind_0.05.bim ${folder}.merge.geno_0.05.ind_0.05.fam ${folder}.merge.geno_0.05.ind_0.05.log
#rm ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.bed ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.bim ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.fam ${folder}.merge.geno_0.05.ind_0.05.maf_0.05.log
