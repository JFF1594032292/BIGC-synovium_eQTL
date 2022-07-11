# -*- coding: utf-8 -*-
#########################################################################
# File Name: 1-split_plink_bed2gen.py
# Created on : 2021-03-22 21:25:19
# Author: JFF
# Last Modified: 2021-03-22 21:25:19
# Description: 将plink bed格式按染色体及5mb 分割生成成 .gen 文件，每条染色体整个按5m分割，防止遗漏SNP。并生成 .strand_g文件(全都按正链)
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import time
import math
import random
from multiprocessing import Pool
geno_file_prefix = "OA_217.TOP.processed.chr1_23.QC"
folder = geno_file_prefix+".split"
int_length = 1000000*5  # 分割长度

if not os.path.exists(folder):
    os.mkdir(folder)

# X染色体分割范围
#ref_folder = "/home/jiangfeng/data/impute2_ref_1000GP_Phase3"  # X染色体
#X_file_list = [ref_folder+"/"+i for i in os.listdir(ref_folder) if re.search("^genetic_map_chrX_.+combined_b37\.txt$", i)]
#d_X_pos = {re.sub("^genetic_map_chrX_|_combined_b37\.txt$", "", os.path.basename(i)): [int(j.strip().split(" ")[0]) for j in open(i) if j.strip().split(" ")[0] != "position"] for i in X_file_list}
#d_X_range = {i: [min(d_X_pos[i]), max(d_X_pos[i])] for i in d_X_pos}  # ["PAR1":[150118,2695340],"nonPAR":[2703391,154929412],"PAR2":[154969038,155235078]]
#直接指定范围
d_X_range= {"PAR1":[1,2699520],"nonPAR":[2699521,154931044],"PAR2":[154931045,155270560]} #hg19 X染色体区域划分，来自 https://www.cog-genomics.org/plink/2.0/data#split_par

# plink get gen
os.system("plink --bfile %s --recode oxford --out %s" % (geno_file_prefix, geno_file_prefix))

# split chr & pos
with open(geno_file_prefix+".gen") as f1:
    d_chr = {}
    for i in f1:
        j = i.strip().split(" ")
        d_chr.setdefault(j[0], {})[int(j[2])] = i  # chr:{pos:line,}
# chromsome size
d_chrsize={}
with open("/home/jiangfeng/data/Genome_ano/chromesize-hg19.filtered") as f1:
    for i in f1:
        i=i.strip().split("\t")
        if i[0]=="chrX":
            i[0]="chr23"
        i[0]=re.sub("^chr","",i[0])
        d_chrsize[i[0]]=int(i[1]) #chr23:155270560
#define task
buffer_length=300000
def task_auto(ch): #常染色体
    chrsize=d_chrsize[ch]
    segment_list=[[i,i+int_length-1] for i in range(1,chrsize,int_length)]
    for segment in segment_list:
        #分割时也要留出buffer 区域的SNP
        start, end = segment
        start_expand = start - buffer_length*2
        end_expand = end + buffer_length*2
        #输出文件
        file_name1 = folder+"/"+"_".join([ch, str(start), str(end)])+".gen"
        file_name2 = folder+"/"+"_".join([ch, str(start), str(end)])+".strand"
        ff1 = open(file_name1, 'w')  # gen 文件
        ff2 = open(file_name2, 'w')  # strand 文件
        for pos in d_chr[ch]:
            if start_expand <= pos <= end_expand:
                ff1.write(d_chr[ch][pos])
                ff2.write(d_chr[ch][pos].strip().split(" ")[2]+"\t+\n")
        ff1.close()
        ff2.close()
    print(ch, "done")

def task_X(ch):  # X染色体
    for region_name in d_X_range:
        region_section = d_X_range[region_name]  # [2703391,154929412]
        segment_list=[[i, min(region_section[1], i+int_length)] for i in range(region_section[0], region_section[1]+1, int_length+1)]
        for segment in segment_list:
            #分割时也要留出buffer 区域的SNP
            start, end = segment
            start_expand = start - buffer_length*2
            end_expand = end + buffer_length*2
            #输出文件
            file_name1 = folder+"/"+"_".join([ch+"-"+region_name, str(start), str(end)])+".gen"
            file_name2 = folder+"/"+"_".join([ch+"-"+region_name, str(start), str(end)])+".strand"
            ff1 = open(file_name1, 'w')  # gen 文件
            ff2 = open(file_name2, 'w')  # strand 文件
            for pos in d_chr[ch]:
                if start_expand <= pos <= end_expand:
                    ff1.write(d_chr[ch][pos])
                    ff2.write(d_chr[ch][pos].strip().split(" ")[2]+"\t+\n")
            ff1.close()
            ff2.close()
        print(ch, "done")
    print(ch, "done")

# start
pool=Pool(23)
for ch in range(1, 23):
    pool.apply_async(task_auto,(str(ch),))
pool.apply_async(task_X,("23",))
pool.close()
pool.join()
