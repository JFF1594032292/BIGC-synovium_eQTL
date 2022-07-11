# -*- coding: utf-8 -*-
#########################################################################
# File Name: 2-impute.py
# Created on : 2021-03-23 11:13:46
# Author: JFF
# Last Modified: 2021-03-23 11:13:46
# Description: 对分割后的文件进行impute
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
t0 = time.time()

folder = "OA_217.TOP.processed.chr1_23.QC.split"
sample_file = re.sub("\.split$", ".sample", folder)
ref_folder = "/home/jiangfeng/data/impute2_ref_1000GP_Phase3"

file_list = os.listdir(folder)
suffix_set = set([i.split(".")[0] for i in file_list if re.search("\.gen$", i)])


def task(file_suffix):
    t1 = time.time()
    start = file_suffix.split("_")[1]
    end = file_suffix.split("_")[2]
    ch = file_suffix.split("_")[0]
    if "23-" not in ch:
        d = {}
        d["m"] = ref_folder+"/"+"genetic_map_chr%s_combined_b37.txt" % ch
        d["h"] = ref_folder+"/"+"1000GP_Phase3_chr%s.hap.gz" % ch
        d["l"] = ref_folder+"/"+"1000GP_Phase3_chr%s.legend.gz" % ch
        d["g"] = folder+"/"+file_suffix+".gen"
        d["strand_g"] = folder+"/"+file_suffix+".strand"
        d["int"] = "%s %s" % (start, end)
        d["o"] = folder+"/"+file_suffix+".imputed"
        os.system("impute2 -m %(m)s -h %(h)s -l %(l)s -g %(g)s -strand_g %(strand_g)s -int %(int)s -Ne 20000 -seed 123456 -k 100 -k_hap 505 -buffer 300 -no_remove -o %(o)s" % d)
    elif "23-" in ch:
        ch = "X_"+ch.split("-")[1]
        d = {}
        d["m"] = ref_folder+"/"+"genetic_map_chr%s_combined_b37.txt" % ch
        d["h"] = ref_folder+"/"+"1000GP_Phase3_chr%s.hap.gz" % ch
        d["l"] = ref_folder+"/"+"1000GP_Phase3_chr%s.legend.gz" % ch
        d["g"] = folder+"/"+file_suffix+".gen"
        d["sample_g"] = sample_file
        d["strand_g"] = folder+"/"+file_suffix+".strand"
        d["int"] = "%s %s" % (start, end)
        d["o"] = folder+"/"+file_suffix+".imputed"
        if "PAR1" in ch or "PAR2" in ch:
            os.system("impute2 -chrX -Xpar -m %(m)s -h %(h)s -l %(l)s -g %(g)s -sample_g %(sample_g)s -strand_g %(strand_g)s -int %(int)s -Ne 20000 -seed 123456 -k 100 -k_hap 505 -buffer 300 -no_remove -o_gz -o %(o)s" % d)
        elif "nonPAR" in ch:
            os.system("impute2 -chrX -m %(m)s -h %(h)s -l %(l)s -g %(g)s -sample_g %(sample_g)s -strand_g %(strand_g)s -int %(int)s -Ne 20000 -seed 123456 -k 100 -k_hap 505 -buffer 300 -no_remove -o_gz -o %(o)s" % d)
    print(file_suffix, "is done. %.5g" % (time.time()-t1))


pool = Pool(28)
pool.map(task, suffix_set)
pool.close()
pool.join()
print("ALL DONE. %.5g" % (time.time()-t0))
