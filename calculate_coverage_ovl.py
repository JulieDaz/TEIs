#-*- coding: utf-8 -*-

import argparse
from collections import defaultdict

import pandas as pd


def ovl_bedtools_intersect(input_f, output_f) :

    with open(input_f, "r") as res, open(f"{output_f}_ovl_coverage.txt", "w") as f :
        f.write(f"chr\tstart_TE\tend_TE\tTE_id\tstart_intron\tend_intron\tintron_id\tovl_len\tTE_len\tTE_cov\tintron_len\tintron_cov\n")

        for line in res : 

            myline = line.strip().split("\t")
            # print(myline)
            chrom = myline[0]
            start_A = int(myline[1])
            end_A = int(myline[2])
            id_A = myline[3]
            start_B = int(myline[5])
            end_B = int(myline[6])
            id_B = myline[7]
            ovl = int(myline[8])

            if start_B != -1 :
                len_A = end_A - start_A
                len_B = end_B - start_B
                # coverage_A = round((ovl/len_A)*100, 2)
                # coverage_B = round((ovl/len_B)*100, 2)
                coverage_A = (ovl/len_A)*100  
                coverage_B = (ovl/len_B)*100  
                # print(f"{chrom}\t{start_A}\t{end_A}\t{id_A}\t{start_B}\t{end_B}\t{id_B}\t{ovl}\t{coverage_A}\t{coverage_B}\n")
                f.write(f"{chrom}\t{start_A}\t{end_A}\t{id_A}\t{start_B}\t{end_B}\t{id_B}\t{ovl}\t{len_A}\t{coverage_A}\t{len_B}\t{coverage_B}\n")


###### MAIN 
parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", help="")
args = parser.parse_args()

input_f = args.i
output_f = input_f.split(".")[0]

ovl_bedtools_intersect(input_f,output_f)
