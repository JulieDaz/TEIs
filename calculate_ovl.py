#-*- coding: utf-8 -*-

import argparse
from collections import defaultdict

import pandas as pd


def ovl_bedtools_intersect(intersect_f, type_TE) :

    num_lines = sum(1 for line in open(intersect_f))

    dict_ovl_stats = {"ovl_all":0,
                    f"{type_TE}inGene":0,
                    f"Genein{type_TE}":0,
                    f"ovl_end{type_TE}":0,
                    f"ovl_start{type_TE}":0,
                    "partial_ovl":0,
                    f"No_{type_TE}-gene_ovl":0}
    list_ovl = []
    with open(intersect_f, "r") as res : 

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
            cov = myline[8]

            if start_B != -1 :

                dict_ovl_stats["ovl_all"] += 1

                # if start_A > start_B and end_A < end_B : 
                #     # check if id_A already in one of sublist
                #     # if not any(id_A in sublist for sublist in list_ovl) :
                #     list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, "{type_TE}inGene"])
                #     dict_ovl_stats["{type_TE}inGene"] += 1
                #     # print(chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, start_A > start_B, end_A < end_B)

                # elif start_A < start_B and end_A > end_B :
                #     # if not any(id_A in sublist for sublist in list_ovl) :
                #     list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, "Genein{type_TE}"])
                #     dict_ovl_stats["Genein{type_TE}"] += 1
                #     # print(chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, start_A < start_B, end_A > end_B)


                # elif start_A < start_B and end_A < end_B :
                #     # if not any(id_A in sublist for sublist in list_ovl) :
                #     list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, "ovl_end{type_TE}"])
                #     dict_ovl_stats["ovl_end{type_TE}"] += 1

                # elif start_A > start_B and end_A > end_B :
                #     # if not any(id_A in sublist for sublist in list_ovl) :
                #     list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, "ovl_start{type_TE}"])
                #     dict_ovl_stats["ovl_start{type_TE}"] += 1

                if start_A > start_B and end_A < end_B : 
                    # check if id_A already in one of sublist
                    dict_ovl_stats[f"{type_TE}inGene"] += 1
                    # if not any(id_A in sublist for sublist in list_ovl) :
                    list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, f"{type_TE}inGene"])
                    # print(chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, start_A > start_B, end_A < end_B)

                elif start_A < start_B and end_A > end_B :
                    dict_ovl_stats[f"Genein{type_TE}"] += 1
                    # if not any(id_A in sublist for sublist in list_ovl) :
                    list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, f"Genein{type_TE}"])
                    # print(chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, start_A < start_B, end_A > end_B)


                elif start_A < start_B and end_A < end_B :
                    dict_ovl_stats[f"ovl_end{type_TE}"] += 1
                    dict_ovl_stats["partial_ovl"] += 1
                    # if not any(id_A in sublist for sublist in list_ovl) :
                    list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, f"ovl_end{type_TE}"])

                elif start_A > start_B and end_A > end_B :
                    dict_ovl_stats[f"ovl_start{type_TE}"] += 1
                    dict_ovl_stats["partial_ovl"] += 1
                    # if not any(id_A in sublist for sublist in list_ovl) :
                    list_ovl.append([chrom, id_A, start_A, end_A, id_B, start_B, end_B, cov, f"ovl_start{type_TE}"])

            else :
                dict_ovl_stats[f"No_{type_TE}-gene_ovl"] += 1

    # print(dict_ovl_stats)
    return list_ovl, dict_ovl_stats, num_lines


def stats(dict_ovl_stats, output_file, type_TE, nb_TE, nb_gene, num_lines) :

    dict_ovl_perc = defaultdict(list)
    dict_ovl_perc[type_TE] = [nb_TE]
    dict_ovl_perc["genes"] = [nb_gene]
    dict_ovl_perc["intersect"] = [num_lines]

    for cat, value in dict_ovl_stats.items() :

        perc_ovl = round(value/num_lines, 3)
        perc_gene = round(value/nb_gene, 3)
        dict_ovl_perc[cat] = [value]
        dict_ovl_perc[f"{cat}_perc_{type_TE}"] = [perc_ovl]
        dict_ovl_perc[f"{cat}_perc_gene"] = [perc_gene]

    # print(dict_ovl_perc)
    dict_pd = pd.DataFrame.from_dict(dict_ovl_perc, orient="columns", )
    dict_pd.to_csv(f"{output_file}_ovl_stats.txt", sep="\t",index=False)


def write_output(list_ovl, output_file, type_TE) :

    list_written_in_bed = []

    with open(f"{output_file}_ovl.txt", "w") as ovl, open(f"{output_file}_ovl.bed", "w") as bed : 

        ovl.write(f"chr\t{type_TE}_id\tstart_{type_TE}\tend_{type_TE}\tgene_id\tstart_gene\tend_gene\tovl_len\tcategory\n")

        for sublist in list_ovl :

            sublist_to_write = "\t".join(str(ele) for ele in sublist)
            ovl.write(f"{sublist_to_write}\n")

            if sublist[-1] == f"{type_TE}inGene" :

                if sublist[1] not in list_written_in_bed :

                    sublist_TE = f"{sublist[0]}\t{sublist[2]}\t{sublist[3]}\t{sublist[1]}"
                    bed.write(f"{sublist_TE}\n")
                    list_written_in_bed.append(sublist[1])


###### MAIN 
parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", help="")
parser.add_argument("-type", help="SV or TE")
parser.add_argument("-transposon", help="Number of SVs in the species")
parser.add_argument("-gene", help="Number of genes in the species")
args = parser.parse_args()

intersect_f = args.i
type_TE = args.type
nb_transposon = int(args.transposon)
nb_gene = int(args.gene)

# intersect_f = "Atrun_SV_introns--intersect.txt"
# output_f = "Atrun_introns_ovl.txt"

output_f = ".".join(intersect_f.split(".")[:-1])
# print(output_f)
# print(nb_transposon, nb_gene)
list_ovl, dico_ovl_stats, nb_lines = ovl_bedtools_intersect(intersect_f, type_TE)
stats(dico_ovl_stats, output_f, type_TE, nb_transposon, nb_gene, nb_lines)
write_output(list_ovl, output_f, type_TE)

