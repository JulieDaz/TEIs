#-*- coding: utf-8 -*-
### script to select the TEs
### change the chr name to match my other files
### change the info col to put in output: keep name, identity and add the TEs class

import argparse
import re

from collections import defaultdict


def parse_renamed_MASiVE(renamed_file) :

    dict_chrom_name = {}
    with open(renamed_file, "r") as renamed : 

        for line in renamed : 

            if not line.startswith("#") :

                myline = line.split("\t")
                orig_name = myline[1].split(" ")[0]
                
                if myline[5] != "" :
                    masive_name = myline[-1].strip()
                    dict_chrom_name[orig_name] = masive_name
                    species = masive_name.split("_")[0]
                    # print(species)

    print(dict_chrom_name)
    return dict_chrom_name, species


def select_TEs_family(TEs_gff_file) : 

    list_TEs_fam = ["CACTA_TIR_transposon", "Copia_LTR_retrotransposon", 
                    "Gypsy_LTR_retrotransposon", "L1_LINE_retrotransposon", 
                    "LINE_element", "LTR_retrotransposon", "Mutator_TIR_transposon", 
                    "PIF_Harbinger_TIR_transposon", "RTE_LINE_retrotransposon", 
                    "Tc1_Mariner_TIR_transposon", "hAT_TIR_transposon", "helitron"]

    dict_TEs = defaultdict(list)

    with open(TEs_gff_file, "r") as TEs_gff : #, open(f"{TEs_gff_file}.TEonly.gff3", "w") as output : 

        for line in TEs_gff : 

            if not line.startswith("#") :

                myline = line.strip().split("\t")

                if len(myline[0].split("_")) > 1 :
                # if species == "Zmays-Oh7B" :
                    chrom = myline[0].split("_")[1]
                else :
                    chrom = myline[0]

                fam_TE = myline[2]

                start = myline[3]
                end = myline[4]
                strand = myline[6]

                if fam_TE in list_TEs_fam :
                
                    info = myline[8]
                    list_info = info.split(";")
                    # print(list_info)
                    for ele in list_info :
                        if "ID=" in ele : 
                            id_TE = ele                        
                        if "Name=" in ele : 
                            name_TE = ele

                    new_info = f"{fam_TE};{id_TE};{name_TE}"
                    dict_TEs.setdefault(chrom, []).append([start, end, new_info, ".", strand])

                    # output.write(line)

    return dict_TEs

def update_chrom_name(dict_chrom_name, dict_TEs, species, type_TEs) :

    dict_TEs_newname = {}

    for orig_chrom, masive_chrom in dict_chrom_name.items() : 
            if orig_chrom in dict_TEs :
                dict_TEs_newname[masive_chrom] = dict_TEs[orig_chrom]

    with open(f"{species}_MaizeGDB_{type_TEs}TEs.bed", "w") as output : 
        for chrom, genes_coord in dict_TEs_newname.items() :

            for gene in genes_coord :
                gene_to_write = "\t".join(str(x) for x in gene)
                # print(f"{chrom}\t{gene_to_write}")
                output.write(f"{chrom}\t{gene_to_write}\n")



###### MAIN

parser = argparse.ArgumentParser(description="")
parser.add_argument("-ren", help="Zm-$1-REFERENCE-NAM-1.0.fa.renamed")
parser.add_argument("-gff", help="Maize TE gff: Zm-$1-REFERENCE-NAM-1.0.TE.gff3")
parser.add_argument("-type", help="intact or all")
# parser.add_argument("-pan", help="Pangenes file")
args = parser.parse_args()

renamed_file = args.ren
gff_TEs_f = args.gff
type_TEs = args.type
# pangene_f = args.pan
# output_f = bed_f.split(".")[0]

dico_chrom_name, species = parse_renamed_MASiVE(renamed_file)
dico_TEs = select_TEs_family(gff_TEs_f)
dico_genes_for_bins = update_chrom_name(dico_chrom_name, dico_TEs, species, type_TEs)


