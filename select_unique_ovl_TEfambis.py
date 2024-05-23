#-*- coding: utf-8 -*-

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_bedtools_intersect_cov(input_file) :

    '''
    Returns a dictionary with the key = TE, and the value is a dict with the key
    = intron_id and the value is the sum of coverages (cov intron + cov TE)
    Also returns a list which has all the data of the bedtools intersect with cov
    '''

    dico_TE_introns_highcov = defaultdict(lambda: defaultdict(float))
    list_file = []
    with open(input_file, "r") as cov_file : 

        for line in cov_file :
            if not line.startswith("chr") :

                myline = line.strip().split("\t")
                list_file.append(myline)
                # print(myline)
                ident_te = myline[3] 
                ident_intron = myline[6]
                cov_te = float(myline[9])
                cov_intron = float(myline[11])

                if cov_te >= 99 and cov_intron >= 99 :

                    sum_cov = cov_te + cov_intron
                    # print(ident_te, ident_intron)
                    dico_TE_introns_highcov[ident_te][ident_intron] = sum_cov

    # print(dico_TE_introns_highcov)
    print(len(dico_TE_introns_highcov))
    # exit()
    return dico_TE_introns_highcov, list_file


def check_intron_present_multiple_time(list_intron_selected, intron, cov) : 

    # loop through sublists
    for sublist in list_intron_selected:
        # if the intron already exists
        if sublist[1] == intron:
            # print(sublist[1], intron)
            # check if the coverage of the new intron is higher
            if cov > sublist[2]:
                # if it's higher we want to add the new intron and remove the old one
                # so I return the whole sublist to easily remove it
                return sublist
            else:
                # if it's not we don't want to add it
                # we don't do anything in select_transcript, already have the best TE/intron
                print(sublist) 
                print(intron, cov)
                return False

    # we couldn't find the new intron in the list of sublist so we want to add it
    # without removing anything
    return True


def select_transcript(dico_TEs_introns_highcov) :

    list_selected_introns = []
    for TEs, introns in dico_TEs_introns_highcov.items() :

        # print(TEs, introns)

        if len(introns) > 1 :
            print("caca")
            ## if several introns for a TE (usually several transcripts of the same gene)
            ## then take the one with highest sum, if identical it takes
            ## the first one or a random one can't remember
            max_intron = max(introns, key=introns.get)
            print(max_intron)
            cov = dico_TEs_introns_highcov[TEs][max_intron]
            
            new_intron = check_intron_present_multiple_time(list_selected_introns, max_intron, cov)
            print(new_intron)
            # first time we encounter this intron
            if new_intron is True:
                # keeping TEs and intron as 1 TE can be overlapped by several TEs
                # and all these TEs will be selected in the next function if the 
                # list only contains intron_id
                # eg: intron1 ovl TE1 and TE2, but TE1 only ovl on 20% and TE2 98%
                # we don't want TE1 in the final output
                list_selected_introns.append([TEs, max_intron, cov])

            # new_intron contains the previous TE/intron that we remove
            # bc we found a better one, so is not "True"
            elif new_intron:
                intron_to_rm = new_intron
                list_selected_introns.append([TEs, max_intron, cov])
                list_selected_introns.remove(intron_to_rm)

        else :
            for intron, value in introns.items() :
                new_intron = check_intron_present_multiple_time(list_selected_introns, intron, value)

                if new_intron is True:
                    list_selected_introns.append([TEs, intron, value])
                elif new_intron:
                    print(new_intron)
                    intron_to_rm = new_intron
                    list_selected_introns.append([TEs, intron, value])
                    list_selected_introns.remove(intron_to_rm)

                    # print("kiki", new_intron, intron, value)

    print(len(list_selected_introns))
    return list_selected_introns


def count_ovl_category(list_file, list_selected_introns) :

    dict_count_cov = defaultdict(lambda: defaultdict(int))

    for sublist in list_file :
        for te, intron, cov in list_selected_introns :

            if te in sublist and intron in sublist :
                # print(sublist)
                # dict_file_updated
                sublist.append("Y")
                # print(sublist)
                cov_te = float(sublist[9])
                cov_intron = float(sublist[11])
                sublist[9] = round(cov_te, 2)
                sublist[11] = round(cov_intron, 2)

                te_type = te.split(";")[0]
                if cov_te >= 95 and cov_intron >= 95 :
                    dict_count_cov[te_type]["95%"] += 1
                if cov_te >= 96 and cov_intron >= 96 :
                    dict_count_cov[te_type]["96%"] += 1
                if cov_te >= 97 and cov_intron >= 97 :
                    dict_count_cov[te_type]["97%"] += 1
                if cov_te >= 98 and cov_intron >= 98 :
                    dict_count_cov[te_type]["98%"] += 1
                if cov_te >= 99 and cov_intron >= 99 :
                    dict_count_cov[te_type]["99%"] += 1

        if sublist[-1] != "Y" :
            sublist.append("N")
            sublist[9] = round(float(sublist[9]), 2)
            sublist[11] = round(float(sublist[11]), 2)

    # print(dict_count_cov)
    return list_file, dict_count_cov


def write_outputs(list_file_updated, dict_count_cov, output_f) :

    # check if list_file_updated is not empty to avoid writing empty files
    if not list_file_updated == [] : 
        species = output_f.split("_")[0]
        with open(f"{output_f}_selected.TEST.txt", "w") as f: 
            
            f.write(f"chr\tstart_TE\tend_TE\tTE_id\tstart_intron\tend_intron\tintron_id\tovl_len\tTE_len\tTE_cov\tintron_len\tintron_cov\tselected\n")

            for sublist in list_file_updated :

                sublist_to_write = "\t".join(str(x) for x in sublist)
                f.write(f"{sublist_to_write}\n")

        dict_pd = pd.DataFrame.from_dict(dict_count_cov, orient = "index")
        dict_pd['TEfam'] = dict_pd.index

        # insert species code in 1st position
        dict_pd.insert(0, "species", species)
        dict_pd.to_csv(f"{output_f}_counts.TEST.txt", sep="\t",index=False)

    else :
        print("No data to write.")

#### MAIN 

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", help="")
args = parser.parse_args()

input_f = args.i
# input_f = "Zmays-B73-v5_TEs_introns--intersect_ovl_coverage.txt"
output_f = Path(input_f).stem
print(output_f)

dict_TE_introns_highcov, list_file = parse_bedtools_intersect_cov(input_f)
list_selected_introns = select_transcript(dict_TE_introns_highcov)
list_file_updated, dict_count_cov = count_ovl_category(list_file, list_selected_introns)
write_outputs(list_file_updated, dict_count_cov, output_f)

