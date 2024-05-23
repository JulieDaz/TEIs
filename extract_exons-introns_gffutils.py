#-*- coding: utf-8 -*-

import argparse
import os
from collections import OrderedDict, defaultdict

import gffutils


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

    return dict_chrom_name, species


def transform_func(f):

    if f.featuretype == 'gene':
        return f
    elif "RNA" in f.featuretype:
        f.attributes['ID'][0] += '_transcript'
    else:
        # assume that anything else is a child of a transcript, so we need
        # to edit the "Parent" attribute
        if 'Parent' in f.attributes:
            f.attributes['Parent'] = [i + '_transcript' for i in f.attributes['Parent']]

    return f


def get_exons_introns(gff) :

    dict_exons = defaultdict(list)
    dict_introns = defaultdict(list)

    # try:
        # create db with the given gff
        # either give a name to the db
        # or :memory: so it uses the RAM
    db = gffutils.create_db(gff, ":memory:", verbose=True, id_spec="ID", transform=transform_func, merge_strategy="create_unique")

    for mrna in db.features_of_type(featuretype="mRNA") :
        # pour especes avec gene id = mrna id
        try :
            # on regarde si le mRNA a un parent
            db[mrna.attributes["Parent"][0]]
        except :
            # si pas de parent
            db = db.add_relation(parent=db[mrna.id.strip("_transcript")], child=mrna, level=1)
            mrna_children = db.children(mrna.id)

            for child in mrna_children:
                db = db.add_relation(parent=db[mrna.id.strip("_transcript")], child=child, level=2)
    # except: 
    #     db = gffutils.interface.FeatureDB("gff.sqlite")

    # list to append the new ids
    features_mod = []

    # loops on gene and children. Return list of gene and children
    for features in db.iter_by_parent_childs(featuretype="gene"):
        # loop on the features

        for feature in features:
            if feature.featuretype == "exon":
                parent = ",".join(feature.attributes['Parent'])
                # condition for Maize, which has Name instead of IDs
                if 'ID' in feature.attributes:
                    ident = ",".join(feature.attributes['ID'])
                    if ident in feature.id :
                        # remove the current feature from the DB
                        # otherwise we have double ids
                        db.delete(feature)
                        # create the new id
                        ident = f"{feature.id}.mod"
                        # add it in the attributes
                        feature.attributes['ID'] = ident
                        features_mod.append(feature)

                    else :
                        print(ident, feature.id)

                elif 'Name' in feature.attributes : 
                    ident = ",".join(feature.attributes['Name'])
                dict_exons[(feature.chrom, parent)].append([feature.start, feature.end, ident, ".", feature.strand])

    # print(len(features_mod))
    # update the db with the new ids
    if features_mod != [] :
        db = db.update(features_mod)

    # creates introns but doesn't work to add it in the 1st db
    intron_db = db.create_introns()

    for intron in intron_db:
        parent = ",".join(intron.attributes['Parent'])

        # condition for Maize, which has Name instead of IDs
        if 'ID' in intron.attributes:
            ident = ",".join(intron.attributes['ID'])
        else : 
            ident = ",".join(intron.attributes['Name']) 
            # in intron.attributes['Name'] exon10 will be before exon 9 -> saw it too late to change

        if intron.start < intron.end :
            # check if current intron is not already in the dict
            # problem with sweetcorn if not done
            if [intron.start, intron.end, ident] not in dict_introns[(intron.chrom, parent)] :
                dict_introns[(intron.chrom, parent)].append([intron.start, intron.end, ident, ".", intron.strand])

    return dict_exons, dict_introns


def update_chrom_name(dict_chrom_name, dict_feature) :

    dict_chrom_updated = {}

    for orig_chrom, masive_chrom in dict_chrom_name.items() : 
        for chrom_parent, coord in dict_feature.items() :

            if orig_chrom == chrom_parent[0]:
                dict_chrom_updated[(masive_chrom, chrom_parent[1])] = dict_feature[(orig_chrom, chrom_parent[1])]

    return dict_chrom_updated


def write_exons_introns(data_to_write, data_type, species) :

    with open(f"{species}_genes-{data_type}.bed", "w") as output : 
        for chrom_parent, coord in data_to_write.items() :

            chrom = chrom_parent[0]
            # parent = chrom_parent[1]
            # print(coord)

            for coord_couple in coord : 
                coord_couple_to_write = "\t".join(str(x) for x in coord_couple)
                # print(f"{chrom}\t{coord_couple_to_write}\t{parent}")
                output.write(f"{chrom}\t{coord_couple_to_write}\n")


############ MAIN 

parser = argparse.ArgumentParser(description="")
parser.add_argument("--ren", help="")
parser.add_argument("--gff", help="")
args = parser.parse_args()

renamed_file = args.ren
gff_file = args.gff

# renamed_file = "/Users/jd493/Documents/Phylogeny/Data/newgenomes/Acer_truncatum/Acer_truncatum_genome.fa.renamed"
# gff_file = "/Users/jd493/Documents/Phylogeny/Data/newgenomes/Acer_truncatum/Acer_truncatum.gff"


dico_chrom_name, species = parse_renamed_MASiVE(renamed_file)
dict_exons, dict_introns = get_exons_introns(gff_file)
dict_exons_newname = update_chrom_name(dico_chrom_name, dict_exons)
dict_introns_newname = update_chrom_name(dico_chrom_name, dict_introns)
write_exons_introns(dict_exons_newname, "exons", species)
write_exons_introns(dict_introns_newname, "introns", species)

# os.remove("tmp.gff")
