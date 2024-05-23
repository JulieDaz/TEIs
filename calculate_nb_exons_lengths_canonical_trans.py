#-*- coding: utf-8 -*-

import argparse

import gffutils


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


def get_transcript_canonical_intronlen(gff) :

    dict_intron_length = {}
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

    dict_size_transcript = {}
    dict_introns_sum_length = {}
    # loops on gene and children. Return list of gene and children
    for features in db.iter_by_parent_childs(featuretype="gene"):

        # loop on the features
        for feature in features:
            if "chr" in feature.chrom : # to get rid of the scf
                if feature.featuretype == "gene":
                    gene_id = feature.id
                    gene_len = feature.end - feature.start
                    dict_size_transcript.setdefault((gene_id, gene_len), [])

                # retrieve transcript id + calculate transcript len
                # and add to dict
                if feature.featuretype == "mRNA":
                    if "canonical_transcript" in feature.attributes :

                        transcript_id = feature.id
                        transcript_len = feature.end - feature.start
                        dict_size_transcript[(gene_id, gene_len)] = [transcript_id, transcript_len]

    # creates introns but doesn't work to add it in the 1st db
    intron_db = db.create_introns()

    for intron in intron_db:
        parent = ",".join(intron.attributes['Parent'])
        length_intron = intron.end - intron.start
        dict_intron_length.setdefault(parent, [])
        dict_intron_length[parent].append(length_intron)

    for transcript, list_intron_length in dict_intron_length.items() :
        sum_introns = sum(list_intron_length)
        dict_introns_sum_length[transcript] = sum_introns

    return dict_size_transcript, db, dict_introns_sum_length


def get_sum_exons_for_canonical_trans(dict_size_transcript, db, dict_introns_sum_length) :

    '''
    Select the longest transcript for each gene and categorise it in 
    3 categories: 1 exon, 2 exons, >2 exons
    Returns a list [gene, transcript, nb_exons, category]
    '''

    list_canonical_transcript = []

    for gene_len, transcript_id_length in dict_size_transcript.items():

        # strip the _transcript to write in the final list
        transcript_id = transcript_id_length[0].strip("_transcript")
        transcript_len = transcript_id_length[1]

        # need to keep the _transcript to extract the nb of exons
        canonical_transcript = transcript_id_length[0]

        nb_exon = len([exon.attributes['Name'][0] for exon in db.children(canonical_transcript, featuretype="exon")])

        if nb_exon == 1 :
            cat = "exon=1"
        elif nb_exon == 2 :
            cat = "exon=2"
        else :
            cat = "exon>2"

        list_exon_start = [exon.start for exon in db.children(canonical_transcript, featuretype="exon")]
        list_exon_end = [exon.end for exon in db.children(canonical_transcript, featuretype="exon")]
        # for exon in list_exon : 
        list_exon_length = [element1 - element2 for (element1, element2) in zip(list_exon_end, list_exon_start)]
        exon_sum = sum(list_exon_length)

        if canonical_transcript in dict_introns_sum_length :
            introns_sum = dict_introns_sum_length[canonical_transcript]

        list_canonical_transcript.append([gene_len[0], transcript_id, nb_exon, cat, gene_len[1], transcript_len, exon_sum, introns_sum])
    # print(list_canonical_transcript)
    return list_canonical_transcript


def write_output(list_longtest_transcript, output) : 

    with open(f"{output}.canonical.exoncat", "w") as f:
        f.write(f"gene_id\ttranscript_id\tnb_exon\tcat_exon\tgene_length\ttranscript_length\texons_sum\tintrons_sum\n")

        for sublist in list_longtest_transcript :
            sublist_to_write = "\t".join(str(f) for f in sublist)
            f.write(f"{sublist_to_write}\n")


#### MAIN 

parser = argparse.ArgumentParser(description="")
parser.add_argument("-gff", help="")
args = parser.parse_args()

gff_file = args.gff


dico_size_transcript, db, dict_introns_sum_length = get_transcript_canonical_intronlen(gff_file)
longest_transcript_list = get_sum_exons_for_canonical_trans(dico_size_transcript, db, dict_introns_sum_length)
write_output(longest_transcript_list, gff_file)
