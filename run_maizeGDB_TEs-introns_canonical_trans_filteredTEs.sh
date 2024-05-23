set -e

# server
maize=$(ls -d Z*/)

for m in $maize
do
    echo $m
    line=$(echo $m | cut -d_ -f3 | cut -d\/ -f1)
    echo $line
    if [ "$line" != "B73-v5" ] ; then
        # bash ~/maize/pipeline_SVs-introns/maizeGDB_TEs-introns.sh $line
        bash ~/maize/pipeline_SVs-introns/maizeGDB_TEs-introns_canonical_trans_filteredTEs.sh $line 100 1000
    fi
done

## for B73
# server
cd Zea_mays_B73-v5

## select and transform TE gff to bed file
# python3 ~/maize/pipeline_SVs-introns/select_TEs_from_NAM_gff.py -ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed -gff Zm-$1-REFERENCE-NAM-1.0.TE.gff3

# create bed for genes -> already done

# create bed for exons/introns -> already done

# create exoncat file that have the info about the a gene, longest transcript, and their lengths + cumulative length of exons, nb of exons, exon category...
python3 ~/maize/pipeline_SVs-introns/calculate_nb_exons_lengths_canonical_trans.py -gff Zm-$1-REFERENCE-NAM-*_Zm*.1.gff3

## filter the TEs
Rscript ~/maize/pipeline_SVs-introns/filtering_TEs.R Zmays-B73-v5_MaizeGDB_allTEs 100 1000

## intersect between TEs and genes
bedtools intersect -a Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000.bed -b Zmays-B73-v5_genes.bed -wao > Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_genes--intersect.txt

nb_TEs=$(wc Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000.bed | awk '{print $1}')
nb_genes=$(wc Zmays-B73-v5_genes.bed | awk '{print $1}')

## categorize the various intersect (SV in Gene, Gene in SV, partial ovl) and create bed of the SVinGene ovl --> /!\ takes about 20 minutes locally
python3 ~/maize/pipeline_SVs-introns/calculate_ovl.py -i Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_genes--intersect.txt -type TE -transposon $nb_TEs -gene $nb_genes

## get intron from longest transcript
Rscript ~/maize/pipeline_SVs-introns/get_introns_longest_trans.R Zmays-B73-v5_genes-introns.bed Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.canonical.exoncat Zmays-B73-v5_genes-introns_canonical_trans.bed

## intersect between the SVinGene ovl and the introns
bedtools intersect -a Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_genes--intersect_ovl.bed -b Zmays-B73-v5_genes-introns_canonical_trans.bed -wo > Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_introns_canonical_trans--intersect.txt

# ## calculate perc of ovl, selection based on threshold
python3 ~/maize/pipeline_SVs-introns/calculate_coverage_ovl.py -i Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_introns_canonical_trans--intersect.txt

# select unique ovl if several ovl in the same intron
python3 ~/maize/pipeline_SVs-introns/select_unique_ovl_TEfam.py -i  Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_introns_canonical_trans--intersect_ovl_coverage.txt
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/select_unique_ovl_TEfam.py -i  Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_introns--intersect_ovl_coverage.txt

# ### Adds category of overlaps between introns and TE in the selected file
python3 ~/maize/pipeline_SVs-introns/add_cov_cat_selected_intron-TE_ovl.py -i Zmays-B73-v5_MaizeGDB_allTEs_filtered100-1000_introns_canonical_trans--intersect_ovl_coverage_selected.txt

echo Maize B73 finish

