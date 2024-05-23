set -e

## server
cd Zea_mays_$1

echo "Starting maize $1"

# local
# cd Zm-$1-REFERENCE-NAM-1.0

## Create all the necessary files
# select and transform TE gff to bed file
python3  ~/maize/pipeline_SVs-introns/select_TEs_from_NAM_gff.py -ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed -gff $1.pseudomolecules-v*.fasta.mod.EDTA.intact.gff3 -type intact

# gene_gff=Zm-$1-REFERENCE-NAM-*1.gff3
## create bed for genes
# python3 ~/maize/pipeline_SVs-introns/prepare_gene_file_Circos-gffutils.py --ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed --gff $gene_gff
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/prepare_gene_file_Circos-gffutils.py --ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed --gff $gene_gff

## create bed for exons/introns
# python3 ~/maize/pipeline_SVs-introns/extract_exons-introns_gffutils.py --ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed --gff $gene_gff
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/extract_exons-introns_gffutils.py --ren Zm-$1-REFERENCE-NAM-1.0.fa.renamed --gff $gene_gff

# create exoncat file that have the info about the a gene, longest transcript, and their lengths + cumulative length of exons, nb of exons, exon category...
# python3  ~/maize/pipeline_SVs-introns/calculate_nb_exons_lengths_canonical_trans.py -gff Zm-$1-REFERENCE-NAM-*_Zm*.1.gff3

# filter the TEs
# Rscript  ~/maize/pipeline_SVs-introns/filtering_TEs.R Zmays-$1_MaizeGDB_intactTEs $2 $3

# # intersect between TEs and genes
# bedtools intersect -a Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3.bed -b Zmays-$1_genes.bed -wao > Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_genes--intersect.txt

# nb_TEs=$(wc Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3.bed | awk '{print $1}')
# nb_genes=$(wc Zmays-$1_genes.bed | awk '{print $1}')

# # categorize the various intersect (SV in Gene, Gene in SV, partial ovl) and create bed of the SVinGene ovl --> /!\ takes about 20 minutes locally
# python3  ~/maize/pipeline_SVs-introns/calculate_ovl.py -i Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_genes--intersect.txt -type TE -transposon $nb_TEs -gene $nb_genes

# # get intron from longest transcript
# Rscript  ~/maize/pipeline_SVs-introns/get_introns_longest_trans.R Zmays-$1_genes-introns.bed Zm-$1-REFERENCE-NAM-1.0_*.gff3.canonical.exoncat Zmays-$1_genes-introns_canonical_trans.bed
# Rscript  ~/maize/pipeline_SVs-introns/get_introns_longest_trans.R Zmays-$1_genes-exons.bed Zm-$1-REFERENCE-NAM-1.0_*.gff3.canonical.exoncat Zmays-$1_genes-exons_canonical_trans.bed

# # intersect between the SVinGene ovl and the introns
# bedtools intersect -a Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_genes--intersect_ovl.bed -b Zmays-$1_genes-introns_canonical_trans.bed -wo > Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_introns_canonical_trans--intersect.txt

# ## calculate perc of ovl, selection based on threshold
# python3  ~/maize/pipeline_SVs-introns/calculate_coverage_ovl.py -i Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_introns_canonical_trans--intersect.txt

# # select unique ovl if several ovl in the same intron
# python3  ~/maize/pipeline_SVs-introns/select_unique_ovl_TEfam.py -i Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_introns_canonical_trans--intersect_ovl_coverage.txt

# # ### Adds category of overlaps between introns and TE in the selected file
# python3  ~/maize/pipeline_SVs-introns/add_cov_cat_selected_intron-TE_ovl.py -i Zmays-$1_MaizeGDB_intactTEs_filtered$2-$3_introns_canonical_trans--intersect_ovl_coverage_selected.txt

## done
echo "Maize $1 is done"

cd ..



# local
# # select and transform TE gff to bed file
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/select_TEs_from_NAM_gff.py -ren ../Data/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.renamed -gff ../Data/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.TE.gff3

# ## create bed for genes -> already done

# ## create bed for exons/introns -> already done

# ## intersect between TEs and genes
# bedtools intersect -a Zmays-B73-v5_MaizeGDB_intactTEs.bed -b Zmays-B73-v5_genes.bed -wao > Zmays-B73-v5_MaizeGDB_intactTEs.bed_genes--intersect.txt

# nb_TEs=$(wc Zmays-B73-v5_MaizeGDB_intactTEs.bed | awk '{print $1}')
# nb_genes=$(wc Zmays-B73-v5_genes.bed | awk '{print $1}')

# # echo $nb_SV
# ## categorize the various intersect (SV in Gene, Gene in SV, partial ovl) and create bed of the SVinGene ovl
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/calculate_ovl.py -i Zmays-B73-v5_MaizeGDB_intactTEs.bed_genes--intersect.txt -type TE -transposon $nb_TEs -gene $nb_genes
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/calculate_ovl.py -i Zmays-B73-v5_MaizeGDB_intactTEs.bed_genes--intersect.txt -type TE -transposon 1267368 -gene 39035

# ## intersect between the SVinGene ovl and the introns
# bedtools intersect -a Zmays-B73-v5_MaizeGDB_intactTEs.bed_genes--intersect_ovl.bed -b Zmays-B73-v5_genes-introns.bed -wo > Zmays-B73-v5_MaizeGDB_intactTEs_introns--intersect.txt

# ## calculate perc of ovl, selection based on threshold
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/calculate_coverage_ovl.py -i Zmays-B73-v5_MaizeGDB_intactTEs_introns--intersect.txt

# ## select unique ovl if several ovl in the same intron
# python3 ~/Documents/Chapter3/Scripts/pipeline_SVs-introns/select_unique_ovl_TEfam.py -i  Zmays-B73-v5_MaizeGDB_intactTEs_introns--intersect_ovl_coverage.txt

