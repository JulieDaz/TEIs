### filter the TEs
### script use in the pipeline to filter the length of TEs

library(tidyverse)

args=commandArgs(trailingOnly = TRUE)

bed_file = args[1] ## don't give the .bed
length_other_TEs = as.numeric(args[2])
length_LTR = as.numeric(args[3])

print(length_other_TEs)
print(length_LTR)
# length_LTR = 500
# length_other_TEs = 100

TEs_bed = read.table(paste(bed_file, ".bed", sep = ""))

# TEs_bed = read.table("Data/Zm-B73-REFERENCE-NAM-5.0/Zmays-B73-v5_MaizeGDB_allTEs.bed") 

names(TEs_bed) = c("chr", "start", "end", "fam_ident")
# head(TEs_bed)
# dim(TEs_bed) # 1267368

TEs_bed_fam <- TEs_bed %>% 
  mutate(TE_fam = str_split(fam_ident, ";", simplify = T)[, 1]) %>%
  mutate(len_TE = end - start) %>% 
  as.data.frame()

copias = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "LTR_retrotransposon")

TEs_bed_fam_filtered = TEs_bed_fam %>% 
  filter(ifelse(TE_fam %in% copias, len_TE >= length_LTR, len_TE >= length_other_TEs))

TEs_bed_fam_filtered_to_write = TEs_bed_fam_filtered %>% 
  select(-c(TE_fam, len_TE))

# dim(TEs_bed_fam_filtered_to_write)

write.table(x = TEs_bed_fam_filtered_to_write, 
            file = paste(bed_file, "_filtered", length_other_TEs, "-", length_LTR, ".bed", sep = ""),
            quote = F,
            row.names = F,
            col.names = F, 
            sep = "\t")



