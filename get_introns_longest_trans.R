library(tidyverse)

suppressMessages(library(tidyverse))

args=commandArgs(trailingOnly = TRUE)

all_introns_file = args[1]
longest_transcript_file = args[2]
output_name = args[3]

# all_introns = read.table("Data/Zm-B73-REFERENCE-NAM-5.0/Zmays-B73-v5_genes-introns.bed",
#                          sep = "\t")

## read the introns bed files
all_introns = read.table(all_introns_file,
                         sep = "\t")

## add col names
colnames(all_introns) = c("chr", "start", "end", "intron_id", "score", "strand")

## get the transcript id, from the introns id
all_introns_mutate <- all_introns %>% 
  mutate(transcript_id = str_split(intron_id, "[.]", simplify = T)[,1])

# longest_transcript = read.table("Data/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.exoncat",
#                                 sep = "\t", header = T)

## read file that contains the longest transcript
longest_transcript = read.table(longest_transcript_file,
                                sep = "\t", header = T)

## join the info with the transcript id with the introns coord
introns_longest_trans = semi_join(all_introns_mutate, longest_transcript) %>% 
  select(-transcript_id) %>% 
  as.data.frame()

write.table(introns_longest_trans, file = output_name, quote = F,
            sep = "\t", row.names = F, col.names = F)


