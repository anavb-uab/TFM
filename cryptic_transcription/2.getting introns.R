library(readxl)
library(dplyr)
library(glue)

#Select introns where cTSS is, divide them in 100 fragments

#Input: 
# 1. Excel file with normalised count data of upregulated transcripts by exons, containing at least: Ensembl_ID (gene), Chr, Start, End, Strand, Length, Transcript_ID, Trasncript_name, Exon, Mean_Unt_expression, Mean_Dox_expression
# 2. Document containing the transcript names of the selected transcripts



expression_data <- read_xlsx("1_upregulated_data_normalized.xlsx")
selected_transcripts <- read.table("2_selected_transcripts.txt")

selected_expression_data <- expression_data %>%
  filter(Transcript_name %in% selected_transcripts$V1)


#Select exons where cryptic trasncription starts 
exons_ctss <- selected_expression_data%>%
  #we calculate Dox_expression - Untr_expression for each exon
  mutate("Dox_Untr_expression" = `Mean H1Xsh +Dox FPKM`-`Mean H1Xsh -Dox FPKM`)%>%
  mutate("Dox_Untr_expression" = ifelse(Dox_Untr_expression == 0, 10^-7, Dox_Untr_expression)) %>% #To avoid indeterminations
  mutate("n-n_1"=ifelse(Exon == 1, NA, Dox_Untr_expression - lag(Dox_Untr_expression)))%>%
  mutate("n/n_1"=ifelse(Exon == 1, NA, Dox_Untr_expression / lag(Dox_Untr_expression))) %>%
  mutate("Dox_Dox-1" = ifelse(Exon == 1, NA, (`Mean H1Xsh +Dox FPKM` - lag(`Mean H1Xsh +Dox FPKM`)))) %>%  
  filter(`n-n_1`>0) %>%
  filter(abs(`n/n_1`)>4) %>%
  group_by(Transcript_id) %>% 
  filter(`n-n_1` == max(`n-n_1`)) %>%
  select(Ensembl_id, Chr, Start, End, Strand, Transcript_id, Transcript_name, Exon) 

write.csv(exons_ctss, "selected_exons_cTSS_candidates.csv")


#-----------------------------------------------------------------------------------
#CHECK EXONs BY EYE 
#-----------------------------------------------------------------------------------

#Obtain introns

#Take previous exon
prev_exons_ctss <- read.csv("selected_exons_cTSS_candidates.csv") %>%
  mutate(prev_exon = Exon-1) %>%
  select(Transcript_id, prev_exon)
prev_exons_ctss <- merge(selected_expression_data, prev_exons_ctss, by.x = c("Transcript_id", "Exon"), by.y = c("Transcript_id", "prev_exon")) %>%
  select(Ensembl_id, Chr, Start, End, Strand, Transcript_id, Transcript_name, Exon)


exons_ctss <- read.csv("selected_exons_cTSS_candidates.csv") %>%
  select(Ensembl_id, Chr, Start, End, Strand, Transcript_id, Transcript_name, Exon)

prev_and_current_exon_ctss <- rbind(exons_ctss, prev_exons_ctss)%>%
  arrange(Transcript_id, Exon)


#Calculate intronic region
introns <- prev_and_current_exon_ctss %>%
  mutate("i_length" = ifelse(Strand == "+", Start - lag(End) -2, lag(Start) - End -2)) %>%
  mutate("i_start" = ifelse(Strand == "+", lag(End)+1, End +1)) %>%
  mutate("i_end" = ifelse(Strand == "+", Start-1, lag(Start)-1)) %>%
  mutate("Exon1" = lag(Exon)) %>%
  mutate("Exon2" = Exon) %>%
  filter(row_number() %% 2 == 0) %>%
  select(-Start, -End, -Exon)

introns_map <- introns %>%
  ungroup()%>%
  select(Chr, Strand, i_start, i_end, i_length, Ensembl_id)
colnames(introns_map) <- c("Chr", "Strand", "Start", "End", "Length", "Ensembl_id")

#Save whole introns
introns_bed <- introns_map %>%
  select(Chr, Start, End, Ensembl_id)
#write.table(introns_bed, file = "whole_introns_CheRNAseq.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#---------------------------------------------------------------------------------
#Ignore first/last 15 bp because they might show peak reads because of being close to exon

introns_map <- introns_map %>%
  mutate(Start = Start+15)%>%
  mutate(End = End-15)


#Divide each intron in 100 fragments to find the one where transcription begins
#Function to divide each intron
divide_intron <- function(df){
  fragment_size <- ceiling((df$Length-30) / 100)
  
  Start_positions <- seq(df$Start, df$End, by = fragment_size)
  End_positions <- c(Start_positions[-1] - 1, df$End)  # Make sure the last fragment ends in End position
  
  fragments <- data.frame(
    Chr = df$Chr,
    Strand = df$Strand,
    Start = Start_positions[1:100],  
    End = End_positions[1:100],       
    Ensembl_id = df$Ensembl_id
  )
  fragments <- fragments %>%
    mutate(Ensembl_id = ifelse(Strand == "+", paste0(Ensembl_id, ".", 1:100), paste0(Ensembl_id, ".", 100:1)))
  return(fragments)
}

introns_map_100 <- introns_map %>%
  filter(Length>1000)

fragments_map <- data.frame(data.frame(Chr=character(), Strand=character(), Start=integer(), End=integer()), Ensembl_id=character())
for (i in 1:nrow(introns_map_100)) {
  div <- divide_intron(introns_map_100[i,])
  fragments_map <- rbind(fragments_map, div)
}

col_order <- c("Ensembl_id", "Chr", "Start", "End", "Strand")
fragments_map <- fragments_map[,col_order]
fragments_map <- na.omit(fragments_map)

colnames(fragments_map)<- c("GeneID", "Chr", "Start", "End", "Strand")

write.table(fragments_map, file = "cTSS_CheRNAseq_candidates_introns_100fragments.SAF", quote = FALSE, row.names = FALSE, col.names = T, sep = "\t")

#------------------------------------------------------------------------------
#Next step: featureCounts