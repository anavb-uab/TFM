##Cryptic transcription detection script - output: exons where cTSS is manifested

library(dplyr)
library(readxl)
library(writexl)



#1. first we have to calculate the mean expression of dox and untr samples in the upregulated excel
#Input: 
# 1. Excel file with normalised count data of upregulated transcripts by exons, containing at least: Ensembl_ID (gene), Chr, Start, End, Strand, Length, Transcript_ID, Trasncript_name, Exon, Mean_Unt_expression, Mean_Dox_expression
# 2. Document with information about differential expression analysis, containing at least: Ensembl_ID (gene), log2FoldChange

#Output: list of selected transcripts - candidates for cryptic transcription

fc_data <- read.table("2_differential_expression_analysis_data.txt", header = TRUE) 

data <- read_excel("1_upregulated_data_normalized.xlsx")
clean_data <- select(data, "Ensembl_id", "Transcript_name", "Exon", "Start", "End", "Strand", `Mean -Dox FPKM`, `Mean +Dox FPKM`) %>%
  filter(Ensembl_id %in% fc_data$Ensembl_ID)
sum(is.na(clean_data)) #Check there are no NAs


## From the upregulated genes we take the longest transcript variant in exon number

higher_exon_num <-clean_data %>%
  group_by(Ensembl_id, Transcript_name)%>%
  summarize(n_exones = n(), .groups = 'drop_last') %>%
  slice_max(n_exones, with_ties = FALSE)

higher_exon_num <- data%>%
  inner_join(higher_exon_num, by = c("Ensembl_id", "Transcript_name"))

#Number of upregulated genes
length(unique(higher_exon_num$Ensembl_id))

#Length of upregulated genes distribution
hist(higher_exon_num$n_exones, breaks = max(higher_exon_num$n_exones))
abline(v=4, col="orchid")



################################################################################
################################################################################

## Preselection

wrk <- higher_exon_num %>%
  # Select transcripts with more than 3 exons
  filter(n_exones >= 4) %>%   
  select(Ensembl_id, Chr, Length,Start, End, Strand, Transcript_id, Transcript_name, Exon, n_exones, `Mean +Dox FPKM`, `Mean -Dox FPKM`) %>%
  #Filter out when dox_exp max is in exon 1
  group_by(Transcript_id) %>%
  mutate("max_dox" = max(`Mean +Dox FPKM`)) %>%
  filter(max_dox != `Mean +Dox FPKM`[Exon==1]) %>%
  select(-max_dox) %>%
  ungroup() %>%
  #Calculate Dox_expression - Untr_expression for each exon
  mutate("Dox_Untr_expression" = `Mean +Dox FPKM`-`Mean -Dox FPKM`)%>%
  mutate("Dox_Untr_expression" = ifelse(Dox_Untr_expression == 0, 10^-7, Dox_Untr_expression))  #To avoid indeterminations

#Number of pre-selected genes
length(unique(wrk$Ensembl_id))

#-------------------------------------------------------------------------------

#Filter StartTogether: 
wrk_StartTog <- wrk %>% 
  group_by(Transcript_id) %>%
  mutate("Max_DOX" = max(`Mean +Dox FPKM`))%>%
  mutate("Dox_untr_exon1" = Dox_Untr_expression[Exon ==1]) %>%
  filter(ifelse(Max_DOX<0.5, Dox_untr_exon1/Max_DOX <= 0.15, Dox_untr_exon1/Max_DOX <= 0.09)) ###FILTER StartTogether

#List of trasncripts and genes that pass filter StartTogether
startTog_Trname <- unique(wrk_StartTog$Transcript_name)
startTog_EnsID <- unique(wrk_StartTog$Ensembl_id)

#List of trasncripts that do not pass filter StartTogether
left_out_filterA <- setdiff(wrk$Transcript_name, wrk_StartTog$Transcript_name)


#-------------------------------------------------------------------------------

#Filter LastVSFirst:

#Prepare the data for the filter
wrk_efirst_elast <- wrk %>%
  group_by(Ensembl_id) %>%
  filter(Exon == 1 | Exon == n_exones | Exon == n_exones-1)%>% #Keep first two and last two exons
  mutate("exon_pos" = ifelse(Exon==1, "First", "Last")) %>%
  ungroup()%>%
  group_by(Transcript_id, exon_pos) %>%
  mutate(average_increment_pos = mean(Dox_Untr_expression)) %>%
  ungroup()

## Create the ratio of Δ(expression last) / Δ(expression first) and their log
wrk_lastVSfirst <- wrk_efirst_elast %>%
  mutate(ratios = average_increment_pos / lag(average_increment_pos)) %>%
  filter(row_number() %% 3 == 2 & average_increment_pos > lag(average_increment_pos)) %>% #∆ last > ∆ first  
  mutate(log_ratios = log(abs(ratios))) %>% 
  filter(abs(log_ratios)>1) #Filter 1 LastVSFirst

#List of trasncripts and genes that pass filter LastVSFirst
lastVSfirst_Trname <- unique(wrk_lastVSfirst$Transcript_name)
lastVSfirst_EnsID <- unique(wrk_lastVSfirst$Ensembl_id)

#List of trasncripts and genes that do not pass filter LastVSFirst
left_out_filter1 <- setdiff(wrk$Transcript_name, wrk_lastVSfirst$Transcript_name)

#-------------------------------------------------------------------------------

#Select the intersection between both filters

intersect_a1_Trname <- intersect(startTog_Trname, lastVSfirst_Trname)
intersect_a1_EnsID <- intersect(startTog_EnsID, lastVSfirst_EnsID)

wrk_filtered <- wrk %>%
  filter(Ensembl_id %in% intersect_a1_EnsID)

#-------------------------------------------------------------------------------

#Order the candidates by FoldChange

fc_filtered <- fc_data %>%
  filter(Ensembl_ID %in% intersect_a1_EnsID) %>%
  select(Ensembl_ID, log2FoldChange)

wrk_filtered <- merge(wrk_filtered, fc_filtered, by.x = "Ensembl_id", by.y = "Ensembl_ID", all.x = TRUE)
wrk_filtered <- wrk_filtered[order(-wrk_filtered$log2FoldChange), ]

#List of trasncript IDs and names that pass both filters
intersect_TrID_sorted <- unique(wrk_filtered$Transcript_id)
intersect_Trname_sorted<- unique(wrk_filtered$Transcript_name)

write.table(intersect_Trname_sorted, "Candidate_trasncripts_cTSS.txt", quote = F, col.names = F)

#-------------------------------------------------------------------------------

#Now plot the exon expression (plot script) to screen by eye and select the candidates of interest

#-------------------------------------------------------------------------------

#Selected candidates GOOD

#CheRNAseq
chH1X_good_ajv <- c('KNG1-001', 'CTNNA3-001', 'RP11-33A14.1-001', 'COL21A1-001', 'EML5-001', 'RTTN-001', 'UGT2B10-001', 'MAPKAP1-201', 'ACE2-201', 'MCU-005', 'PPP1R1B-013', 'DGKB-001', 'AC003958.2-002', 'GUCY1A2-002', 'NLGN1-001')
chH12_good_ajv <- c('HCAR1-201', 'VCAN-001', 'KCNH6-002', 'CCDC39-002', 'LRRC36-001', 'KIFAP3-201', 'CRYM-001', 'EFHB-002', 'PLB1-012', 'AL450307.1-201', 'CNTNAP3-001', 'RP11-986E7.7-001', 'INTS4L2-002', 'AXDND1-004', 'TENM4-001', 'SRD5A2-201', 'UMODL1-002', 'DUOX1-006', 'ATP1A4-001', 'CD4-001', 'LMX1A-001', 'KRT81-001', 'ADAMTS14-001', 'IL13RA2-001', 'TAT-001', 'PNLDC1-001', 'ST6GALNAC5-002') 
chH14_good_ajv <- c('ARPP21-001', 'PAPOLG-001', 'B4GALT2-001', 'SETD3-001', 'DENND1B-008')

#RNAseq
rnH1X_good_ajv <- c('KNG1-001', 'CTNNA3-001', 'COL21A1-001', 'EML5-001', 'RTTN-001', 'MAPKAP1-201', 'NLGN1-001', 'ATF7IP-001', 'STXBP3-001', 'TRIM24-001')
rnH12_good_ajv <- c('AXDND1-004', 'CRYM-001', 'KCNH6-002', 'KIFAP3-201', 'RP5-1198O20.4-001', 'TGM3-001', 'GRB14-001', 'TFAP2B-201')
rnH14_good_ajv <- c('AADAC-001', 'ARHGAP26-001', 'ARPP21-001', 'B4GALT2', 'CATSPERB-001', 'CLIP4-001', 'DGKB-001', 'HGD-001', 'ITGB8-005', 'LRRTM4-001', 'PCDH15-007', 'DENND1B-008', 'PAPOLG-001', 'SETD3-001')
rnH13_good_ajv <- c('BMPR1B-001', 'CMYA5-001', 'FGD4-004', 'FLRT2-005', 'GBP2-001', 'MX1-001')

mult_ISG_good_ajv <- c('APOL3-002', 'CFB-201', 'HLA-F-001', 'MX1-001', 'CSF1-201', 'ETV6-001', 'SNTB2-001')
mult_noISG_good_ajv <- c('BEND7-001', 'CNTN6-001', 'CRYAB-015', 'CUBN-001', 'DAPK1-201', 'DPP10-002', 'FMO6P-003', 'PTPRB-001', 'THEMIS-001', 'ZNF326-001')

list_all_variants <- ls(pattern = "good_ajv$")

for (name in list_all_variants) {
  data <- get(name)
  
  # Crear nombre de archivo
  file_name <- paste0("cTSS_candidates_", name, ".txt")
  write.table(data, file = file_name, row.names = F, col.names = F, quote = FALSE)
  
}

#Selected candidates MEDIUM

#CheRNAseq
chH1X_medium_ajv <- c("RAMP2-001", "DDR2-001", "CNTD1-001", "STMN3-001", "CTB-32H22.1-001", "CD3E-001", "ALOX15B-001", "BLM-001", "CGNL1-001","TRIB3-201", "DMGDH-001", "GNG7-003", "RP11-986E7.7-001", "STXBP3-001", "TRIM24-001", "KLHL24-003", "ASNS-201", "TSC22D3-005", "RAB11FIP4-001")
chH12_medium_ajv <- c("DNAH3", "AC017050.1-001" , "DNAJB13-001", "TUBA4B-001", "TTC18-001", "NTM-201", "STON2-002", "CSMD3-005", "LRRC23-001", "CAPSL-001", "ZMYND10-001", "RP11-53M11.3-001", "PPIL6-003", "RP11-163N6.2-003", "DZIP1L-001", "AMPH-001", "TPPP-001", "CCDC158-002", "FAM81A-001", "KIAA1456-002", "ULK4-001", "AKAP14-004", "CCDC157-001", "C6-001", "ZNF487-201", "KCNMB1-001", "RP11-188C12.3-001", "ICAM2-004", "GPR142-001", "CSGALNACT1-001", "SUSD4-203", "RP11-379F4.4-003", "PTPRB-001", "NPHS2-001", "SLC51B-001", "CXCR2-002", "TBATA-001", "NEU4-002", "PTPRT-005")
chH14_medium_ajv <- c("CNTN6-001", "TAF5L-004", "LUZP2-001", "CDHR4-001", "CLIP4-001", "DCDC1-010", "DNAH3-001", "LAMP3-001", "RP11-267N12.3-009", "RP6-109B7.3-002", "SCG5-001", "TAPT1-AS1-004", "VCAN-001")

#RNAseq
rnH1X_medium_ajv <- c("DDR2-001", "MCU-005", "ATF2-003", "CLPB-006", "FCHO2-001", "PLEKHA5-004")
rnH12_medium_ajv <- c("CCDC108-001", "DDAH2-001", "GPD1-001", "HCAR1-201", "KLK10-001", "KMO-001", "NR5A2-001", "NTM-201", "FRMPD3-201", "TTBK1-001", "AZGP1-001", "C17orf72-001", "STK39-001")
rnH14_medium_ajv <- c("OXCT1-001", "PDE3A-001", "SCG5-001", "EPHA4-001", "WARS2-001", "ZFAT-002")
rnH13_medium_ajv <- c('CECR7-001', 'CPNE4-001', 'KIAA1217-004', 'TBC1D8B-001', 'TNIK-001', 'APC2-001', 'BRSK2-201', 'CC2D2A-002', 'DACH1-002', 'RBMS2-001', 'TRPC1-002', 'TLSP-002', 'VWDE-005')

mult_ISG_medium_ajv <- c("ASPHD2-001", "CFB-001", "CIITA-001", "HLA-DPA1-001", "HLA-F-001", "IFI27-013", "STAT4-001", "XAF1-001", "APOL2-001", "KCNG1-001", "LGALS3BP-001", "MAP1B-001", "RNF144A-001", "SLC15A3-001", "TRANK1-201")
mult_noISG_medium_ajv <- c("CCDC28A-001", "CDYL2-002", "CHRNB2-001", "CSAG2-001", "CSAG3-001", "EFTUD1-002", "FAM83A-202", "FLIP1L-001", "FRAS1-201", "IGSF9B-002", "KRT83-001", "MTSS1-005", "PODXL-005", "RXRG-001", "TMC7-001", "TNFRSF11A-001", "TOX3-002")

list_all_med_variants <- ls(pattern = "medium_ajv$")

for (name in list_all_med_variants) {
  data <- get(name)
  
  # Crear nombre de archivo
  file_name <- paste0("cTSS_candidates_", name, ".txt")
  write.table(data, file = file_name, row.names = F, col.names = F, quote = FALSE)
  
}

