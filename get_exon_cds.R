library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(seqinr)
library(AnnotationHub)

MASTER_FOLDER <- '.'
HOMO_SAPIENS_SL_PATH <- file.path(MASTER_FOLDER, 'genomics_ensembl_database.csv')
genes_database <- read.csv(HOMO_SAPIENS_SL_PATH)


### Function ###
geneNameTogeneIdENSG <- function(gene_name) {
  gene_id_ <- genes(edb, 
                    filter = list(
                      GeneNameFilter(gene_name),
                      GeneIdFilter("ENSG", condition = "startsWith"),
                      TxIsCanonicalFilter(1, "==")
                    ), 
                    return.type = "DataFrame")
  if (nrow(gene_id_) == 0) {
    return ("")
  }
  gene_id <- gene_id_[gene_id_$seq_coord_system == 'chromosome',][[1]]
  return (gene_id)
}

geneNameTogeneIdENST <- function(gene_name, 
                                  edb_transcripts) {
  "
  Function to change gene name to gene id
  
  Args:
  gene_name -- Gene name to be querried
  edb_transcripts 
  
  Returns:
  gene_id (str) -- in ENST format 
  "
  gene_id_query <- as.data.frame(edb_transcripts[edb_transcripts$gene_name == gene_name])
  enst_gene_id <- gene_id_query[1, 'tx_id_version']
  return (enst_gene_id)
}

get3UTRinfo <- function(gene_id_enst,
                        txdb,
                        threeUTRDatabase) {
  "
  Get 3' UTR information as a dataframe from an ENST gene id 
  
  Args:
  gene_id_enst -- Gene id in ENST format
  txdb -- Transcript database
  threeUTRDatabase -- 3' UTR database 
  
  Returns:
  A dataframe consists of information on 3' UTR for the given gene name
  "
  threeUTRDatabase <- as.data.frame(threeUTRsByTranscript(txdb,
                                                          use.names=TRUE))
  querried_3utr_info <- threeUTRDatabase[threeUTRDatabase$group_name == gene_id_enst,]
  querried_3utr_info_df <- as.data.frame(querried_3utr_info)
  return (querried_3utr_info_df)
} 

get5UTRinfo <- function(gene_id_enst,
                        txdb,
                        fiveUTRDatabase) {
  "
  Get 5' UTR information as a dataframe from an ENST gene id 
  
  Args:
  gene_id_enst -- Gene id in ENST format
  txdb -- Transcript database
  fiveUTRDatabase -- 5' UTR database 
  
  Returns:
  A dataframe consists of information on t' UTR for the given gene name
  "
  fiveUTRDatabase <- as.data.frame(fiveUTRsByTranscript(txdb,
                                                        use.names=TRUE))
  querried_5utr_info <- fiveUTRDatabase[fiveUTRDatabase$group_name == gene_id_enst,]
  querried_5utr_info_df <- as.data.frame(querried_5utr_info)
  return (querried_5utr_info_df)
}

############### Main ################
# Human genome library hg38
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotation hub
ah <- AnnotationHub::AnnotationHub(cache = './cache')
qr <- ah['AH116291']

# edb
edb <- qr[[1]]
edb <- filter(edb, filter=~tx_is_canonical == TRUE)

# DNA genome 2 bit file
dna <- getGenomeTwoBitFile(edb)

# Transcripts
txs <- transcripts(edb, columns = c("tx_id", "tx_biotype", "tx_id_version", "gc_content", "gene_name", "gene_id"))

# 3' UTR database 
utr3 <- as.data.frame(threeUTRsByTranscript(txdb,
                                            use.names=TRUE))
# 5' UTR database 
utr5 <- as.data.frame(fiveUTRsByTranscript(txdb,
                                           use.names=TRUE))

gene_name_example <- "EXOSC4"
# Get ENSG gene id from gene name
gene_id_ensg <- geneNameTogeneIdENSG(gene_name_example)
# Get ENST gene id from gene name 
gene_id_enst <- geneNameTogeneIdENST(gene_name_example, 
                                     txs)

# Get 3' UTR information as a dataframe 
utr3_info_df <- get3UTRinfo(gene_id_enst,
                            txdb,
                            utr3)
utr3_info_df

# Get 5' UTR information as a dataframe 
utr5_info_df <- get5UTRinfo(gene_id_enst,
                            txdb,
                            utr5)

# Get exon 
exons_seq <- exons(edb,
                   columns=listColumns(edb),
                   return.type = "GRanges",
                   filter=AnnotationFilterList(
                         GeneIdFilter(gene_id_esng),
                         TxIsCanonicalFilter(1, "==")
                   )
)

exon_cds_seq_df <- BSgenome::getSeq(dna, 
                                    exons_seq,
                                    start=utr5_info_df[1, 'end'],
                                    end=utr3_info_df[1, 'start']) 
exon_cds_seq_df
