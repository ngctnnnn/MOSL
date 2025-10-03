# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(seqinr)
library(AnnotationHub)
# choosebank("genbank")


MASTER_FOLDER <- '.'
HOMO_SAPIENS_SL_PATH <- file.path(MASTER_FOLDER, 'genomics_ensembl_database.csv')
genes_database <- read.csv(HOMO_SAPIENS_SL_PATH)

## Human genome library hg38
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## 3 UTR
utr3 <- as.data.frame(threeUTRsByTranscript(txdb,
                                      use.names=TRUE))
utr3

## 5 UTR
utr5 <- as.data.frame(fiveUTRsByTranscript(txdb,
                                            use.names=TRUE))
utr5

## Annotation hub
ah <- AnnotationHub::AnnotationHub(cache = './cache')
qr <- ah['AH116291']
qr

## edb
edb <- qr[[1]]
edb <- filter(edb, filter=~tx_is_canonical == TRUE)

## DNA genome 2 bit file
dna <- getGenomeTwoBitFile(edb)

## Transcripts 
txs <- transcripts(edb, columns = c("tx_id", "tx_biotype", "tx_id_version", "gc_content", "gene_name", "gene_id"))
txs

geneNameTogeneId <- function(gene_name, 
                             edb_transcripts) {
  "
  Function to change gene name to gene id
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
  "
  fiveUTRDatabase <- as.data.frame(fiveUTRsByTranscript(txdb,
                                                          use.names=TRUE))
  querried_5utr_info <- fiveUTRDatabase[fiveUTRDatabase$group_name == gene_id_enst,]
  querried_5utr_info_df <- as.data.frame(querried_5utr_info)
  return (querried_5utr_info_df)
}

########################
# Main
txs <- transcripts(edb, 
                   columns = c("tx_id", 
                               "tx_biotype", 
                               "tx_id_version", 
                               "gc_content", 
                               "gene_name", 
                               "gene_id")
                   )
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

gene_name_example <- "EXOSC4"

utr3 <- as.data.frame(threeUTRsByTranscript(txdb,
                                            use.names=TRUE))
utr5 <- as.data.frame(fiveUTRsByTranscript(txdb,
                                            use.names=TRUE))

gene_id_enst <- geneNameTogeneId(gene_name_example, txs)

utr3_info_df <- get3UTRinfo(gene_id_enst,
                             txdb,
                             utr3)

utr5_info_df <- get5UTRinfo(gene_id_enst,
                            txdb,
                            utr5)
utr5_info_df
