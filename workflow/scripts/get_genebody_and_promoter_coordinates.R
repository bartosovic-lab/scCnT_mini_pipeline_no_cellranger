library(GenomicRanges)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(argparse)

# Define command line arguments
parser <- ArgumentParser()
parser$add_argument("--out", help="Output file", required=TRUE)
parser$add_argument("--genome", help="Genome version [currently available only hg38]", default="hg38")
args <- parser$parse_args()




CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

annotations                <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations)     <- paste0('chr', seqlevels(annotations))
annotations                <- Signac::Extend(annotations, upstream = 2000)
annotations.collapsed      <- CollapseToLongestTranscript(annotations) 
annotations.collapsed$name <- annotations.collapsed$gene_name
rtracklayer::export(con=args$out, obj = annotations.collapsed, format = "bed")
saveRDS(annotations.collapsed, file = gsub(pattern = ".bed", replacement = ".rds", x = args$out))

