library(biomaRt)

#use melphalan_genes_with_drivers.txt has updated genes
#mel-focused melphalan_genes_updated_03152026.txt
genes <- readLines("melphalan_genes_with_drivers.txt")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

# Keep only autosomes, add chr prefix, sort
#coords <- coords[coords$chromosome_name %in% as.character(1:22), ]
#coords$chromosome_name <- paste0("chr", coords$chromosome_name)
#coords <- coords[order(as.numeric(sub("chr", "", coords$chromosome_name)), coords$start_position), ]

#alternatively keep all chromosomes
#autosomes was to format the bed file for monopogen, but I need all for my WES data

##debugging previous code, need to version control this
#keep std chroms (1-22, X, Y, MT), do not want patch etc
standard_chrs <- c(as.character(1:22), "X", "Y", "MT")
coords <- coords[coords$chromosome_name %in% standard_chrs, ]

#remove dups
coords <- coords[!duplicated(coords$hgnc_symbol), ]

#manually add gstt1 because biomart is given alt chrom location
gstt1_row <- data.frame(
  hgnc_symbol = "GSTT1",
  chromosome_name = "22",
  start_position = 270308,
  end_position = 278486
)
coords <- rbind(coords, gstt1_row)

#all chroms
coords$chromosome_name <- paste0("chr", coords$chromosome_name)
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrMT")
coords$chromosome_name <- factor(coords$chromosome_name, levels = chr_order)
coords <- coords[order(coords$chromosome_name, coords$start_position), ]
coords$chromosome_name <- as.character(coords$chromosome_name)

# Write BED file
write.table(
  coords[, c("chromosome_name", "start_position", "end_position", "hgnc_symbol")],
  file = "melphalan_genes_complete_final.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Report any missing genes
missing <- setdiff(genes, coords$hgnc_symbol)
cat("Genes not found:", paste(missing, collapse=", "), "\n")
