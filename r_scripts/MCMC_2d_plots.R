options(stringsAsFactors = FALSE)

# Install dependencies
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("ggrastr")
safelyLoadAPackageInCRANorBioconductor("ggh4x")


wd <- commandArgs(TRUE)[1]
# wd <- "/scratch/ldelisle/baredSC_Darbellay/"
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# directory <- "baredSC_2d"
# The table is:
# input\tgene\txmax
table.fn <- commandArgs(TRUE)[3]
# table.fn <- "/home/ldelisle/softwares/scriptsAroundBaredSCForDarbellayEtAl2024/tables/table_baredSC_2d.txt"
# output prefix
output.prefix <- commandArgs(TRUE)[4]
# output.prefix <- "/home/ldelisle/softwares/scriptsAroundBaredSCForDarbellayEtAl2024/plots/2d"

my.order <- c("P", "PM", "C", "IM", "D", "T", "PO", "MCT")

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h = F,
                       col.names = c("full.path.input", "genex", "geney", "xmax", "ymax"))
my.table$i <- rownames(my.table)

my.table$input <- basename(my.table$full.path.input)

my.table$genes <- paste0(my.table$genex, "VS", my.table$geney)

# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf2d_flat.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf2d_flat.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "genes", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "genes")]), unique(meta.data[, c("i", "genes")]), by = "i")
if (!all(temp$genes.x == temp$genes.y)) {
  # The gene had "_" in its name:
  meta.data$genes <- my.table$genes[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("file", "genes")], 1, function(v){
    return(paste(strsplit(gsub("_pdf2d_flat.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$genes, "__", meta.data$info)
}
# Process genes
meta.data <- cbind(meta.data, matrix(unlist(strsplit(meta.data$genes, "VS")), ncol = 2, byrow = T))
colnames(meta.data)[-1:0 + ncol(meta.data)] <- c("genex", "geney")
meta.data$input <- my.table$input[match(meta.data$i, my.table$i)]
# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$value <- df$mean
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Get the corr
meta.data$file.cor <- gsub("pdf2d_flat", "corr", meta.data$file)
corr <- do.call(rbind, lapply(meta.data$file.cor, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file.cor <- fn
  return(df)
}))
# Add the meta data
corr <- merge(corr, meta.data)
# A label is formatted:
corr$label <- paste0(round(corr$mean, 2), "[-", round(corr$mean, 2) - round(corr$low, 2), "]",
                     "^{+", round(corr$high, 2) - round(corr$mean, 2), "}")
# For the p-value, a superior value is given (still only mean + sd)
corr$p.label <- sapply(corr$pval + corr$error, function(v){paste0("p<", format(v, digits = 2))})

cluster.name <- gsub("cluster_", "", gsub(".txt$", "", unique(my.table$input)))
names(cluster.name) <- unique(my.table$input)
pdfs$cluster <- cluster.name[pdfs$input]
corr$cluster <- cluster.name[corr$input]
names(cluster.name) <- toupper(sapply(sapply(strsplit(cluster.name, " |-"), substr, 1, 1), paste, collapse = ""))
pdfs$cluster <- factor(pdfs$cluster, levels = cluster.name[my.order])
corr$cluster <- factor(corr$cluster, levels = cluster.name[my.order])

my.corr <- subset(corr, cluster %in% c("Chondrocytes", "Perichondrium"))

g <- ggplot(pdfs, aes(x, y)) +
  geom_tile_rast(aes(fill = log(1 + value))) +
  geom_text(
    data = my.corr,
    aes(label = label),
    x = 1, y = 3.5,
    size = 2.5,
    parse = T
  ) +
  geom_text(
    data = my.corr,
    aes(label = p.label),
    x = 1, y = 2.7,
    size = 2.5
  ) +
  facet_wrap2(cluster ~ .) +
  theme_minimal() +
  ylab(unique(pdfs$geney)) +
  xlab(unique(pdfs$genex)) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 2.5)) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
    legend.position = "none",
    axis.title = element_text(face = "italic")
  )

dir.create(dirname(output.prefix), showWarnings = FALSE)
ggsave(paste0(output.prefix, ".pdf"), g,  width = 6.5, height = 3.7)
