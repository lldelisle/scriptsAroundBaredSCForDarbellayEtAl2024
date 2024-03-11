options(stringsAsFactors = F)

# Install dependencies
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("ggh4x")
safelyLoadAPackageInCRANorBioconductor("grid")

wd <- commandArgs(TRUE)[1]
# wd <- "/scratch/ldelisle/baredSC_Darbellay/"
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# directory <- "baredSC_1d"
# The table is:
# input\tgene\txmax
table.fn <- commandArgs(TRUE)[3]
# table.fn <- "/home/ldelisle/softwares/scriptsAroundBaredSCForDarbellayEtAl2024/tables/table_baredSC_1d.txt"
# output prefix
output.prefix <- commandArgs(TRUE)[4]
# output.prefix <- "/home/ldelisle/softwares/scriptsAroundBaredSCForDarbellayEtAl2024/plots/clusters"

rows.to.use <- c(1:16)
my.order <- c("P", "PM", "C", "IM", "D", "T", "PO", "MCT")

setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h = FALSE,
                       col.names = c("full.path.input", "gene", "xmax"))[rows.to.use, ]
my.table$i <- rownames(my.table)
my.table$input <- basename(my.table$full.path.input)
# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "gene", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "gene")]), unique(meta.data[, c("i", "gene")]), by = "i")
if (!all(temp$gene.x == temp$gene.y)) {
  # The gene had "_" in its name:
  meta.data$gene <- my.table$gene[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("file", "gene")], 1, function(v){
    return(paste(strsplit(gsub("_pdf.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)
}
# I put the order as in my.table
my.table$gene <- factor(my.table$gene, levels = unique(my.table$gene))
meta.data$gene <- factor(meta.data$gene, levels = levels(my.table$gene))

# Subset the meta.data
meta.data <- subset(meta.data, i %in% my.table$i)

# Add the input
meta.data$input <- my.table$input[match(meta.data$i, my.table$i)]

# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)


short.name <- toupper(sapply(sapply(strsplit(gsub(".txt", "", gsub("cluster_", "", unique(my.table$input))), " |-"), substr, 1, 1), paste, collapse = ""))
names(short.name) <- unique(my.table$input)
pdfs$short.name <- factor(short.name[pdfs$input], levels = my.order)
size.text <- 6
plt <- ggplot(pdfs, aes(x = x, color = short.name, fill = short.name)) +
  geom_path(aes(y = mean), linewidth = 0.1) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, color = NA) +
  geom_ribbon(aes(ymin = 0, ymax = low), alpha = 0.3, color = NA) +
  facet_grid2(gene + short.name ~ .,
              labeller = labeller(gene = c('eGFP-SV40pA' = "EGFP", 'Col2a1' = "Col2a1")),
              strip = strip_nested(text_y = c(list(element_text(size = size.text, face = "bold.italic"), element_text(size = size.text, face = "bold")),
                                              rep(list(element_text(angle = 0, face = "bold", size = size.text, hjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 1, unit = "mm"))), 16)),
                                   background_y = c(rep(list(element_part_rect(side = "r", fill = NA, colour = "black")), 2),
                                                    rep(list(element_rect(colour = NA, fill = NA)), 16)),
                                   size = "variable"),
              scales = "free", switch = "y") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title.x = element_text(size = size.text, face = "bold"),
        axis.text.x = element_text(size = size.text, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(.3, "mm"),
        strip.placement = "outside") +
  xlab("Expression level") + 
  scale_x_continuous(position = "top", breaks = c(0, 2, 4, 6)) +
  scale_y_continuous(breaks = 0) +
  scale_color_discrete("cluster", type = scales::hue_pal()(length(my.order))) +
  scale_fill_discrete("cluster", type = scales::hue_pal()(length(my.order))) +
  geom_path(aes(y = - mean), linewidth = 0.1) +
  geom_ribbon(aes(ymin = -low, ymax = -high), alpha = 0.5, color = NA) +
  geom_ribbon(aes(ymin = 0, ymax = -low), alpha = 0.3, color = NA) +
  coord_cartesian(ylim = c(-1.2, 1.2), xlim = c(0, 6), expand = F)

gt <- ggplot_gtable(ggplot_build(plt))
gt$heights[which(as.character(gt$heights) == "0.3mm")[8]] <- unit(1, "mm")
grid.draw(gt)
dir.create(dirname(output.prefix), showWarnings = FALSE)
ggsave(paste0(output.prefix, "_fixed_ylim_symetric.pdf"), gt, width = 66, height = 66, units = "mm")
