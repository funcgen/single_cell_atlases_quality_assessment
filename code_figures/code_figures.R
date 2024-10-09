### SCRIPT TO REPRODUCE THE FIGURES OF THE MANUSCRIPT ###

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DropletQC)
library(Matrix)
library(scCustomize)
library(viridis)
library(cowplot)
# Set color palette
pal <- viridis(n = 10, option = "D")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")


### 10x GENOMICS TABULA MURIS KIDNEY DATASET ###
# Load dataset
seurat_kidney <- readRDS("../datasets/kidney1_seurat.RDS")
seurat_kidney$lognCount <- log10(seurat_kidney@meta.data$nCount_RNA)

# Identify empty droplets
nf.umi4_5 <- data.frame(
  nf = seurat_kidney@meta.data$nuclear_fraction[seurat_kidney$orig.ident == "run4_5"], 
  umi = seurat_kidney@meta.data$nCount_RNA[seurat_kidney$orig.ident == "run4_5"]
)
nf.umi4_6 <- data.frame(
  nf = seurat_kidney@meta.data$nuclear_fraction[seurat_kidney$orig.ident == "run4_6"], 
  umi = seurat_kidney@meta.data$nCount_RNA[seurat_kidney$orig.ident == "run4_6"]
)
nf.umi7_5 <- data.frame(
  nf = seurat_kidney@meta.data$nuclear_fraction[seurat_kidney$orig.ident == "run7_5"], 
  umi = seurat_kidney@meta.data$nCount_RNA[seurat_kidney$orig.ident == "run7_5"]
)
res4_5 <- identify_empty_drops(nf_umi = nf.umi4_5, include_plot = TRUE)
res4_6 <- identify_empty_drops(nf_umi = nf.umi4_6, include_plot = TRUE)
res7_5 <- identify_empty_drops(nf_umi = nf.umi7_5, include_plot = TRUE)
res <- c(res4_5$cell_status, res4_6$cell_status, res7_5$cell_status)
seurat_kidney$cell_status <- res
seurat_kidney@meta.data$cell_status <- factor(seurat_kidney@meta.data$cell_status, levels = c("empty_droplet", "cell"))


# Introduce metadata variables
seurat_kidney$introns_fraction <- seurat_kidney@meta.data$nuclear_fraction
seurat_kidney$Malat1 <- seurat_kidney@assays$RNA@data["Malat1",]
seurat_kidney$Neat1 <- seurat_kidney@assays$RNA@data["Neat1",]
seurat_kidney$authors_annotations <- seurat_kidney$cell_ontology_class

# Take the metadata
metadata_kidney <- seurat_kidney@meta.data



### MACPARLAND LIVER DATASET ###
# Load dataset
seurat_liver <- readRDS(file = "../datasets/liver1_seurat.RDS")
seurat_liver$lognCount <- log10(seurat_liver@meta.data$nCount_RNA)

# Identify empty droplets
nf.umi1 <- data.frame(
  nf = seurat_liver@meta.data$nuclear_fraction[seurat_liver$orig.ident == "patient1"], 
  umi = seurat_liver@meta.data$nCount_RNA[seurat_liver$orig.ident == "patient1"]
)
nf.umi2 <- data.frame(
  nf = seurat_liver@meta.data$nuclear_fraction[seurat_liver$orig.ident == "patient2"], 
  umi = seurat_liver@meta.data$nCount_RNA[seurat_liver$orig.ident == "patient2"]
)
nf.umi3 <- data.frame(
  nf = seurat_liver@meta.data$nuclear_fraction[seurat_liver$orig.ident == "patient3"], 
  umi = seurat_liver@meta.data$nCount_RNA[seurat_liver$orig.ident == "patient3"]
)
nf.umi4 <- data.frame(
  nf = seurat_liver@meta.data$nuclear_fraction[seurat_liver$orig.ident == "patient4"], 
  umi = seurat_liver@meta.data$nCount_RNA[seurat_liver$orig.ident == "patient4"]
)
nf.umi5 <- data.frame(
  nf = seurat_liver@meta.data$nuclear_fraction[seurat_liver$orig.ident == "patient5"], 
  umi = seurat_liver@meta.data$nCount_RNA[seurat_liver$orig.ident == "patient5"]
)
res1 <- identify_empty_drops(nf_umi = nf.umi1, include_plot = TRUE, nf_rescue = 0.02)
res2 <- identify_empty_drops(nf_umi = nf.umi2, include_plot = TRUE)
res3 <- identify_empty_drops(nf_umi = nf.umi3, include_plot = TRUE)
res4 <- identify_empty_drops(nf_umi = nf.umi4, include_plot = TRUE)
res5 <- identify_empty_drops(nf_umi = nf.umi5, include_plot = TRUE, nf_rescue = 0.02)
res <- c(res1$cell_status, rep("cell", nrow(res2)), res3$cell_status, res4$cell_status, res5$cell_status)

seurat_liver$cell_status <- res
seurat_liver@meta.data$cell_status <- factor(seurat_liver@meta.data$cell_status, levels = c("empty_droplet", "cell"))


# Introduce metadata variables
seurat_liver$introns_fraction <- seurat_liver@meta.data$nuclear_fraction
seurat_liver$MALAT1 <- seurat_liver@assays$RNA@data["MALAT1",]
seurat_liver$NEAT1 <- seurat_liver@assays$RNA@data["NEAT1", ]
seurat_liver$authors_annotations <- seurat_liver$CellType

# Take the metadata
metadata_liver <- seurat_liver@meta.data



### HCA RETINA DATASET ###
# Load dataset
seurat_retina <- readRDS(file = "../datasets/retina1_seurat.RDS")
seurat_retina$lognCount <- log10(seurat_retina@meta.data$nCount_RNA)

# Identify empty droplets
nf.umi1 <- data.frame(
  nf = seurat_retina@meta.data$nuclear_fraction[seurat_retina$orig.ident == "sample1"], 
  umi = seurat_retina@meta.data$nCount_RNA[seurat_retina$orig.ident == "sample1"]
)
nf.umi2 <- data.frame(
  nf = seurat_retina@meta.data$nuclear_fraction[seurat_retina$orig.ident == "sample2"], 
  umi = seurat_retina@meta.data$nCount_RNA[seurat_retina$orig.ident == "sample2"]
)
nf.umi3 <- data.frame(
  nf = seurat_retina@meta.data$nuclear_fraction[seurat_retina$orig.ident == "sample3"], 
  umi = seurat_retina@meta.data$nCount_RNA[seurat_retina$orig.ident == "sample3"]
)
nf.umi4 <- data.frame(
  nf = seurat_retina@meta.data$nuclear_fraction[seurat_retina$orig.ident == "sample4"], 
  umi = seurat_retina@meta.data$nCount_RNA[seurat_retina$orig.ident == "sample4"]
)
nf.umi5 <- data.frame(
  nf = seurat_retina@meta.data$nuclear_fraction[seurat_retina$orig.ident == "sample5"], 
  umi = seurat_retina@meta.data$nCount_RNA[seurat_retina$orig.ident == "sample5"]
)
res1 <- identify_empty_drops(nf_umi = nf.umi1, include_plot = TRUE)
res2 <- identify_empty_drops(nf_umi = nf.umi2, include_plot = TRUE)
res3 <- identify_empty_drops(nf_umi = nf.umi3, include_plot = TRUE, nf_rescue = 0.07)
res4 <- identify_empty_drops(nf_umi = nf.umi4, include_plot = TRUE, nf_rescue = 0.07)
res5 <- identify_empty_drops(nf_umi = nf.umi5, include_plot = TRUE, nf_rescue = 0.06)
res <- c(res1$cell_status, res2$cell_status, res3$cell_status, res4$cell_status, res5$cell_status)

seurat_retina$cell_status <- res
seurat_retina@meta.data$cell_status <- factor(seurat_retina@meta.data$cell_status, levels = c("empty_droplet", "cell"))

# Introduce metadata variables
seurat_retina$introns_fraction <- seurat_retina@meta.data$nuclear_fraction
seurat_retina$MALAT1 <- seurat_retina@assays$RNA@data["MALAT1",]
seurat_retina$NEAT1 <- seurat_retina@assays$RNA@data["NEAT1", ]
seurat_retina$authors_annotations <- seurat_retina$cell.id.orig

# Take the metadata
metadata_retina <- seurat_retina@meta.data



### SMART-SEQ2 TABULA MURIS KIDNEY DATASET ###
# Load dataset
seurat_kidney_smart <- readRDS(file = "../datasets/kidney1_smart_seurat.RDS")
seurat_kidney_smart$lognCount <- log10(seurat_kidney_smart@meta.data$nCount_RNA)

# Introduce metadata variables
seurat_kidney_smart$nuclear_fraction <- seurat_kidney_smart$introns_percent / 100
seurat_kidney_smart$introns_fraction <- seurat_kidney_smart$nuclear_fraction
seurat_kidney_smart$Malat1 <- seurat_kidney_smart@assays$RNA@data["Malat1",]
seurat_kidney_smart$Neat1 <- seurat_kidney_smart@assays$RNA@data["Neat1",]
seurat_kidney_smart$authors_annotations <- seurat_kidney_smart$cell_ontology_class

# Identify empty droplets
nf.umi <- data.frame(nf = seurat_kidney_smart$nuclear_fraction, umi = seurat_kidney_smart$nCount_RNA)
res <- identify_empty_drops(nf_umi = nf.umi, include_plot = TRUE, nf_rescue = 0.07)
seurat_kidney_smart$cell_status <- res$cell_status
seurat_kidney_smart@meta.data$cell_status <- factor(seurat_kidney_smart@meta.data$cell_status, levels = c("empty_droplet", "cell"))

# Take the metadata
metadata_kidney_smart <- seurat_kidney_smart@meta.data




##### TABULA SAPIENS ######

### TABULA SAPIENS ALL ###
## TO RUN IN AN HPC ##
# Load dataset
seurat_sapiens_all <- readRDS("../datasets/sapiens_all_seurat.RDS")

# Introduce metadata variables
seurat_sapiens_all$MALAT1 <- seurat_sapiens_all@assays$RNA@data["ENSG00000251562",]
seurat_sapiens_all$authors_annotations <- seurat_sapiens_all$cell_type


### TABULA SAPIENS KIDNEY DATASET ###
# Load dataset
seurat_sapiens_kidney <- readRDS(file = "../datasets/sapiens_kidney_seurat.RDS")

# Introduce metadata variables
seurat_sapiens_kidney$MALAT1 <- seurat_sapiens_kidney@assays$RNA@data["ENSG00000251562", ]
metadata_sapiens_kidney <- seurat_sapiens_kidney@meta.data

# Add metadata
seurat_sapiens_kidney$authors_annotations <- seurat_sapiens_kidney$cell_type

### TABULA SAPIENS KIDNEY DATASET ###
# Load dataset
seurat_sapiens_kidney <- readRDS(file = "../datasets/sapiens_kidney_seurat.RDS")

# Introduce metadata variables
seurat_sapiens_kidney$MALAT1 <- seurat_sapiens_kidney@assays$RNA@data["ENSG00000251562", ]
metadata_sapiens_kidney <- seurat_sapiens_kidney@meta.data

# Add metadata
seurat_sapiens_kidney$authors_annotations <- seurat_sapiens_kidney$cell_type



### TABULA SAPIENS HEART DATASET ###
# Load dataset
seurat_sapiens_heart <- readRDS(file = "../datasets/sapiens_heart_seurat.RDS")

# Introduce metadata variables
seurat_sapiens_heart$MALAT1 <- seurat_sapiens_heart@assays$RNA@data["ENSG00000251562", ]
seurat_sapiens_heart$authors_annotations <- seurat_sapiens_heart$cell_type

# Take the metadata
metadata_sapiens_heart <- seurat_sapiens_heart@meta.data


### TABULA SAPIENS PROSTATE ###
# Load dataset
seurat_sapiens_prostate <- readRDS(file = "../datasets/sapiens_prostate_seurat.RDS")

# Introduce metadata variables
seurat_sapiens_prostate$MALAT1 <- seurat_sapiens_prostate@assays$RNA@data["ENSG00000251562", ]
seurat_sapiens_prostate$authors_annotations <- seurat_sapiens_prostate$cell_type

# Take the metadata
metadata_sapiens_prostate <- seurat_sapiens_prostate@meta.data





##### TABULA MURIS SENIS #####

### TABULA MURIS SENIS ALL ###
## TO RUN IN AN HPC ##

# Load dataset
seurat_muris_all <- readRDS("../../TABULA_MURIS_SENIS/data/tabula_muris_senis.rds")

# Introduce metadata variables
seurat_muris_all$Malat1 <- seurat_muris_all@assays$RNA@data["ENSG00000251562",]
seurat_muris_all$authors_annotations <- seurat_muris_all$cell_type


### TABULA MURIS SENIS PANCREAS 10x ### 
seurat_muris_pancreas_10x <- readRDS("../data/seurat_muris_pancreas_10x.rds")
seurat_muris_pancreas_10x$Malat1 <- seurat_muris_pancreas_10x@assays$RNA@data["ENSMUSG00000092341",]


### TABULA MURIS SENIS KIDNEY 10x ###
seurat_muris_kidney_10x <- readRDS("../../TABULA_MURIS_SENIS/data/seurat_muris_kidney_10x.rds")
seurat_muris_kidney_10x$Malat1 <- seurat_muris_kidney_10x@assays$RNA@data["ENSMUSG00000092341",]


### TABULA MURIS SENIS LIVER 10x ###
seurat_muris_liver_10x <- readRDS("../../TABULA_MURIS_SENIS/data/seurat_muris_liver_10x.rds")
seurat_muris_liver_10x$Malat1 <- seurat_muris_liver_10x@assays$RNA@data["ENSMUSG00000092341",]



### TABULA MURIS SENIS PANCREAS smart ### 
seurat_muris_pancreas_smart <- readRDS("../../TABULA_MURIS_SENIS/data/seurat_muris_pancreas_smart.rds")
seurat_muris_pancreas_smart$Malat1 <- seurat_muris_pancreas_smart@assays$RNA@data["ENSMUSG00000092341",]


### TABULA MURIS SENIS KIDNEY smart ###
seurat_muris_kidney_smart <- readRDS("../../TABULA_MURIS_SENIS/data/seurat_muris_kidney_smart.rds")
seurat_muris_kidney_smart$Malat1 <- seurat_muris_kidney_smart@assays$RNA@data["ENSMUSG00000092341",]


### TABULA MURIS SENIS LIVER smart ###
seurat_muris_liver_smart <- readRDS("../../TABULA_MURIS_SENIS/data/seurat_muris_liver_smart.rds")
seurat_muris_liver_smart$Malat1 <- seurat_muris_liver_smart@assays$RNA@data["ENSMUSG00000092341",]














### CREATING A COMMON DATA.FRAME ###

# Select the variable of interest to include in the common metadata data.frame
metadata_kidney_simpler <- metadata_kidney[, c("cell_status", "Malat1", "Neat1", "introns_fraction", "authors_annotations")]
metadata_kidney_simpler$tissue <- "kidney_10x_TM"
metadata_liver_simpler <- metadata_liver[, c("cell_status", "MALAT1", "NEAT1", "introns_fraction", "authors_annotations")]
colnames(metadata_liver_simpler) <- c("cell_status", "Malat1", "Neat1", "introns_fraction", "authors_annotations")
metadata_liver_simpler$tissue <- "liver"
metadata_retina_simpler <- metadata_retina[, c("cell_status", "MALAT1", "NEAT1", "introns_fraction", "authors_annotations")]
colnames(metadata_retina_simpler) <- c("cell_status", "Malat1", "Neat1", "introns_fraction", "authors_annotations")
metadata_retina_simpler$tissue <- "retina"
metadata_kidney_smart_simpler <- metadata_kidney_smart[, c("cell_status", "Malat1", "Neat1", "introns_fraction", "authors_annotations")]
colnames(metadata_kidney_smart_simpler) <- c("cell_status", "Malat1", "Neat1", "introns_fraction", "authors_annotations")
metadata_kidney_smart_simpler$tissue <- "kidney_smart_seq2_TM"

# Create the common data.frame
metadata_all <- list(metadata_kidney_simpler, metadata_liver_simpler, metadata_retina_simpler, metadata_kidney_smart_simpler)
metadata_df <- do.call("rbind", metadata_all)


### TABLE WITH CORRELATIONS BETWEEN MALAT1 AND INTRONS FRACTION ###
cor_vector <- c()
for (i in 1:length(unique(metadata_df$tissue))) {
  cor_vector <- c(cor_vector, cor(
    metadata_df$Malat1[metadata_df$tissue == unique(metadata_df$tissue)[i]],
    metadata_df$introns_fraction[metadata_df$tissue == unique(metadata_df$tissue)[i]],
    method = "pearson"
  ))
}
cor_df <- data.frame(
  dataset = c("kidney_10x", "liver", "retina", "kidney_smart_seq2"),
  cor = cor_vector
)
cor_df






### KIDNEY2 DATASET ###
# Load dataset
seurat_kidney2 <- readRDS(file = "../malat1_datasets/kidney.RDS")


# Introduce metadata variables
seurat_kidney2$Malat1 <- seurat_kidney2@assays$RNA@data["ENSMUSG00000092341",]
seurat_kidney2$Neat1 <- seurat_kidney2@assays$RNA@data["ENSMUSG00000092274", ]
seurat_kidney2$authors_annotations <- seurat_kidney2$author_cell_type

# Add Malat1_status variable
seurat_kidney2$Malat1_status <- NA
seurat_kidney2$Malat1_status <- ifelse(
  seurat_kidney2$Malat1 < 3.5,
  "Cytosolic debris (sc)", "Cells (sc)"
)
seurat_kidney2$Malat1_status <- ifelse(
  seurat_kidney2$suspension_type == "nucleus", 
  "Single-nucleus (sn)", seurat_kidney2$Malat1_status
)


seurat_kidney2$Malat1_status <- factor(seurat_kidney2$Malat1_status, levels = c("Cytosolic debris (sc)", "Cells (sc)", "Single-nucleus (sn)"))



# Take metadata
metadata_kidney2 <- seurat_kidney2@meta.data

# Split by type of suspension (sc and sn)
seurat_kidney2_sc <- subset(seurat_kidney2, subset = suspension_type == "cell")
seurat_kidney2_sn <- subset(seurat_kidney2, subset = suspension_type == "nucleus")
metadata_kidney2_sc <- seurat_kidney2_sc@meta.data
metadata_kidney2_sn <- seurat_kidney2_sn@meta.data















################
### FIGURE 1 ###
################

# Figure 1a
plot1 <- DimPlot(
  object = seurat_kidney,
  group.by = "authors_annotations",
  repel = TRUE,
  label = TRUE,
  label.size = 4,
) + NoLegend() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Original authors annotations")

plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney, 
  features = c("introns_fraction"),
  max.cutoff = c(
    quantile(seurat_kidney@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal
) + 
  theme(
    legend.position = c(0.03, 0.9),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")

plot3 <- DimPlot(
  object = seurat_kidney,
  group.by = "cell_status"
) + 
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Barcode status")

plota <- ggpubr::ggarrange(plot1, plot2, plot3, ncol=3, nrow=1, common.legend = FALSE)
plota <- ggpubr::annotate_figure(plota, top = ggpubr::text_grob("Tabula Muris kidney 10x", color = "black", face = "bold", 
                                                                size = 25
                                                                )
                                      )

# Figure 1b
plot1 <- DimPlot(
  object = seurat_liver,
  group.by = "authors_annotations",
  repel = TRUE,
  label = TRUE,
  label.size = 4,
) + NoLegend() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Original authors annotations")

plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_liver, 
  features = c("introns_fraction"),
  max.cutoff = c(
    quantile(seurat_liver@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal
) +
  theme(
    legend.position = c(0.05, 0.9),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") + 
  labs(title = "Introns fraction")
plot3 <- DimPlot(
  object = seurat_liver,
  group.by = "cell_status",
  repel = TRUE
) + 
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Barcode status")
plotb <- ggpubr::ggarrange(plot1, plot2, plot3, ncol=3, nrow=1, common.legend = FALSE)
plotb <- ggpubr::annotate_figure(plotb, top = ggpubr::text_grob("MacParland liver", color = "black", face = "bold", 
                                                                size = 25
                                                                )
)

# Figure 1c
plot1 <- DimPlot(
  object = seurat_retina,
  group.by = "authors_annotations",
  repel = TRUE,
  label = TRUE,
  label.size = 4,
) + NoLegend() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
        ) +
  labs(title = "Original authors annotations")

plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_retina, 
  features = c("introns_fraction"),
  max.cutoff = c(
    quantile(seurat_retina@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal
) + 
  theme(
    legend.position = c(0.05, 0.9),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plot3 <- DimPlot(
  object = seurat_retina,
  group.by = "cell_status"
) + 
  theme(
    legend.position = c(0.05, 0.96),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Barcode status")
plotc <- ggpubr::ggarrange(plot1, plot2, plot3, ncol=3, nrow=1, common.legend = FALSE)
plotc <- ggpubr::annotate_figure(plotc, top = ggpubr::text_grob("Lukowski retina", color = "black", face = "bold", 
                                                                size = 25
                                                                )
)

# Figure 1d
plot1 <- DimPlot(
  object = seurat_kidney_smart,
  group.by = "authors_annotations",
  repel = TRUE,
  label = TRUE,
  label.size = 4,
) + NoLegend() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Original authors annotations")

plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney_smart, 
  features = c("introns_fraction"),
  max.cutoff = c(
    quantile(seurat_kidney_smart@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal
) + 
  theme(
      legend.position = c(0.05, 0.9),
      axis.ticks = element_blank(),
      axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plot3 <- DimPlot(
  object = seurat_kidney_smart,
  group.by = "cell_status"
) + 
  theme(
      legend.position = c(0.05, 0.95),
      axis.ticks = element_blank(),
      axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Barcode status")
plotd <- ggpubr::ggarrange(plot1, plot2, plot3, ncol=3, nrow=1, common.legend = FALSE)
plotd <- ggpubr::annotate_figure(plotd, top = ggpubr::text_grob("Tabula Muris Smart-seq2 kidney", color = "black", face = "bold", 
                                                                size = 25
                                                                ))

plot_fig1 <- ggpubr::ggarrange(plota, plotb, plotc, plotd, ncol=1, nrow=4, common.legend = TRUE, labels = "AUTO")

ggsave(
  filename = "../results/Fig1.pdf",
  plot = plot_fig1,
  width = 17,
  height = 20
)
ggsave(
  filename = "../results/Fig1.jpeg",
  plot = plot_fig1,
  width = 17,
  height = 20,
  bg="white"
)












################
### FIGURE 2 ###
################

# Load necessary libraries
library(ggplot2)
library(dplyr)

metadata_df_kidney_10x = subset(metadata_df, subset = tissue == "kidney_10x_TM")
plot1 = ggplot(metadata_df_kidney_10x, aes(x = authors_annotations, y = introns_fraction)) +
  geom_violin(fill = "white", color = "black") +
  # geom_boxplot(width = 0.1, fill = "white") +  # Boxplot inside violin
  geom_jitter(width = 0.1, size = 0.5, aes(color = cell_status)) +  # Add jittered points
  labs(title = "Tabula Muris kidney 10x", 
       x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    legend.title=element_blank(),
    plot.title = element_text(size = 18),     # Title size
    axis.title.x = element_text(size = 15),   # X-axis title size
    axis.title.y = element_text(size = 15),   # Y-axis title size
    axis.text.y = element_text(size = 12),    # Y-axis text (tick labels) size
    legend.text = element_text(size = 15) 
  )


metadata_df_liver = subset(metadata_df, subset = tissue == "liver")
plot2 = ggplot(metadata_df_liver, aes(x = authors_annotations, y = introns_fraction)) +
  geom_violin(fill = "white", color = "black") +
  # geom_boxplot(width = 0.1, fill = "white") +  # Boxplot inside violin
  geom_jitter(width = 0.1, size = 0.5, aes(color = cell_status)) +  # Add jittered points
  labs(title = "MacParland liver", 
       x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    legend.title=element_blank(),
    plot.title = element_text(size = 18),     # Title size
    axis.title.x = element_text(size = 15),   # X-axis title size
    axis.title.y = element_text(size = 15),   # Y-axis title size
    axis.text.y = element_text(size = 12),    # Y-axis text (tick labels) size
    legend.text = element_text(size = 15) 
  )

metadata_df_retina = subset(metadata_df, subset = tissue == "retina")
plot3 = ggplot(metadata_df_retina, aes(x = authors_annotations, y = introns_fraction)) +
  geom_violin(fill = "white", color = "black") +
  # geom_boxplot(width = 0.1, fill = "white") +  # Boxplot inside violin
  geom_jitter(width = 0.1, size = 0.5, aes(color = cell_status)) +  # Add jittered points
  labs(title = "Lukowski retina", 
       x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    legend.title=element_blank(),
    plot.title = element_text(size = 18),     # Title size
    axis.title.x = element_text(size = 15),   # X-axis title size
    axis.title.y = element_text(size = 15),   # Y-axis title size
    axis.text.y = element_text(size = 12),    # Y-axis text (tick labels) size
    legend.text = element_text(size = 15) 
  )

metadata_df_kidney_smart_seq2 = subset(metadata_df, subset = tissue == "kidney_smart_seq2_TM")
plot4 = ggplot(metadata_df_kidney_smart_seq2, aes(x = authors_annotations, y = introns_fraction)) +
  geom_violin(fill = "white", color = "black") +
  # geom_boxplot(width = 0.1, fill = "white") +  # Boxplot inside violin
  geom_jitter(width = 0.1, size = 0.5, aes(color = cell_status)) +  # Add jittered points
  labs(title = "Tabula Muris Smart-seq2 kidney", 
       x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    legend.title=element_blank(),
    plot.title = element_text(size = 18),     # Title size
    axis.title.x = element_text(size = 15),   # X-axis title size
    axis.title.y = element_text(size = 15),   # Y-axis title size
    axis.text.y = element_text(size = 12),    # Y-axis text (tick labels) size
    legend.text = element_text(size = 15) 
  )


plota <- ggpubr::ggarrange(
  plot1, plot2, plot3, plot4, 
  ncol=2, nrow=2, 
  align = "h",
  common.legend = TRUE,
  legend = "right",
  labels = "AUTO"
  )
plota <- ggpubr::annotate_figure(plota, top = ggpubr::text_grob(
  "Introns-free cells distribution", color = "black", face = "bold", size = 25
  ), 
  left = ggpubr::text_grob("Intronic fraction", rot = 90)
)

plota

ggsave(
  filename = "../results/Fig2.pdf",
  plot = plota,
  width = 17,
  height = 20
)

ggsave(
  filename = "../results/Fig2.jpeg",
  plot = plota,
  width = 17,
  height = 20,
  bg = "white"
)









################
### FIGURE 3 ###
################

# Figure 3a
plota1 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_kidney@assays$RNA@data["Malat1",], 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Malat1")

plota2 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney,
  features = "introns_fraction",
  max.cutoff = c(
    quantile(seurat_kidney@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plota <- ggpubr::ggarrange(plota1, plota2, ncol=2, nrow=1, common.legend = FALSE)
plota <- ggpubr::annotate_figure(plota, top = ggpubr::text_grob("Tabula Muris kidney 10x", color = "black", face = "bold", 
                                                                size = 25
))

# Figure 3b
plotb1 <- FeaturePlot_scCustom(
  seurat_object = seurat_liver,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_liver@assays$RNA@data["MALAT1",], 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Malat1")

plotb2 <- FeaturePlot_scCustom(
  seurat_object = seurat_liver,
  features = "introns_fraction",
  max.cutoff = c(
    quantile(seurat_liver@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plotb <- ggpubr::ggarrange(plotb1, plotb2, ncol=2, nrow=1, common.legend = FALSE)
plotb <- ggpubr::annotate_figure(plotb, top = ggpubr::text_grob("MacParland liver", color = "black", face = "bold", 
                                                                size = 25
))


# Figure 3c
plotc1 <- FeaturePlot_scCustom(
  seurat_object = seurat_retina,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_retina@assays$RNA@data["MALAT1",], 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Malat1")

plotc2 <- FeaturePlot_scCustom(
  seurat_object = seurat_retina,
  features = "introns_fraction",
  max.cutoff = c(
    quantile(seurat_retina@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plotc <- ggpubr::ggarrange(plotc1, plotc2, ncol=2, nrow=1, common.legend = FALSE)
plotc <- ggpubr::annotate_figure(plotc, top = ggpubr::text_grob("Lukowski retina", color = "black", face = "bold", 
                                                                size = 25
))


# Figure 3d
plotd1 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney_smart,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_kidney_smart@assays$RNA@data["Malat1",], 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Malat1")

plotd2 <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney_smart,
  features = "introns_fraction",
  max.cutoff = c(
    quantile(seurat_kidney_smart@meta.data$nuclear_fraction, 0.95)
  ),
  colors_use = pal,
  na_cutoff = c(0)
) +
  theme(
    legend.position = c(0.05, 0.95),
    axis.ticks = element_blank(),
    axis.text = element_blank()
    # axis.title.x = element_blank()
  ) +
  xlab("") +
  ylab("") +
  labs(title = "Introns fraction")
plotd <- ggpubr::ggarrange(plotd1, plotd2, ncol=2, nrow=1, common.legend = FALSE)
plotd <- ggpubr::annotate_figure(plotd, top = ggpubr::text_grob("Tabula Muris Smart-seq2 kidney", color = "black", face = "bold", 
                                                                size = 25
))




# Figure 3e
plote <- ggplot(
  data = metadata_df, 
  mapping = aes(x = cell_status, y = Malat1, fill = cell_status)
) +
  geom_boxplot() +
  xlab("DropletQC barcode classification") +
  ylab("Malat1 normalized counts")  +
  facet_grid(. ~ tissue) +
  labs(fill = "Barcode status") +
  theme_bw() + 
  theme(legend.position = "none")
plot_fig3 <- ggpubr::ggarrange(plota, plotb, plotc, plotd, plote, ncol=1, nrow=5, common.legend = FALSE, labels = "AUTO")

plot_fig3

ggsave(
  filename = "../results/Fig3.pdf",
  plot = plot_fig3,
  width = 15,
  height = 22
)


ggsave(
  filename = "../results/Fig3.jpeg",
  plot = plot_fig3,
  width = 15,
  height = 22,
  bg = "white"
)





################
### FIGURE 4 ###
################

# Re-annotated cells
metadata_kidney2 <- seurat_kidney2@meta.data

# Create new column
metadata_kidney2$reduced_annotation <- NA 
metadata_kidney2$reduced_annotation <- ifelse(
  metadata_kidney2$authors_annotations %in% 
    c(
      "PTS1", "PTS2", "PTS3", "PTS3T2",
      "MTAL", 
      "CTAL", 
      "ICB",
      "DCT", 
      "CNT", 
      "PC"
    ),
  as.character(metadata_kidney2$authors_annotations), "Other cell type"
)

seurat_kidney2@meta.data <- metadata_kidney2



# Figure 4a
plota <- DimPlot(
  object = seurat_kidney2,
  group.by = "reduced_annotation",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
  # cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "darkgrey", "#F564E3")
) + NoLegend() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(
    title = "Authors annotations"
    )

# Figure 4b
plotb <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney2_sc,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_kidney2@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  # num_columns = 3,
  colors_use = pal,
  na_cutoff = 0
) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
  ) +
  labs(title = "Malat1 | single-cells") +
  xlab("") +
  ylab("")

# Figure 4c
plotc <- FeaturePlot_scCustom(
  seurat_object = seurat_kidney2_sn,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_kidney2@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  # num_columns = 3,
  colors_use = pal,
  na_cutoff = 0
) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Malat1 | single-nuclei") +
  xlab("") +
  ylab("")

plot_fig4 <- ggpubr::ggarrange(plota, plotb, plotc, ncol=3, nrow=1, common.legend = FALSE, labels = "AUTO")
plot <- ggpubr::annotate_figure(
  plot_fig4, top = ggpubr::text_grob("Kidney Novella-Rausell", color = "black", face = "bold", size = 20)
)
plot


ggsave(
  filename = "../results/Fig4.pdf",
  plot = plot,
  width = 15,
  height = 8
)


ggsave(
  filename = "../results/Fig4.jpeg",
  plot = plot,
  width = 15,
  height = 8,
  bg = "white"
)





################
### FIGURE 5 ###
################

metadata_kidney2$mitochondrial_ratio <- metadata_kidney2$pct_counts_mt/100
metadata_kidney2_sc$mitochondrial_ratio <- metadata_kidney2_sc$pct_counts_mt/100
metadata_kidney2_sn$mitochondrial_ratio <- metadata_kidney2_sn$pct_counts_mt/100

# Figure 5a
plot5a <- ggplot(
  data = metadata_kidney2, mapping = aes(x = Malat1)
) +
  geom_histogram(bins = 60, fill = "steelblue", color = "black") +
  facet_wrap(~ suspension_type, labeller = labeller(
    suspension_type = c("cell" = "single-cell", "nucleus" = "single-nucleus")
    )
  ) +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
  geom_vline(xintercept = 3.5) +
  xlab("Malat1 expression") +
  ylab("Number of cells") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 15, l = 0)),
    strip.text = element_text(
      size = 20)
  )

# Figure 5b
plot5b <- ggplot(data = metadata_kidney2, mapping = aes(x = Malat1_status, y = mitochondrial_ratio, fill = Malat1_status)) +
  geom_boxplot(show.legend = FALSE) +
  xlab("") +
  ylab("Mitochondrial ratio") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
    strip.text = element_text(
      size = 20
    )
  ) +
  ylim(c(0, 1))

plot_fig5 <- ggpubr::ggarrange(plot5a, plot5b, ncol=1, nrow=2, common.legend = FALSE, labels = "AUTO")
plot <- ggpubr::annotate_figure(
  plot_fig5, top = ggpubr::text_grob("Kidney Novella-Rausell", color = "black", face = "bold", size = 20)
)

plot

ggsave(
  filename = "../results/Fig5.pdf",
  plot = plot,
  width = 15,
  height = 15
)

ggsave(
  filename = "../results/Fig5.jpeg",
  plot = plot,
  width = 15,
  height = 15,
  bg = "white"
)


################
### FIGURE 6 ###
################

plot6a <- ggplot(data = metadata_kidney, mapping = aes(x = introns_fraction, y = log10(nCount_RNA), color = percent.mt)) +
  geom_point(size = 2) +
  ylab("log10(nUMI)") +
  xlab("Intron fraction") +
  xlim(c(0, 0.6)) +
  labs(title = "Tabula Muris 10x kidney", color = "% Mitochondrial") +
  theme_bw() +
  scale_color_gradient(low="blue", high="red") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 20),
    legend.title=element_text(size=15)
    
  )
plot6b <- ggplot(data = metadata_liver, mapping = aes(x = introns_fraction, y = log10(nCount_RNA), color = percent.mt)) +
  geom_point(size = 2) + 
  ylab("log10(nUMI") +
  xlab("Intron fraction") +
  labs(title = "MacParland 10x liver", color = "% Mitochondrial") +
  theme_bw() +
  scale_color_gradient(low="blue", high="red") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 20),
    legend.title=element_text(size=15)
    
  )
plot6c <- ggplot(data = metadata_retina, mapping = aes(x = introns_fraction, y = log10(nCount_RNA), color = percent.mt)) +
  geom_point(size = 2) +
  ylab("log10(nUMI") +
  xlab("Intron fraction") +
  labs(title = "Lukowski 10x retina", color = "% Mitochondrial") +
  theme_bw() +
  scale_color_gradient(low="blue", high="red") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 20),
    legend.title=element_text(size=15)
    
  )
plot6d <- ggplot(data = metadata_kidney_smart, mapping = aes(x = introns_fraction, y = log10(nCount_RNA), color = percent.mt)) +
  geom_point(size = 2) +
  ylab("log10(nUMI") +
  xlab("Intron fraction") +
  labs(title = "Tabula Muris Smart-seq2 kidney", color = "% Mitochondrial") +
  theme_bw() +
  scale_color_gradient(low="blue", high="red") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 20),
    legend.title=element_text(size=15)
    
  )

plot_fig6 <- ggpubr::ggarrange(
  plot6a, plot6b, plot6c, plot6d, ncol=2, nrow=2, common.legend = TRUE, labels = "AUTO", legend = "bottom"
  )
plot_fig6


ggsave(
  filename = "../results/Fig6.pdf",
  plot = plot_fig6,
  width = 25,
  height = 15
)

ggsave(
  filename = "../results/Fig6.jpeg",
  plot = plot_fig6,
  width = 25,
  height = 15,
  bg = "white"
)




################
### FIGURE 7 ###
################


### TO RUN IN AN HPC ###

# Figure 7a
plot1 <- DimPlot(
  object = seurat_sapiens_all,
  group.by = "tissue_in_publication",
  label = FALSE,
  reduction = "umap",
  raster = FALSE
) +
labs(title = "Tissue in publication") +
theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
)
plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_sapiens_all,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_sapiens_all@assays$RNA@data["ENSG00000251562",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  raster = FALSE,
  na_cutoff = 0
) +
labs(title = "MALAT1") +
theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank()
)
plot7a <- ggpubr::ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = FALSE)
plot7a <- ggpubr::annotate_figure(plot7a, top = ggpubr::text_grob("Tabula Sapiens | All", color = "black", face = "bold", size = 20))



# Figure 7b 
plot1 <- DimPlot(
  object = seurat_sapiens_kidney,
  group.by = "authors_annotations",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
)
plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_sapiens_kidney,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_sapiens_kidney@assays$RNA@data["ENSG00000251562",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
) +
xlab("") +
ylab("")
plot7b <- ggpubr::ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = FALSE)
plot7b <- ggpubr::annotate_figure(plot7b, top = ggpubr::text_grob("Tabula Sapiens | Kidney", color = "black", face = "bold", size = 20))



# Figure 7c
plot1 <- DimPlot(
  object = seurat_sapiens_heart,
  group.by = "authors_annotations",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations") + 
    theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
)
plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_sapiens_heart,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_sapiens_heart@assays$RNA@data["ENSG00000251562",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
) +
    theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
) +
xlab("") + 
ylab("")
plot7c <- ggpubr::ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = FALSE)
plot7c <- ggpubr::annotate_figure(plot7c, top = ggpubr::text_grob("Tabula Sapiens | Heart", color = "black", face = "bold", size = 20))



# Figure 7d
plot1 <- DimPlot(
  object = seurat_sapiens_prostate,
  group.by = "authors_annotations",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
)
plot2 <- FeaturePlot_scCustom(
  seurat_object = seurat_sapiens_prostate,
  features = "MALAT1",
  max.cutoff = c(
    quantile(seurat_sapiens_prostate@assays$RNA@data["ENSG00000251562",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
) +
xlab("") + 
ylab("")
plot6d <- ggpubr::ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = FALSE)
plot6d <- ggpubr::annotate_figure(plot7d, top = ggpubr::text_grob("Tabula Sapiens | Prostate", color = "black", face = "bold", size = 20))

plot_fig7 <- ggpubr::ggarrange(plot7a, plot7b, plot7c, plot7d, ncol=1, nrow=4, common.legend = FALSE, labels = "AUTO")


ggsave(
  filename = "../results/Fig7.png",
  plot = plot_fig7,
  width = 15,
  height = 20
)







### TO RUN IN AN HPC ###

######################
### FIGURE SUPPL 1 ###
######################

# Figure supp 1a
plot_supp1 <- DimPlot(
  object = seurat_muris_all,
  group.by = "tissue",
  label = FALSE,
  reduction = "umap",
  raster = FALSE
) +
labs(title = "Tissue in publication")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_all,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_all@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  raster = FALSE,
  na_cutoff = 0
) +
labs(title = "Malat1")
plot1a <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot1a <- ggpubr::annotate_figure(plot1a, top = ggpubr::text_grob("Tabula Muris | All", color = "black", face = "bold", size = 14))


# Figure supp 1b 10x
plot_supp1 <- DimPlot(
  object = seurat_muris_pancreas_10x,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_pancreas_10x,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_pancreas_10x@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot1b <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot1b <- ggpubr::annotate_figure(plot1b, top = ggpubr::text_grob("Tabula Muris | Pancreas", color = "black", face = "bold", size = 14))


ggsave(
  filename = "../results/fig_supp1b.png",
  plot = plot7b,
  width = 20,
  height = 15
)

# Figure supp 1c 10x
plot_supp1 <- DimPlot(
  object = seurat_muris_kidney_10x,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_kidney_10x,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_kidney_10x@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot1c <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot1c <- ggpubr::annotate_figure(plot1c, top = ggpubr::text_grob("Tabula Muris | Kidney", color = "black", face = "bold", size = 14))


# Figure supp 1d 10x
plot_supp1 <- DimPlot(
  object = seurat_muris_liver_10x,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_liver_10x,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_liver_10x@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot1d <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot1d <- ggpubr::annotate_figure(plot1d, top = ggpubr::text_grob("Tabula Muris | Liver", color = "black", face = "bold", size = 14))

plot_fig_supp1 <- ggpubr::ggarrange(plot1a, plot1b, plot1c, plot1d, ncol=1, nrow=4, common.legend = FALSE, labels = "AUTO")




######################
### FIGURE SUPPL 2 ###
######################

# Figure supp 1b smart-seq2 
plot_supp1 <- DimPlot(
  object = seurat_muris_pancreas_smart,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_pancreas_smart,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_pancreas_smart@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot2b <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot2b <- ggpubr::annotate_figure(plot2b, top = ggpubr::text_grob("Tabula Muris | Pancreas", color = "black", face = "bold", size = 14))


# Figure supp 1c smart-seq2
plot_supp1 <- DimPlot(
  object = seurat_muris_kidney_smart,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_kidney_smart,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_kidney_smart@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot2c <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot2c <- ggpubr::annotate_figure(plot2c, top = ggpubr::text_grob("Tabula Muris | Kidney", color = "black", face = "bold", size = 14))


# Figure supp 1d smart-seq2
plot_supp1 <- DimPlot(
  object = seurat_muris_liver_smart,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + NoLegend() +
  labs(title = "Authors annotations")
plot_supp2 <- FeaturePlot_scCustom(
  seurat_object = seurat_muris_liver_smart,
  features = "Malat1",
  max.cutoff = c(
    quantile(seurat_muris_liver_smart@assays$RNA@data["ENSMUSG00000092341",], 0.95)
  ),
  colors_use = pal,
  reduction = "umap",
  na_cutoff = 0
)
plot2d <- ggpubr::ggarrange(plot_supp1, plot_supp2, ncol=2, nrow=1, common.legend = FALSE)
plot2d <- ggpubr::annotate_figure(plot2d, top = ggpubr::text_grob("Tabula Muris | Liver", color = "black", face = "bold", size = 14))


plot_fig7 <- ggpubr::ggarrange(plot2a, plot2b, plot2c, plot2d, ncol=1, nrow=4, common.legend = FALSE, labels = "AUTO")

ggsave(
  filename = "../results/suppl_Fig1_smart.png",
  plot = plot_fig7,
  width = 20,
  height = 30
)
