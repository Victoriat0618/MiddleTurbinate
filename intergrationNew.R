library(Seurat)
library(tidyverse)
library(sleepwalk)
library(SCINA)
library(BiocManager)
library(DESeq2)

#importdata
SEPscL310_counts <-Read10X(data.dir ="/home/bzawisza/seurat/SEPscL310/filtered_feature_bc_matrix")

#loading Seurat Objects
SEPscL310_seuratobj<-CreateSeuratObject(counts = SEPscL310_counts, project = "SEPscL310", min.cells = 3, min.features = 200)
SEPscL310_seuratobj@meta.data$orig.ident <- as.factor("SEPscL310")

##QC##
#nCount_RNA the total number of reads in the dataset
#nFeature_RNA the number of observed genes
#number of genes from mitochondrial sources

SEPscL310_seuratobj[["percent.mt"]] <- PercentageFeatureSet(SEPscL310_seuratobj, pattern = "^MT-")

##calculate percentage of largest gene
#calculate what percentage of the data comes from the single most observed gene
#having a high proportion of your data dominated by a single gene would be a concerning sign
SEPscL310_seuratobj$Percent.Largest.Gene <- apply(
  SEPscL310_seuratobj@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
) 

head(SEPscL310_seuratobj$Percent.Largest.Gene)

###Normalisation, Selection and Scaling
##Normalize USING log normalize in Seurat
##"SCT","CLR"
SEPscL310_seuratobj <- NormalizeData(SEPscL310_seuratobj, normalization.method = "CLR")

##get a list of the most highly expressed genes
##identify any potential nuclear expressed transcript which tends to persist when cells have lysed and the cytoplasm has gone
##or any high amounts of ribosomal proteins
gene.expression<- apply(SEPscL310_seuratobj@assays$RNA@data,1,mean)

gene.expression<- sort(gene.expression, decreasing = TRUE) 

head(gene.expression, n=50)

##look at how well the data have been normalised by looking at any housekeeping genes
ggplot(mapping = aes(SEPscL310_seuratobj@assays$RNA@data["GAPDH",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("GAPDH expression")

##look a bit wider and pick the first 100 cells and look at the distributions of their expression values
as_tibble(
  SEPscL310_seuratobj@assays$RNA@data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,6), xlim=c(0,3))

##overall metrics after normalization
##compare the quantitative value at the 95th percentile to the mean expression
normalisation.qc<- tibble(
  pc95 = apply(SEPscL310_seuratobj[["RNA"]]@data,2,quantile,0.95),
  measured = apply(SEPscL310_seuratobj[["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
)

normalisation.qc %>% 
  ggplot(aes(x=measured,y=pc95))+
  geom_point()+
  ggtitle("Normalisation of data")


###Cell Cycle Scoring
SEPscL310_seuratobj <- CellCycleScoring(SEPscL310_seuratobj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
SEPscL310_seuratobj[[]]

as_tibble(SEPscL310_seuratobj[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

as_tibble(SEPscL310_seuratobj[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

###Gene selection
##selecting the 500 most variable genes
##calculate a normalised intensity for each gene, and can then select the top ‘n’ most variable features
SEPscL310_seuratobj<- FindVariableFeatures(
  SEPscL310_seuratobj, 
  selection.method = "vst", 
  nfeatures=500
)

##access to variability information
variance.data<- as_tibble(HVFInfo(SEPscL310_seuratobj),rownames = "Gene") 

variance.data %>% 
  mutate(hypervariable=Gene %in% VariableFeatures(SEPscL310_seuratobj)
  ) -> variance.data

head(variance.data, n=10)

##graph of the variance vs mean and highlight the selected genes so we can see whether we think we’re likely to capture what we need
variance.data %>% 
  ggplot(aes(log(mean),log(variance),color=hypervariable)) + 
  geom_point() + 
  scale_color_manual(values=c("black","red"))

###Scaling
SEPscL310_seuratobj<- ScaleData(SEPscL310_seuratobj,features=rownames(SEPscL310_seuratobj))


###Dimensionality Reduction
##PCA
SEPscL310_seuratobj <- RunPCA(SEPscL310_seuratobj,features=VariableFeatures(SEPscL310_seuratobj))


DimPlot(SEPscL310_seuratobj,reduction="pca")
##looking at different PC
DimPlot(SEPscL310_seuratobj,reduction="pca", dims=c(9,10))

##elbow plot which simply quantitates the amount of variance captured in the different PCs
ElbowPlot(SEPscL310_seuratobj)

###Defining Cell Clusters
##finds the ‘k’ nearest neighbours to each cell and makes this into a graph. 
##It then looks for highly inter-connected subgraphs within the graph and uses these to define clusters.
SEPscL310_seuratobj <- FindNeighbors(SEPscL310_seuratobj,dims=1:15)
SEPscL310_seuratobj@graphs$RNA_snn[1:10,1:10]
SEPscL310_seuratobj <- FindClusters(SEPscL310_seuratobj,resolution = 0.5)

##plot the cluster
DimPlot(SEPscL310_seuratobj,reduction="pca",label = TRUE)+ggtitle("PC1 vs PC2 with Clusters")

##if you want to look at different PCs
DimPlot(SEPscL310_seuratobj,reduction="pca", dims=c(3,4), label=TRUE)+ggtitle("PC3 vs PC4 with Clusters")

###look to see if the clusters are being influenced by any of the QC metrics
##number of reads
VlnPlot(SEPscL310_seuratobj,features="nCount_RNA")

##number of genes
VlnPlot(SEPscL310_seuratobj,features="nFeature_RNA")

##cell cycle
SEPscL310_seuratobj@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")   


###Finding Markers for each Cluster
SEPscL310_seuratobj.markers <- FindAllMarkers(SEPscL310_seuratobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SEPscL310_seuratobj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write_csv(SEPscL310_seuratobj.markers, file = "/home/bzawisza/seurat/SEPscL310/DEXanalysisNewIntegration.csv")

###We can extract from this list the most upregulated gene from each cluster
SEPscL310_seuratobj.markers %>%
  group_by(cluster) %>%
  slice(1) %>%
  pull(gene) -> best.wilcox.gene.per.cluster

best.wilcox.gene.per.cluster

VlnPlot(SEPscL310_seuratobj,features=best.wilcox.gene.per.cluster)

###identify the unique genes in each clusters
##look at the power value, 1 =perfectly separating,0 = random
SEPscL310_seuratobj.markers.ROC <-FindAllMarkers(SEPscL310_seuratobj, test.use = "roc", only.pos = TRUE)

##plot out genes that you think it might be unique to the clusters
#VlnPlot(SEPscL310_seuratobj,features="##gene of interest##")


###Automated Cell Type Annotation
as.data.frame(SEPscL310_seuratobj@assays$RNA[,]) -> scina.data
load(system.file('extdata','example_signatures.RData', package = "SCINA"))

signatures

SCINA(
  scina.data,
  signatures, 
  max_iter = 100, 
  convergence_n = 10, 
  convergence_rate = 0.999, 
  sensitivity_cutoff = 0.9, 
  rm_overlap=TRUE, 
  allow_unknown=TRUE
) -> scina.results

SEPscL310_seuratobj$scina_labels <- scina.results$cell_labels

##plot cell type annotation
DimPlot(SEPscL310_seuratobj,reduction = "pca", pt.size = 1, label = TRUE, group.by = "scina_labels", label.size = 5)
