# runSeuratPipeline
# Dane Vassiliadis
# 8-6-20

# Pipeline to automate basic 10X sample processing and visualisation using
# Seurat


#-----------------------
# Required libs
#-----------------------
require(Seurat)
require(sctransform)
require(DropletUtils)
require(PCAtools)
require(biomaRt)

#-----------------------
# Main function
#-----------------------
runSeuratPipeline <- function(counts.dir, outdir = "./", norm = "sctransform",
  mitoCutoff = 20,featuresLower = 100L,
  featuresUpper = 10000L,countsLower = 1000L,
  countsUpper = 100000L, dims = 1:30, verbose = T,
  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
  samplename = "seurat_pipeline_default",
  regressVars = c("percent.mito", "S.Score", "G2M.Score"),
  cellcycleGenes = NULL, species = "mouse", project = "scRNAseq",
  meta = NULL, seed.use  = 10101){

  #-----------------------
  # Run pipeline
  #-----------------------
  if(is.null(counts.dir)){
    message("no input directory given. quitting")
    quit(status = 1)
  } else {
    counts.dir = file.path(counts.dir)
  }

  message("-------------------------")
  message("Seurat scRNAseq pipeline")
  message("-------------------------")
  message(paste("Processing:", samplename))
  message(paste("Input directory:", counts.dir))
  message(paste("Output directory:", outdir))
  message(paste("Species:", species))
  message(paste("Project:", project))
  message(paste("Mitochondrial % cutoff:", mitoCutoff))
  message(paste("# features cutoff:", featuresLower, featuresUpper))
  message(paste("# read counts cutoff:", countsLower, countsUpper))
  message(paste("Variables to regress during normalisation:", toString(regressVars)))
  message(paste("PCA dimensions to use:", toString(dims)))
  message(paste("Clustering resolutions to use:", toString(resolution)))
  #message(paste("Metadata:", meta))
  message(paste("Random seed:", seed.use))

  # setup output dirs
  output.dir <- file.path(outdir)
  dir.create(output.dir, showWarnings = F)

  plot.dir <- file.path(outdir, "plots")
  dir.create(plot.dir, showWarnings = F)

  qc.dir <- file.path(outdir, "qc")
  dir.create(qc.dir, showWarnings = F)

  Seurat.obj.dir <- file.path(outdir, "seurat_obj")
  dir.create(Seurat.obj.dir, showWarnings = F)

  message("")
  message("Importing count data")
  counts <- Read10X(counts.dir, strip.suffix = T)
  sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = project, meta.data = meta)

  #-----------------------
  # Sample QC
  #-----------------------
  message("Performing sample QC")
  mito.pattern <- "^mt-"
  if(species == "human" | species == "Human" | species == "Hs" | species == "hs"){
    mito.pattern <- "^MT-"
  }

  # calculate mito percentage
  message("Calculating mito %")
  mito.genes <- grep(pattern = mito.pattern, x = rownames(x = sample@assays$RNA), value = TRUE)
  percent.mito <- Matrix::colSums(sample@assays$RNA[mito.genes, ])/Matrix::colSums(sample@assays$RNA)*100
  sample <- AddMetaData(object = sample, metadata = percent.mito, col.name = "percent.mito")

  # plot QC metrics
  message("Plotting QC metrics - Pre-filter")
  pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_pre-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T))
  dev.off()

  # subset based on mito count and feature cutoffs
  message("Filtering low quality cells")
  message("Starting dimensions")
  print(dim(sample))
  sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoCutoff)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA > countsLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])

  message("")
  message("Filtered dimensions")
  print(dim(sample))

  # plot QC metrics
  message("Plotting QC metrics - Post-filter")
  pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T))
  dev.off()

  # Feature Scatter plots
  pdf(file.path(qc.dir,paste(samplename, "_count_mito_QCscatter-postfilter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  dev.off()
  pdf(file.path(qc.dir,paste(samplename, "_count_feature_QCscatter-post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()

  #-----------------------
  # Normalisation #1
  #-----------------------
  message(paste0("Normalising counts using: ", norm))
  # SCTransform workflow
  message("Regressing following variables: percent.mito")

  if(norm == "sctransform") {
    sample <- SCTransform(sample, vars.to.regress = "percent.mito", verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use = seed.use)
  }

  # log normalisation workflow
  if(norm == "lognorm") {
    sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                              vars.to.regress = regress.vars,
                              verbose = F, conserve.memory = F)
  }

  #----------------------------------------------
  # Seurat workflow without cell cycle regression
  #----------------------------------------------
  message("Seurat workflow without cell cycle regression")

  # run PCA
  message("Running PCA")
  sample <- RunPCA(sample, npcs = 100, verbose = T, seed.use = seed.use)

  # Run tSNE & UMAP
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))

  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

  message("Finding nearest neighbors")
  sample <- FindNeighbors(object = sample, verbose = T, dims = dims)

  message("Running cluster sequence")
  # run a sequence of cluster resolutions
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

  message("Saving normalised dataset pre cell cycle regression")
  # save normalised and clustered dataset
  saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_normalised_no_cc_regression.rds", sep = "")))

  #----------------------------------------------
  # Seurat workflow with cell cycle regression
  #----------------------------------------------

  #-----------------------
  # Cell cycle scoring
  #-----------------------

  message("Performing Cell cycle scoring")
  # perform cell cycle analysis
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

  if (species == "mouse"){
    message("Species is mouse. Converting human to mouse gene symbols")
    convertHumanGeneList <- function(x){
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
        values = x , mart = human, attributesL = c("mgi_symbol"),
        martL = mouse, uniqueRows=T)
      humanx <- unique(genesV2[, 2])
      return(humanx)
    }
    s.genes <- convertHumanGeneList(s.genes)
    g2m.genes <- convertHumanGeneList(g2m.genes)
    message("mouse cell cycle genes to use: S Phase")
    print(s.genes)
    message("mouse cell cycle genes to use: G2M Phase")
    print(g2m.genes)
  }

  sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, search = F)

  # add in cell cycle difference scoring
  sample$CC.Difference <- sample$S.Score - sample$G2M.Score

  # run PCA on cell cycle components
  sample <- RunPCA(sample, features = c(s.genes, g2m.genes))

  pdf(file.path(plot.dir,paste(samplename, "_PCA_cell_cycle_anno.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5))
  dev.off()

  #-----------------------
  # Normalisation #2
  #-----------------------
  message("Normalisation #2")
  message(paste("Normalising counts using:", norm))
  # SCTransform workflow
  message("Regressing following variables:")

  regressVars <- str_split(regressVars, ',')[[1]]
  print(regressVars)

  if(norm == "sctransform") {
    sample <- SCTransform(sample, vars.to.regress = regressVars, verbose = T, conserve.memory = F, return.only.var.genes = F)
  }

  # log normalisation workflow
  if(norm == "lognorm") {
    sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                              vars.to.regress = regress.vars,
                              verbose = F, conserve.memory = F)
  }

  #--- Dimensional reduction & clustering
  message("Running dimensional reduction algos")

  # always calculate 100 PCs
  sample <- RunPCA(object = sample, verbose = T, npcs = 100, seed.use = seed.use)
  message("Plotting...")

  # Elbow plot
  percent.var <- Stdev(sample)
  elbow.dim <- PCAtools::findElbowPoint(percent.var)
  message(paste("Predicted elbow PC =", elbow.dim))
  pdf(file.path(plot.dir,paste(samplename,"_elbow_plot.pdf", sep = '')), width = 8, height = 6, useDingbats = F)
  print(ElbowPlot(sample, ndims = 50) +
    geom_vline(aes(xintercept = elbow.dim), color = "red") +
    labs(title = paste("Elbow plot -",samplename), subtitle = paste("Predicted elbow PC =", elbow.dim)))
  dev.off()

  # Run tSNE & UMAP on 1 to elbow dimensions
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="umap_elbow")
  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="tsne_elbow")

  # find neighbors based on 1 to elbow dimensions
  sample <- FindNeighbors(object = sample, verbose = T, dims = 1:elbow.dim)

  # Run tSNE & UMAP on 1 to dims
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))
  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

  # run a sequence of cluster resolutions for elbow dim
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

  # save normalised and clustered dataset
  message("Saving normalised dataset post cell cycle regression")
  saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_cc_regressed_normalised.rds", sep = "")))

  message("performing final plots on normalised dataset")

  # PCA plot reduction on clusters
  p <- DimPlot(object = sample, reduction = "pca", label = T, label.size = 5, dims = c(1,2)) +
    NoLegend() +
    ggtitle(paste('PCA -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_PCA.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # UMAP reduction on clusters - elbow
  p <- DimPlot(object = sample, reduction = "umap_elbow", label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('UMAP -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_UMAP_elbow.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # t-SNE reduction on clusters - elbow
  p <- DimPlot(object = sample, reduction = "tsne_elbow", label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('t-SNE -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_t-SNE_elbow.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # # UMAP reduction on clusters - dims
  # p <- DimPlot(object = sample, reduction = paste("umap",tail(dims, 1)), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
  #   NoLegend() +
  #   ggtitle(paste('UMAP -', samplename))
  # p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  # ggsave(file.path(plot.dir,paste(samplename, paste("_UMAP_",tail(dims, 1)), sep=''), ".pdf", sep = '')), width = 9, height = 7)

  # # t-SNE reduction on clusters - dims
  # p <- DimPlot(object = sample, reduction = paste("tsne",dims), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
  #   NoLegend() +
  #   ggtitle(paste('t-SNE -', samplename))
  # p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  # ggsave(file.path(plot.dir,paste(samplename, paste("_t-SNE_",dims), sep=''), ".pdf", sep = '')), width = 9, height = 7)

  message(paste("Finished", samplename, "analysis. Enjoy your data!"))

  sessionInfo()
}


#-----------------------
# Main function with Seurat obj input
#-----------------------
runSeuratPipelineObj <- function(obj, outdir = "./", norm = "sctransform",
  mitoCutoff = 20,featuresLower = 100L,
  featuresUpper = 10000L,countsLower = 1000L,
  countsUpper = 100000L, dims = 1:30, verbose = T,
  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
  samplename = "seurat_pipeline_default",
  regressVars = c("percent.mito", "S.Score", "G2M.Score"),
  cellcycleGenes = NULL, species = "mouse", project = "scRNAseq",
  meta = NULL, seed.use  = 10101){

  #-----------------------
  # Run pipeline
  #-----------------------
  if(is.null(obj)){
    message("no input object given. quitting")
    quit(status = 1)
  }

  message("-------------------------")
  message("Seurat scRNAseq pipeline")
  message("-------------------------")
  message(paste("Processing:", samplename))
  message("Input sample:")
  obj
  message(paste("Output directory:", outdir))
  message(paste("Species:", species))
  message(paste("Project:", project))
  message(paste("Mitochondrial % cutoff:", mitoCutoff))
  message(paste("# features cutoff:", featuresLower, featuresUpper))
  message(paste("# read counts cutoff:", countsLower, countsUpper))
  message(paste("Variables to regress during normalisation:", toString(regressVars)))
  message(paste("PCA dimensions to use:", toString(dims)))
  message(paste("Clustering resolutions to use:", toString(resolution)))
  #message(paste("Metadata:", meta))
  message(paste("Random seed:", seed.use))

  # setup output dirs
  output.dir <- file.path(outdir)
  dir.create(output.dir, showWarnings = F)

  plot.dir <- file.path(outdir, "plots")
  dir.create(plot.dir, showWarnings = F)

  qc.dir <- file.path(outdir, "qc")
  dir.create(qc.dir, showWarnings = F)

  Seurat.obj.dir <- file.path(outdir, "seurat_obj")
  dir.create(Seurat.obj.dir, showWarnings = F)

  message("")
  message("Importing count data")
  #counts <- Read10X(counts.dir, strip.suffix = T)
  #sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = project, meta.data = meta)
  sample <- obj

  #-----------------------
  # Sample QC
  #-----------------------
  message("Performing sample QC")
  mito.pattern <- "^mt-"
  if(species == "human" | species == "Human" | species == "Hs" | species == "hs"){
    mito.pattern <- "^MT-"
  }

  # calculate mito percentage
  message("Calculating mito %")
  mito.genes <- grep(pattern = mito.pattern, x = rownames(x = sample@assays$RNA), value = TRUE)
  percent.mito <- Matrix::colSums(sample@assays$RNA[mito.genes, ])/Matrix::colSums(sample@assays$RNA)*100
  sample <- AddMetaData(object = sample, metadata = percent.mito, col.name = "percent.mito")

  # plot QC metrics
  message("Plotting QC metrics - Pre-filter")
  pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_pre-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T))
  dev.off()

  # subset based on mito count and feature cutoffs
  message("Filtering low quality cells")
  message("Starting dimensions")
  print(dim(sample))
  sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoCutoff)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA > countsLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])

  message("")
  message("Filtered dimensions")
  print(dim(sample))

  # plot QC metrics
  message("Plotting QC metrics - Post-filter")
  pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T))
  dev.off()

  # Feature Scatter plots
  pdf(file.path(qc.dir,paste(samplename, "_count_mito_QCscatter-postfilter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  dev.off()
  pdf(file.path(qc.dir,paste(samplename, "_count_feature_QCscatter-post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()

  #-----------------------
  # Normalisation #1
  #-----------------------
  message(paste0("Normalising counts using: ", norm))
  # SCTransform workflow
  message("Regressing following variables: percent.mito")

  if(norm == "sctransform") {
    sample <- SCTransform(sample, vars.to.regress = "percent.mito", verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use = seed.use)
  }

  # log normalisation workflow
  if(norm == "lognorm") {
    sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                            vars.to.regress = regress.vars,
                            verbose = F, conserve.memory = F)
  }

  #----------------------------------------------
  # Seurat workflow without cell cycle regression
  #----------------------------------------------
  message("Seurat workflow without cell cycle regression")

  # run PCA
  message("Running PCA")
  sample <- RunPCA(sample, npcs = 100, verbose = T, seed.use = seed.use)

  # Run tSNE & UMAP
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))

  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

  message("Finding nearest neighbors")
  sample <- FindNeighbors(object = sample, verbose = T, dims = dims)

  message("Running cluster sequence")
  # run a sequence of cluster resolutions
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

  message("Saving normalised dataset pre cell cycle regression")
  # save normalised and clustered dataset
  saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_normalised_no_cc_regression.rds", sep = "")))

  #----------------------------------------------
  # Seurat workflow with cell cycle regression
  #----------------------------------------------

  #-----------------------
  # Cell cycle scoring
  #-----------------------

  message("Performing Cell cycle scoring")
  # perform cell cycle analysis
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

  if (species == "mouse"){
    message("Species is mouse. Converting human to mouse gene symbols")
    convertHumanGeneList <- function(x){
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                       values = x , mart = human, attributesL = c("mgi_symbol"),
                       martL = mouse, uniqueRows=T)
      humanx <- unique(genesV2[, 2])
      return(humanx)
    }
    s.genes <- convertHumanGeneList(s.genes)
    g2m.genes <- convertHumanGeneList(g2m.genes)
    message("mouse cell cycle genes to use: S Phase")
    print(s.genes)
    message("mouse cell cycle genes to use: G2M Phase")
    print(g2m.genes)
  }

  sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, search = F)

  # add in cell cycle difference scoring
  sample$CC.Difference <- sample$S.Score - sample$G2M.Score

  # run PCA on cell cycle components
  sample <- RunPCA(sample, features = c(s.genes, g2m.genes))

  pdf(file.path(plot.dir,paste(samplename, "_PCA_cell_cycle_anno.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
  print(DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5))
  dev.off()

  #-----------------------
  # Normalisation #2
  #-----------------------
  message("Normalisation #2")
  message(paste("Normalising counts using:", norm))
  # SCTransform workflow
  message("Regressing following variables:")

  regressVars <- str_split(regressVars, ',')[[1]]
  print(regressVars)

  if(norm == "sctransform") {
    sample <- SCTransform(sample, vars.to.regress = regressVars, verbose = T, conserve.memory = F, return.only.var.genes = F)
  }

  # log normalisation workflow
  if(norm == "lognorm") {
    sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                            vars.to.regress = regress.vars,
                            verbose = F, conserve.memory = F)
  }

  #--- Dimensional reduction & clustering
  message("Running dimensional reduction algos")

  # always calculate 100 PCs
  sample <- RunPCA(object = sample, verbose = T, npcs = 100, seed.use = seed.use)
  message("Plotting...")

  # Elbow plot
  percent.var <- Stdev(sample)
  elbow.dim <- PCAtools::findElbowPoint(percent.var)
  message(paste("Predicted elbow PC =", elbow.dim))
  pdf(file.path(plot.dir,paste(samplename,"_elbow_plot.pdf", sep = '')), width = 8, height = 6, useDingbats = F)
  print(ElbowPlot(sample, ndims = 50) +
          geom_vline(aes(xintercept = elbow.dim), color = "red") +
          labs(title = paste("Elbow plot -",samplename), subtitle = paste("Predicted elbow PC =", elbow.dim)))
  dev.off()

  # Run tSNE & UMAP on 1 to elbow dimensions
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="umap_elbow")
  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="tsne_elbow")

  # find neighbors based on 1 to elbow dimensions
  sample <- FindNeighbors(object = sample, verbose = T, dims = 1:elbow.dim)

  # Run tSNE & UMAP on 1 to dims
  message("Running UMAP")
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))
  message("Running TSNE")
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

  # run a sequence of cluster resolutions for elbow dim
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

  # save normalised and clustered dataset
  message("Saving normalised dataset post cell cycle regression")
  saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_cc_regressed_normalised.rds", sep = "")))

  message("performing final plots on normalised dataset")

  # PCA plot reduction on clusters
  p <- DimPlot(object = sample, reduction = "pca", label = T, label.size = 5, dims = c(1,2)) +
    NoLegend() +
    ggtitle(paste('PCA -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_PCA.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # UMAP reduction on clusters - elbow
  p <- DimPlot(object = sample, reduction = "umap_elbow", label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('UMAP -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_UMAP_elbow.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # t-SNE reduction on clusters - elbow
  p <- DimPlot(object = sample, reduction = "tsne_elbow", label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('t-SNE -', samplename))

  pdf(file.path(plot.dir,paste(samplename, "_t-SNE_elbow.pdf", sep = '')), width = 9, height = 7)
  print(AugmentPlot(p, dpi = 600, width = 8, height = 6))
  dev.off()

  # # UMAP reduction on clusters - dims
  # p <- DimPlot(object = sample, reduction = paste("umap",tail(dims, 1)), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
  #   NoLegend() +
  #   ggtitle(paste('UMAP -', samplename))
  # p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  # ggsave(file.path(plot.dir,paste(samplename, paste("_UMAP_",tail(dims, 1)), sep=''), ".pdf", sep = '')), width = 9, height = 7)

  # # t-SNE reduction on clusters - dims
  # p <- DimPlot(object = sample, reduction = paste("tsne",dims), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
  #   NoLegend() +
  #   ggtitle(paste('t-SNE -', samplename))
  # p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  # ggsave(file.path(plot.dir,paste(samplename, paste("_t-SNE_",dims), sep=''), ".pdf", sep = '')), width = 9, height = 7)

  message(paste("Finished", samplename, "analysis. Enjoy your data!"))

  sessionInfo()
}













