# single-cell custom functions


#' plotBarcodesPerCell
#'
#' plot the number of barcodes detected per cell.
#' For use with barcode labelled cells. also useful for single cell CRISPR screens
#'
#' @param obj a Seurat object
#' @param samplename sample name
#'
#' @return Returns a plot of the number of barcodes detected per cell
#'
#' @import Seurat
#' @export
#'

plotBarcodesPerCell <- function(obj, samplename = "Seurat.obj"){
  meta.data <- obj@meta.data
  counts <- c()
  # [TO-DO] use apply here on rows of metadata with a custom function to speedup
  for (i in 1:nrow(meta.data)){
    cell <- meta.data[i,]
    if (cell$barcode != "not.detected"){
      num <- length(unlist(strsplit(as.character(cell$barcode), ";")))
      counts <- c(counts,num)
    }
  }
  print(table(counts))
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = counts), binwidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste("Num Barcodes per Cell:", samplename))
  print(p)
}

#' convertHumanGeneList
#'
#' convert human gene list to mouse gene symbols
#'
#' @param genes a list of genes to convert
#'
#' @return Returns a list of converted gene symbols.
#'
#' @import biomaRt
#' @export
#'

convertHumanGeneList <- function(genes){
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = genes , mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows=T)

  humanx <- unique(genesV2[, 2])

  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#' convertGeneListFormat
#'
#' convert gene list format between ENSEMBL ENTREZ and HGNC sybmol
#'
#' @param genes a list of gene symbols
#' @param species species of the gene symbols. One of "Hs" (human) or "Mm" (mouse)
#' @param from gene list format to convert from
#' @param to gene list format to convert to
#'
#' @return Returns a data-frame containing the gene symbols converted to new format.
#'
#' @export
#'

# convert between different gene name formats using bitr
convertGeneListFormat <- function(genes, species = "Hs", from = "SYMBOL", to = "ENTREZID"){

  # import Org.Db
  if (species == "Hs"){
    org.Hs.eg.db::org.Hs.eg.db
    db <- org.Hs.eg.db::org.Hs.eg.db
  }
  if (species == "Mm"){
    org.Mm.eg.db::org.Mm.eg.db
    db <- org.Mm.eg.db::org.Mm.eg.db
  }
  df <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = db, drop = F)
  return(df)
}


#' findDoubletsByBarcode
#'
#' Find cells that have a unique combination of barcodes. Probably doublets
#'
#' @param obj a Seurat object
#' @param threshold number of times barcode pair must be observed for cell to be considered not a doublet
#'
#' @return Returns a Seurat object containing doubletBarcode metadata field
#'
#' @export
#'

findDoubletsByBarcode <- function(obj, threshold = 2){
  doublets.by.barcode <- c()

  # assert barcodes are present in object
  if(is.null(obj$barcode)){
    message("no barcode field identified in object metadata. Stopping")
    break()
  }

  # identify cells that have more than 1 barcode
  multbarcodes<-obj$barcode[grep(";",obj$barcode)]
  sortbc <- c()
  for (i in 1:length(multbarcodes)){
    barcodes <- str_split(multbarcodes[i], pattern = ";")
    barcodes <- lapply(barcodes,sort)
    barcodes <- paste(unlist(barcodes),collapse=";")
    names(barcodes) <- names(multbarcodes[i])
    sortbc <- c(sortbc,barcodes)

  }

  #combination of barcodes that only appear once
  sortbc <- as.data.frame(sortbc)
  doubletbc <- sortbc[which(!(duplicated(sortbc) | duplicated(sortbc, fromLast=TRUE))),,drop=F]

  # find barcodes that

  #insert metadata into seruat obj
  obj$doubletBarcode <- ifelse(rownames(obj@meta.data) %in% rownames(doubletbc), yes = "doublet", no = "singlet")
  unknown <- which(obj$barcode == "not.detected")
  obj$doubletBarcode[unknown] <- "unknown"

  return(obj)
}

#' runClusterProfiler
#'
#' run clusterProfiler on a set of dge results from Seurat::FindMarkers
#'
#' @param dge a Seurat object
#' @param lfc.threshold fold change threshold
#' @param padj.threshold adjusted p value threshold
#' @param OrgDb org db object for relevant species. org.Mm.eg.db or org.Hs.eg.db
#' @param category GO category to test. e.g "BP", "MF", or "CC"
#' @param sample name of sample gene set
#' @param outdir path to desired output directory
#'
#' @return Returns plots and a table of clusterProfiler results
#'
#' @export
#'

runClusterProfiler <- function(dge, lfc.threshold = 0.2, padj.threshold = 0.1, OrgDb = org.Mm.eg.db,
                               category = "BP", sample = "gene-set", outdir = NULL){
  message(paste("Running ClusterProfiler on", sample))
  if (category == "BP"){
    message("Examining Biological Process categories")
  }
  if (category == "MF"){
    message("Examining Molecular Function categories")
  }
  if (category == "CC"){
    message("Examining Cellular Component categories")
  }
  geneset.up <- rownames(dge[which(dge$avg_logFC > lfc.threshold),])
  geneset.dn <- rownames(dge[which(dge$avg_logFC < -lfc.threshold),])

  # convert gene symbols to Entrez id
  gene.df.up <- clusterProfiler::bitr(geneset.up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  gene.df.dn <- clusterProfiler::bitr(geneset.dn, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)

  # get enriched GO terms
  ego.up <- clusterProfiler::enrichGO(gene = gene.df.up$ENTREZID,
                     OrgDb = OrgDb, ont = category,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  head(ego.up, 20)

  ego.dn <- clusterProfiler::enrichGO(gene = gene.df.dn$ENTREZID,
                     OrgDb = OrgDb, ont = category,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  head(ego.dn, 20)

  print(clusterProfiler::barplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
  print(clusterProfiler::barplot(ego.dn, showCategory=15) + ggtitle(paste("Downregulated", category, "GO terms in", sample)))

  print(clusterProfiler::emapplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
  print(clusterProfiler::emapplot(ego.dn, showCategory=15)+ ggtitle(paste("Downregulated", category, "GO terms in", sample)))

  #ego3 <- gseGO(geneList = gene.df.up,
  #              OrgDb = OrgDb,
  #              ont  = category,
  #              nPerm  = 1000,
  #              minGSSize = 50,
  #              maxGSSize = 500,
  #              pvalueCutoff = 0.05,
  #              verbose = T,
  #              keyType = "ENTREZID")

  #print(ego3)

  if (!is.null(outdir)){
    plot.dir <- file.path(outdir)
    pdf(file.path(plot.dir,paste(sample,"_clusterProfiler_output.pdf", sep = '')), useDingbats = F)
    print(clusterProfiler::barplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
    print(clusterProfiler::barplot(ego.dn, showCategory=15) + ggtitle(paste("Downregulated", category, "GO terms in", sample)))
    print(clusterProfiler::emapplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
    print(clusterProfiler::emapplot(ego.dn, showCategory=15)+ ggtitle(paste("Downregulated", category, "GO terms in", sample)))
    dev.off()
  }
}

#' findDoubletsBySimulation
#'
#' Find cells in a Seurat / SCE object that may be doublets based on a few different simulation methods from scds and scran packages
#'
#' @param obj a Seurat or SingleCellExperiment object
#' @param ntop number of top variable genes to use for simulations
#' @param mads median absolute deviations to use as threshold for doublet calls
#' @param seed random seed to use
#'
#' @return Returns a Seurat object containing doubletBarcode metadata field
#'
#' @export
#'

findDoubletsBySimulation <- function(obj, ntop = 1000, mads = 3, seed = 10101){
  message("Running doublet detection on input single cell object using scds and scran methods")

  # setup random seed
  set.seed(seed)

  # identify single cell object class
  sample <- obj
  if (class(sample) == "Seurat"){
    message("Seurat input detected, converting to SingleCellExperiment object")
    require(Seurat)
    sample.sce <- Seurat::as.SingleCellExperiment(sample)
  }

  if (class(sample) == "SingleCellExperiment"){
    message("SCE input detected, continuing")
    require(SingleCellExperiment)
    sample.sce <- obj
  }

  # Run scds doublet detection
  message("")
  message("#- Annotate doublets using co-expression based doublet scoring:")
  sample.sce <- cxds(sce = sample.sce, retRes = T, verb = T, estNdbl = T, ntop = ntop)

  message("")
  message("#- Annotate doublets using simulation approach:")
  sample.sce <- bcds(sample.sce, retRes = T, verb = T, ntop = ntop, varImp = T, estNdbl = T)

  message("")
  message("#- Annotate doublets using hybrid co-expression and simulation approach:")
  sample.sce <- cxds_bcds_hybrid(sample.sce, verb = T, estNdbl = T)

  # Run scran simulation approach
  message("")
  message("#- Annotate doublets using scran doublet simulation approach:")
  dbl.dens <- doubletCells(sample.sce, d = ncol(reducedDim(sample.sce)))
  summary(dbl.dens)
  sample.sce$DoubletScore <- log10(dbl.dens+1)
  sample.sce$dbl.dens <- dbl.dens

  # plot distribution of scores
  par(mfcol=c(1,4))
  boxplot(sample.sce$cxds_score ~ sample.sce$orig.ident, main="cxds (scds co-expression)")
  boxplot(sample.sce$bcds_score ~ sample.sce$orig.ident, main="bcds (scds simulation)")
  boxplot(sample.sce$hybrid_score ~ sample.sce$orig.ident, main="hybrid (scds hybrid)")
  boxplot(sample.sce$DoubletScore ~ sample.sce$orig.ident, main="scran simulation")

  # plot
  if ("UMAP" %in% SingleCellExperiment::reducedDimNames(sample.sce)){
    plotReducedDim(sample.sce, dimred = "UMAP", colour_by = "cxds_call")
    plotReducedDim(sample.sce, dimred = "UMAP", colour_by = "hybrid_call")
    plotReducedDim(sample.sce, dimred = "UMAP", colour_by = "bcds_call")
    plotReducedDim(sample.sce, dimred = "UMAP", colour_by = "DoubletScore")
  }

  # summary tables of outliers based on mads
  table(isOutlier(sample.sce$DoubletScore, type = "higher", nmads = mads))
  table(isOutlier(sample.sce$cxds_score, type = "higher", nmads = mads))
  table(isOutlier(sample.sce$bcds_score, type = "higher", nmads = mads))
  table(isOutlier(sample.sce$hybrid_score, type = "higher", nmads = mads))

  # calculate outliers based on mads
  message("")
  message("Calling outliers based on mads")
  sample.sce$doubletscore_outlier <- isOutlier(sample.sce$DoubletScore, type = "higher", nmads = mads)
  sample.sce$cxds_outlier <- isOutlier(sample.sce$cxds_score, type = "higher", nmads = mads)
  sample.sce$bcds_outlier <- isOutlier(sample.sce$bcds_score, type = "higher", nmads = mads)
  sample.sce$hybrid_outlier <- isOutlier(sample.sce$hybrid_score, type = "higher", nmads = mads)

  plotColData(sample.sce, x = "doubletscore_outlier", y = "DoubletScore", colour_by = "DoubletScore")
  plotColData(sample.sce, x = "orig.ident", y = "DoubletScore", colour_by = "DoubletScore")

  if (class(sample) == "Seurat"){
    message("Merging doublet calls into Seurat object")
    sample$cxds_score <- sample.sce$cxds_score
    sample$cxds_call <- sample.sce$cxds_call
    sample$bcds_score <- sample.sce$bcds_score
    sample$bcds_call <- sample.sce$bcds_call
    sample$hybrid_score <- sample.sce$hybrid_score
    sample$hybrid_call <- sample.sce$hybrid_call
    sample$DoubletScore <- sample.sce$DoubletScore

    sample$doubletscore_outlier <- sample.sce$doubletscore_outlier
    sample$cxds_outlier <- sample.sce$cxds_outlier
    sample$bcds_outlier <- sample.sce$bcds_outlier
    sample$hybrid_outlier <- sample.sce$hybrid_outlier

    message("doublet annotation complete!")
    return(sample)
  } else {
    return(sample.sce)
  }
  table(rownames(sample@meta.data) == rownames(colData(sample.sce)))
}


#' run_edgeR_LRT
#'
#' Run edge R for 2 groups Ctrl vs Test using the Likelihood Ratio Test (LRT)
#' Adapted from Sonesson et al. 2019?? Nat comms
#' https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRLRT.R
#' This is a simple A vs B test. See other functions for more complex GLMs
#'
#' @param counts raw counts from single cell object, rows as features and columns as cells
#' @param group which factor in object metadata to use for testing. Should only have 2 levels
#' @param plots Logical. print plots
#'
#' @return Returns a Seurat object containing doubletBarcode metadata field
#'
#' @import edgeR
#'
#' @export

run_edgeR_LRT <- function(counts, group, plots = F){

  message("Running edgeR LRT")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- edgeR::DGEList(counts, group = group)
    message("Calculate normalisation factors")
    dge <- edgeR::calcNormFactors(dge)
    message("Generate model matrix")
    design <- edgeR::model.matrix(~ group)
    print(design)
    message("Estimate dispersions")
    dge <- edgeR::estimateDisp(dge, design = design)
    message("Fit GLM model")
    fit <- edgeR::glmFit(dge, design = design)
    message("Likelihood ratio test")
    lrt <- edgeR::glmLRT(fit)
    message("Export top tags")
    tt <- edgeR::topTags(lrt, n = Inf)
  })

  # plots
  if(isTRUE(plots)){
    message("Plotting, this could take a while...")
    edgeR::plotBCV(dge)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(group)), pch = 19)
    edgeR::plotSmear(lrt)
    message("Plotting complete")
  }

  # export results
  message("generate results")
  results <- list(session_info = session_info,
                  timing = timing,
                  tt = tt$table,
                  df = data.frame(PValue = tt$table$PValue,
                                  FDR = tt$table$FDR,
                                  logFC = tt$table$logFC,
                                  logCPM = tt$table$logCPM,
                                  rownames = rownames(tt$table)))
  message("complete")
  return(results)
}
