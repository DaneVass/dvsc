# single-cell custom functions



# custom functions to output rasterised versions of plots for ease of use in illustrator. this was so incredibly frustrating
# Using the ggpubr package. just print out the legend. Useful for AugmentPlot plots because the legends are lost
getlegend <- function(plot){
  require(ggpubr)
  # Extract the legend. Returns a gtable
  legend <- get_legend(plot)
  # Convert to a ggplot and print
  as_ggplot(legend)
}

# these functions are to be used under pdf statements
plotAugmentDim <- function(obj, reduction = "umap", label = F, label.size = 6, group.by = NULL, title = NULL, dpi = 600, pt.size = 1, cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 2, split.by = NULL){
  p <- DimPlot(obj, reduction = reduction, label = label, label.size = label.size, group.by = group.by, pt.size = pt.size, cells.highlight = cells.highlight, cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, split.by = split.by) +
    ggtitle(title)
  print(AugmentPlot(p, dpi = dpi))
  print(getlegend(p))
}

plotAugmentFeature <- function(obj, reduction = "umap", features = NULL, label = F, label.size = 6, group.by = NULL, title = NULL, dpi = 600, pt.size = 1, order = T, ncol = NULL){
  p <- FeaturePlot(obj, reduction = reduction, features = features, label = label, label.size = label.size, pt.size = pt.size, order = order, ncol = ncol) +
    ggtitle(title)
  print(AugmentPlot(p, dpi = dpi))
  print(getlegend(p))
}

plotAugmentScatter <-


plotElbow <- function(obj, outdir = NULL, samplename = "sample"){
  require(PCAtools)
  require(Seurat)
  percent.var <- Stdev(obj)
  elbow.dim <- PCAtools::findElbowPoint(percent.var)
  p <- Seurat::ElbowPlot(obj, ndims = 50) +
    geom_vline(aes(xintercept = elbow.dim), color = "red") +
    labs(title = paste("Elbow plot -", samplename), subtitle = paste("Predicted elbow PC =", elbow.dim))

  if(is.null(outdir)){
    print(p)
  } else {
    ggsave(file.path(outdir,paste(prefix,".pdf", sep = '')), width = 8, height = 6)
  }
}

match.barcodes.to.cells <- function(all.cells, cells.w.barcode.df){
    # matches cells in a single cell experiment to detected DNA barcodes.
    # all cells is a list of all cells in the experiment
    # cells.w.barcode.df is a dataframe with cell id and barcode id as columns
    # dataframe returned will have all cells matched to a barcode
    # if there is no barcode matchable to a cell "not.detected" is returned
    # for cells that have multiple detected barcodes each barcode is returned separated by ';'
    cell.barcode.annotation <- data.frame()
    for (id in all.cells){
        keep <- which(cells.w.barcode.df$Cell.10X.Barcode == id)
        df <- cells.w.barcode.df[keep,]
        unique.barcodes <- length(unique(df$referenceID))
        if (unique.barcodes == 0){
            df.2 <- data.frame(cell.id = id, barcode = "not.detected")
            cell.barcode.annotation <- rbind(df.2, cell.barcode.annotation)
        }
        if (unique.barcodes == 1){
            df.2 <- data.frame(cell.id = id, barcode = df[1,2])
            cell.barcode.annotation <- rbind(df.2, cell.barcode.annotation)
        }
        if (unique.barcodes > 1){
            barcodes <- paste(unique(df$referenceID), collapse = ";")
            df.2 <- data.frame(cell.id = id, barcode = barcodes)
            cell.barcode.annotation <- rbind(df.2, cell.barcode.annotation)
        }
    }
    return(cell.barcode.annotation)
}

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
  p <- ggplot() + geom_histogram(aes(x = counts), binwidth = 1) +
    theme_bw() +
    ggtitle(paste("Num Barcodes per Cell:", samplename))
  print(p)
}

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = x , mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows=T)

  humanx <- unique(genesV2[, 2])

  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

# convert between different gene name formats using bitr
convertGeneListFormat <- function(x, species = "Hs", from = "SYMBOL", to = "ENTREZID", ...){
  require(clusterProfiler)
  require(bitr)

  # import Org.Db
  if (species == "Hs"){
    library(org.Hs.eg.db)
    db <- org.Hs.eg.db
  }
  if (species == "Mm"){
    library(org.Mm.eg.db)
    db <- org.Mm.eg.db
  }
  df <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = F)
  return(df)
}

findDoubletsByBarcode.2 <- function(obj){
  message("Finding doublets by Barcode")
  doublets.by.barcode <- c()

  # assert barcodes are present in object
  if(is.null(obj$barcode)){
    message("no barcode field identified in object metadata. Stopping")
    break()
  }

  # identify cells that have more than 2 barcodes for automatic filtering
  barcodes <- str_split(obj$barcode, pattern = ";", simplify = F)
  selection <- which(lapply(barcodes, length) > 2)
  selection <- rownames(obj@meta.data[selection,])

  # filter cells that have 1 or more than two barcodes
  barcodes <- as.data.frame(str_split(obj$barcode, pattern = ";", simplify = T), row.names = rownames(obj@meta.data))[,c(1,2)]
  barcodes$V2 <- gsub(barcodes$V2, pattern = "^$|^ $", replacement = NA)
  barcodes <- na.omit(barcodes)

  # find cells that have unique combinations of barcodes in any order
  # these are likely doublets that can be filtered
  for (i in 1:nrow(barcodes)){
    cell.fwd <- paste(barcodes[i,]$V1, barcodes[i,]$V2, sep = ";")
    cell.rev <- paste(barcodes[i,]$V2, barcodes[i,]$V1, sep = ";")
    same.cells <- which(obj$barcode == cell.fwd | obj$barcode == cell.rev)
    if (length(obj$barcode[same.cells]) == 1){
      cell.to.filter <- names(obj$barcode[same.cells])
      doublets.by.barcode <- c(doublets.by.barcode, cell.to.filter)
    }
  }
  # add cells with 3 or more barcodes to the other doublets
  doublets.by.barcode <- c(doublets.by.barcode, selection)
  return(doublets.by.barcode)
}

# Enids method - way better than mine
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



runClusterProfiler <- function(dge, lfc.threshold = 0.2, padj.threshold = 0.1, OrgDb = org.Mm.eg.db, category = "BP", sample = "gene-set", outdir = NULL){
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
  gene.df.up <- bitr(geneset.up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  gene.df.dn <- bitr(geneset.dn, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)

  # get enriched GO terms
  ego.up <- enrichGO(gene = gene.df.up$ENTREZID,
                     OrgDb = OrgDb, ont = category,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  head(ego.up, 20)

  ego.dn <- enrichGO(gene = gene.df.dn$ENTREZID,
                     OrgDb = OrgDb, ont = category,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  head(ego.dn, 20)

  print(barplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
  print(barplot(ego.dn, showCategory=15) + ggtitle(paste("Downregulated", category, "GO terms in", sample)))

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
    print(barplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
    print(barplot(ego.dn, showCategory=15) + ggtitle(paste("Downregulated", category, "GO terms in", sample)))
    print(clusterProfiler::emapplot(ego.up, showCategory=15) + ggtitle(paste("Upregulated", category, "GO terms in", sample)))
    print(clusterProfiler::emapplot(ego.dn, showCategory=15)+ ggtitle(paste("Downregulated", category, "GO terms in", sample)))
    dev.off()
  }
}

findDoubletsBySimulation <- function(obj, ntop = 1000, mads = 3, seed = 10101){
  message("Running doublet detection on input single cell object using scds and scran methods")
  require(scds)
  require(scater)
  require(rsvd)
  require(Rtsne)
  require(cowplot)
  require(BiocSingular)
  require(scran)

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

  # scds doublet detection
  message("")
  message("#- Annotate doublets using co-expression based doublet scoring:")
  sample.sce <- cxds(sce = sample.sce, retRes = T, verb = T, estNdbl = T, ntop = ntop)

  message("")
  message("#- Annotate doublets using simulation approach:")
  sample.sce <- bcds(sample.sce, retRes = T, verb = T, ntop = ntop, varImp = T, estNdbl = T)

  message("")
  message("#- Annotate doublets using hybrid co-expression and simulation approach:")
  sample.sce <- cxds_bcds_hybrid(sample.sce, verb = T, estNdbl = T)

  # scran simulation approach
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

# edge R for 2 groups Ctrl vs Test
run_edgeR_LRT <- function(counts, group, plots = F){
  suppressPackageStartupMessages(require(edgeR))
  # this is a simple A vs B test. See other functions for more complex GLMs
  # adapted from
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRLRT.R
  message("Running edgeR LRT")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(counts, group = group)
    message("Calculate normalisation factors")
    dge <- calcNormFactors(dge)
    message("Generate model matrix")
    design <- model.matrix(~ group)
    print(design)
    message("Estimate dispersions")
    dge <- estimateDisp(dge, design = design)
    message("Fit GLM model")
    fit <- glmFit(dge, design = design)
    message("Likelihood ratio test")
    lrt <- glmLRT(fit)
    message("Export top tags")
    tt <- topTags(lrt, n = Inf)
  })

  # plots
  if(isTRUE(plots)){
    message("Plotting, this could take a while...")
    plotBCV(dge)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(group)), pch = 19)
    plotSmear(lrt)
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
