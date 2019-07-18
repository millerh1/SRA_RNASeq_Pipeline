DGE_pipeline <- function(
  countDir,
  # countFileType = c("salmon"),
  colData,
  tx2gene,
  contrast, 
  outdir,
  toRun = c("Volcano", "circos", "GSEA", "Heatmap", "plotCounts", "dimReduce"),
  TERM2GENE = NULL,
  species = "human",
  mart.res = NULL,
  numTopHeat = 30,
  numTopPlotCount = 50,
  numCircosGenes = 2000,
  heatLegend = FALSE,
  dimReduceMethod = c("PCA", "distanceHeatmap"),
  returnDataOnly = FALSE,
  mart.chr = NULL,
  extraAnnoRows = NULL,
  dds = NULL
  ) {
  
  # # # Bug testing
  # # dds <- dds
  # load("Data/tx2gene.RData")
  # countDir <- "Results/Salmon.out/"
  # colData <- as.data.frame(colData(ddsEvH))
  # contrast <- c("category", "EWS", "MSC")
  # outdir <- "Results/DGE/EWS_vs_Healthy/EWS_vs_MSC/"
  # tx2gene <- tx2gene
  # pathwaysToAdd <- NULL
  # load("Data/hsapiens_simple_TERM2GENE.rda")
  # TERM2GENE <- hsapiens_simple_TERM2GENE
  # toRun <- c("Volcano", "GSEA", "Heatmap")
  
  # Libraries required
  require(EnhancedVolcano)
  require(clusterProfiler)
  require(pheatmap)
  require(circlize)
  require(data.table)
  require(ggpubr)
  require(biomaRt)
  require(RColorBrewer)
  require(gplots)
  require(DESeq2)
  
  if (! dir.exists(outdir) & ! returnDataOnly) {
    dir.create(outdir)
  }
  
  resList <- list()
  
  
  # Read in data
  contrastCol <- colnames(colData)[which(colnames(colData) == contrast[1])]
  # Wrangle colData
  colData <- colData[which(colData[, contrastCol] %in% c(contrast[2], contrast[3])),]
  colData[, contrastCol] <- factor(colData[, contrastCol], levels = c(contrast[3], contrast[2]))
  
  if (is.null(dds)) {
    print("Loading count data ... ")
    
    # Get count files
    samples <- rownames(colData)
    countFiles <- file.path(countDir, samples, "quant.sf")
    names(countFiles) <- samples
    # Load with Tximport
    txi <- tximport(files = countFiles, 
                    type="salmon", tx2gene = tx2gene,
                    ignoreTxVersion = T)
    # Make DESeq dataset
    dds <- DESeqDataSetFromTximport(txi = txi, colData = colData, 
                                    design = formula(paste0("~", contrastCol)))
    mn <- min(table(colData[, contrastCol]))
    keep <- rowSums(counts(dds)) >= mn
    dds <- dds[keep,]
    print("DESeq analysis ...")
    timestamp()
    dds <- DESeq(dds)
    print("DONE")
    timestamp()
    
    if (! returnDataOnly) {
      print("Saving dds object...")
      save(dds, file = file.path(outdir, "dds.RData"))
    }
    resList[["dds"]] <- dds
  } else {
    # Does DESeq need to be run still on user-supplied dds?
    gg <- resultsNames(dds)
    if (! length(gg)) {
      print("No results found in dds object -- running DESeq2...")
      timestamp()
      dds <- DESeq(dds)
      print("DONE")
      timestamp()
    }
    
    if (! returnDataOnly) {
      print("Saving dds object...")
      save(dds, file = file.path(outdir, "dds.RData"))
    }
    resList[["dds"]] <- dds
    
  }
  
  # Downstream of DESeq
  print("Gathering results from dds object ... ")
  #load("Results/DGE/EWS_vs_EWS/CHLA10_vs_CHLA9/dds.RData")
  # contrast <- c("Cell", "CHLA10", "CHLA9")
  res <- results(dds, contrast = contrast)
  # Change dds for plotting purposes
  dds <- dds[,which(colData(dds)[,contrastCol] %in% c(contrast[2], contrast[3]))]
  colData(dds)[,contrastCol] <- droplevels(colData(dds)[,contrastCol])
  
  resdf <- as.data.frame(res)
  resdf$geneName <- rownames(resdf)
  resdf <- resdf[,c(7,1:6)]
  resdf$padj[which(resdf$padj == 0)] <- .Machine$double.xmin
  resdf$GSEA <- -log10(resdf$padj) * sign(resdf$log2FoldChange)
  
  titleStr <- paste0(contrast[2], " vs. ", contrast[3])
  print("DONE")
  
  
  # Volcano
  if ("Volcano" %in% toRun) {
    print("Generating volcano plot ... ")
    pord <- resdf$padj
    pord <- -log10(pord)
    pord <- pord[order(pord, decreasing = T)]
    pval <- pord[150]
    ford <- resdf$log2FoldChange
    ford <- abs(ford)
    ford <- ford[order(ford, decreasing = T)]
    fval <- ford[600]
    ev <- EnhancedVolcano(resdf, lab = resdf$geneName, x = "log2FoldChange", 
                    title = paste0(titleStr, " Volcano Plot"),
                    y = "padj", pCutoff = 10^-pval, FCcutoff = fval)
    if (! returnDataOnly) {
      ggsave(ev, filename = file.path(outdir, "Volcano.png"), height = 7, width = 8)
    }
    resList[["Volcano"]] <- ev
    print("DONE")
  }
  
  # Circos
  if ("circos" %in% toRun) {
    print("Generating circos plot ... ")
    cytopath <- file.path(tempdir(), 'ideoband.txt')
    cytopath2 <- paste0(cytopath, ".gz")
    if (species == "human") {
      cytobandURL <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz"
    }
    cmd <- paste0("wget -O ", cytopath2, " ", cytobandURL, " && gunzip ", cytopath2)
    system(cmd)
    cytoBand <- cytopath
    resdf2 <- resdf[, c(1, 8)]
    resdf2 <- unique(resdf2)
    pord <- abs(resdf2$GSEA)
    pord <- pord[order(pord, decreasing = T)]
    pval <- pord[numCircosGenes]
    
    res.df.sig <- resdf2[which(abs(resdf2$GSEA) >= pval),]

    if (is.null(mart.chr)) {
      ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
      mart.chr <- getBM(attributes = c("chromosome_name", "start_position", "end_position",
                                       "strand",
                                       "external_gene_name"),
                        mart = ensembl)
    }
    
    
    good.seq <- c(1:22, "X", "Y")
    mart.chr <- mart.chr[which(mart.chr$chromosome_name %in% good.seq),]
    mart.chr$chromosome_name <- paste0("chr", mart.chr$chromosome_name)
    colnames(mart.chr) <- c("seqname", "start", "end", "strand", "geneName")
    
    res.df.sig <- merge(x = mart.chr, y = res.df.sig, by = "geneName")
    res.df.sig2 <- res.df.sig[,c(2:4, 1, 6)]
    res.df.sig2 <- res.df.sig2[order(res.df.sig2$seqname, res.df.sig2$start),]
    qt <- quantile(abs(res.df.sig2$GSEA))
    topCol <- as.numeric(qt["75%"]) # use abs quartiles to determine color scale
    # topCol <- max(abs(res.df.sig2$GSEA)) * .85
    col_fun <- colorRamp2(c((-1 * topCol), 0, topCol), c("green", "black", "red"))
    
    resList[["circosData"]] <- res.df.sig2
    png(file.path(outdir, "DGE_circos.png"), height = 12, width = 12, units = c("in"), res = 600)
    circos.par(track.margin = c(.01, .01))
    circos.initializeWithIdeogram(cytoband = cytoBand, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
    circos.genomicHeatmap(res.df.sig2, col = col_fun, numeric.column = 5)
    title(main = paste0("DGE circos ", titleStr))
    mtext(paste0("Top ", numCircosGenes, " significant genes"), side = 1)
    dev.off()
  }
  
  # Heatmap
  if ("Heatmap" %in% toRun) {
    print("Generating heatmap ... ")
    cts <- counts(dds, normalized = T)
    colData <- colData(dds)
    choose <- resdf[order(resdf$GSEA, decreasing = T),]
    choose <- choose[which(! is.na(choose$GSEA)),]
    n <- length(choose$geneName)
    l <- n-floor(numTopHeat/2)
    choose <- choose[c(1:floor(numTopHeat/2), l:n),]
    genes <- unique(as.character(choose$geneName))
    cts <- cts[genes,]
    colData <- as.data.frame(colData[which(as.character(colData[, contrast[1]]) %in% c(contrast[2], contrast[3])),])
    colDataPlt <- colData[, contrastCol, drop = F]
    if (! is.null(extraAnnoRows)) {
      colDataPlt <- colData[, c(contrastCol, extraAnnoRows), drop = F]
    }
    cts <- cts[,which(colnames(cts) %in% rownames(colDataPlt))]
    cts2 <- log2(cts + 1)
    cols <- c("lightcoral", "lightblue")
    names(cols) <- c(contrast[2], contrast[3])
    ph <- pheatmap(mat = cts2, annotation_col = colDataPlt, 
                   scale = "row", color = greenred(100), legend = heatLegend,
                   show_colnames = F, main = paste0(titleStr, " DGE Heat Map\n") )
    if (! returnDataOnly) {
      ggsave(ph, filename = file.path(outdir, "Heatmap.png"), height = 8)
    }
    resList[["Heatmap"]] <- ph
    print("DONE")
  }
  
  # GSEA
  if ("GSEA" %in% toRun) {
    print("GSEA analysis ...")
    ranks <- resdf$GSEA
    names(ranks) <- resdf$geneName
    ranks <- ranks[which(! is.na(ranks))]
    ranks <- ranks[which(! duplicated(names(ranks)))]
    ranks <- ranks[order(ranks, decreasing = T)]
    GSEA <- myGSEA(ranks = ranks, returnDataOnly = returnDataOnly,
                   TERM2GENE = TERM2GENE, 
                   outDir = outdir, 
                   plotFile = "GSEA", 
                   Condition = titleStr)
    resList[["GSEA"]] <- GSEA
    print("DONE")
  }
  
  # plotCounts
  if ("plotCounts" %in% toRun) {
    print("Plotting counts for top DGEs...")
    resdf2 <- resdf[order(resdf$GSEA, decreasing = T),]
    m <- floor(numTopPlotCount/2)
    n <- length(resdf2$geneName)
    geneToPlot <- resdf2[c(1:m, (n-m+1):n ) ,]
    k <- length(geneToPlot$geneName)
    plotCounts <- list()
    if (! returnDataOnly) {
      dir.create(file.path(outdir, "plotCounts"))
    }
    for (i in 1:k) {
      cat("\n", i, " of ", k)
      gene <- geneToPlot$geneName[i]
      subTitle <- paste0("FDR: ", signif(geneToPlot$padj[i], 3), 
                         "\nLog2FoldChange: ", signif(geneToPlot$log2FoldChange[i], 3))
      # counts <- plotCounts(dds, gene = gene, intgroup = "condition", returnData = T)
      counts <- plotCounts(dds, gene = gene, intgroup = contrast[1], returnData = T)
      colnames(counts)[2] <- "condition"
      counts$count <- log2(counts$count + 1)
      gbp <- ggboxplot(data = counts, x = "condition", y = "count", add = "jitter", 
                       ylab = "Normalized read counts", title = gene, subtitle = subTitle,
                       xlab = FALSE, fill = "condition") + rremove("legend")
      plotCounts[[i]] <- gbp
      names(plotCounts)[i] <- gene
      if (! returnDataOnly) {
        ggsave(gbp, filename = file.path(outdir, "plotCounts", paste0(gene, ".png")))
      }
    }
    resList[["plotCounts"]] <- plotCounts
  }
  
  # dimReduce
  if ("dimReduce" %in% toRun) {
    if ("PCA" %in% dimReduceMethod) {
      print("Calculating dimension reduction by PCA ... ")
      vsd <- vst(dds)
      resList[["vsd"]] <- vsd
      pca <- plotPCA(vsd, intgroup = contrast[1])
      resList[["pca"]] <- pca
      if (! returnDataOnly) {
        ggsave(pca, filename = file.path(outdir, "PCA.png"))
      }
      print("DONE")
    } 
    if ("distanceHeatmap" %in% dimReduceMethod) {
      print("Calculating sample-to-sample distance from normalized read counts ... ")
      cts3 <- counts(dds, normalized = T)
      sampleDists <- dist(t(log2(cts3+1)))
      sampleDistMatrix <- as.matrix(sampleDists)
      colnames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      anns <- as.data.frame(colData(dds))
      # annoCol <- anns[,"condition", drop = F]
      annoCol <- anns[,contrast[1], drop = F]
      if (! is.null(extraAnnoRows)) {
        print("Adding extra columns to heatmap")
        print(extraAnnoRows)
        
        annoCol <- anns[,c(contrast[1], extraAnnoRows), drop = F]
        print(annoCol)
      } 
      ph <- pheatmap(sampleDistMatrix, annotation_row = annoCol,
                     clustering_distance_rows=sampleDists, show_rownames = F,
                     main = "Dim Reduction plot of gene expression\n",
                     clustering_distance_cols=sampleDists,
                     col=colors)
      
      resList[["dimReducePlot"]] <- ph
      resList[["dimReduceDistData"]] <- sampleDistMatrix
      resList[["dimReduceAnno"]] <- annoCol
      
      if (! returnDataOnly) {
        ggsave(ph, filename = file.path(outdir, "dimReducePlot.png"), height = 7, width = 7)
      }
      print("DONE")
    }
    
    
  }
  
  # Final output
  resdf <- resdf[which(! is.na(resdf$GSEA)),]
  cat("\nReturning results...\n")
  resdf <- merge(x = mart.res, y = resdf, by = "geneName")
  resdf <- resdf[order(resdf$padj),]
  if (! returnDataOnly) {
    fwrite(resdf, file = file.path(outdir, "DGE_results.csv"), sep = ",")
  }
  
  resList[["DGE_Result"]] <- resdf
  
  return(resList)
  
}



DTEE_pipeline <- function(
  colData,
  TEcountFolder,
  outdir,
  contrast,
  dds = NULL, # Supply dds object to skip DESeq() step
  CTSobj = NULL, # If dds is NULL -- supply CTSobj to skip obtaining count matrix from cntTable files
  toRun = c("Heatmap", "Volcano", "plotCounts", "dimReduce"),
  ens2symb = NULL,
  compRanksGSEA = NULL,
  numTopHeat = 30,
  numTopPlotCount = 50,
  returnDataOnly = FALSE,
  extraAnnoRows = NULL # Choose extra annotation rows to add to heatmap
) {
  # colData <- as.data.frame(ddsEvH@colData)
  # TEcountFolder <- "Results/TEcount.out"
  # outdir <- "Results/DTEE/EWS_vs_Healthy/EWS_vs_MSC"
  # contrast <- c("category", "EWS", "MSC")
  # compRanksGSEA <- GSEA_comp 
  # ens2symb <- NULL
  # toRun <- c("Heatmap", "Volcano", "plotCounts", "dimReduce")
  # CTSobj <- NULL
  # # dds <- NULL
  # numTopHeat = 30
  # numTopPlotCount = 50
  # returnDataOnly <- F
  
  # PIPELINE START 
  # Libraries
  require(DESeq2)
  require(EnhancedVolcano)
  require(ggpubr)
  require(biomaRt)
  require(VennDiagram)
  require(pheatmap)
  require(RColorBrewer)
  require(gplots)
  
  # Initialize results object
  resList <- list()
  
  # Make outdir
  if (! dir.exists(outdir) & ! returnDataOnly) {
    dir.create(outdir)
  }
  
  # Wrangle colData
  contrastCol <- colnames(colData)[which(colnames(colData) == contrast[1])]
  # Wrangle colData
  colData <- colData[which(colData[, contrastCol] %in% c(contrast[2], contrast[3])),]
  colData[, contrastCol] <- factor(colData[, contrastCol], levels = c(contrast[3], contrast[2]))
  
  # Get ens2symb object if does not exist
  if (is.null(ens2symb)) {
    print("No ens2symb object -- gathering from biomaRt...")
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens2symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
    
    print("DONE")
  } 
  
  colnames(ens2symb) <- c("gene.TE", "geneName")
  
  if (is.null(dds)) {
    # Make sure files exist
    files_cnt <- file.path(TEcountFolder, paste0(rownames(colData), ".cntTable"))
    names(files_cnt) <- rownames(colData)
    if(! all(file.exists(files_cnt))) {
      filesBad <- files_cnt[which(! file.exists(files_cnt))]
      msg <- paste0("File does not exist: ", filesBad)
      stop(msg)
    }
    
    
    # Get count data for genes and TEs
    if (is.null(CTSobj)) {
      print("Assembling TE count dataframe...")
      n <- length(files_cnt)
      for (i in 1:n) {
        cat("\n", i, " of ", n)
        sampleName <- names(files_cnt)[i]
        fileFrame <- read.table(files_cnt[i], header = T, stringsAsFactors = F)
        
        # Get ensembl genes and convert to gene symbols
        geneFrame <- fileFrame[which(substr(fileFrame$gene.TE, 1, 4) == "ENSG"),]
        geneFrame$gene.TE <- substr(geneFrame$gene.TE, 1, 15)
        
        geneFrame <- merge(x = ens2symb, y = geneFrame, by = "gene.TE")
        geneFrame <- geneFrame[,c(2, 3)]
        geneFrame <- unique(geneFrame)
        geneFrame2 <- aggregate(geneFrame[,2], 
                                by = list(geneFrame$geneName), FUN = sum)
        colnames(geneFrame2) <- c("geneName", sampleName)
        
        # Get TE frame
        teFrame <- fileFrame[which(substr(fileFrame$gene.TE, 1, 4) != "ENSG"),]
        colnames(teFrame) <- c("teName", sampleName)
        
        if (i == 1) {
          geneCTS <- geneFrame2
          teCTS <- teFrame
        } else {
          geneCTS <- merge(x = geneCTS, y = geneFrame2, by = "geneName")
          teCTS <- merge(x = teCTS, y = teFrame, by = "teName")
        }
      }
      rownames(geneCTS) <- geneCTS$geneName
      geneCTS <- geneCTS[,c(-1)]
      rownames(teCTS) <- teCTS$teName
      teCTS <- teCTS[,c(-1)]
      CTS <- rbind(geneCTS, teCTS)
      if (! returnDataOnly) {
        print("Saving count data..")
        save(CTS, file = file.path(outdir, "CTS.RData"))
      }
      
      print("DONE")
    } else {
      print("Loading count data from CTS_list object...")
      CTS <- CTSobj
      print("DONE")
    }
    
    # DESeq2
    dds <- DESeqDataSetFromMatrix(countData = CTS,
                                  colData = colData,
                                  design = formula(paste0("~", contrastCol)))
    mn <- min(table(colData[,contrastCol]))
    keep <- rowSums(counts(dds)) >= mn
    cat(paste0("\n",length(which(keep)), " genes passed pre-filtering step!\n"))
    dds <- dds[keep,]
    print("DESeq analysis ...")
    timestamp()
    dds <- DESeq(dds)
    print("DONE")
    timestamp()
    resList[["dds"]] <- dds
    if (! returnDataOnly) {
      save(dds, file = file.path(outdir, "dds.RData"))
    }
  } else {
    print("Loading data from dds object...")
    resList[["dds"]] <- dds
    if (! returnDataOnly) {
      save(dds, file = file.path(outdir, "dds.RData"))
    }
  }
  
  
  print("Gathering DESeq results ... ")
  # res <- DESeq2::results(dds, contrast = c("condition", contrast[2], contrast[3]))
  res <- DESeq2::results(dds, contrast = contrast)
  resdf <- as.data.frame(res)
  resdf$geneName <- rownames(resdf)
  resdf <- resdf[,c(7,1:6)]
  resdf$padj[which(resdf$padj == 0)] <- .Machine$double.xmin
  resdf$GSEA <- -log10(resdf$padj) * sign(resdf$log2FoldChange)
  resdfGene <- resdf[which(resdf$geneName %in% ens2symb$geneName),]
  if (! is.null(compRanksGSEA)) {
    print("Testing validity of TEcount DGE method against user-supplied ranks...")
    ranks <- resdfGene$GSEA
    names(ranks) <- resdfGene$geneName
    ranks <- ranks[which(! is.na(ranks))]
    ranks <- ranks[which(! duplicated(names(ranks)))]
    ranks <- ranks[order(ranks, decreasing = T)]
    oList <- list("TEcount Genes" = names(ranks), "Comparison Genes" =  names(compRanksGSEA))
    
    x <- calculate.overlap(oList)
    if (! returnDataOnly) {
      venn.diagram(x = oList, filename = file.path(outdir, "TEcount_genes_compared.png"), margin = .05,
                   fill = c("firebrick", "lightblue"))
    }
    names(x) <- c(names(oList), c("Overlap"))
    resList[["overlapList"]] <- x
    ranks2 <- ranks[which(names(ranks) %in% x$Overlap)]
    compRanksGSEA2 <- compRanksGSEA[which(names(compRanksGSEA) %in% x$Overlap)]
    ranks2 <- ranks2[order(names(ranks2))]
    compRanksGSEA2 <- compRanksGSEA2[order(names(compRanksGSEA2))]
    all(names(ranks2) == names(compRanksGSEA2)) # Should be TRUE
    resList[["corTest"]] <- cor.test(x = ranks2, y = compRanksGSEA2)
    df <- data.frame(row.names = names(ranks2), "TEcountRanks" = ranks2, "ComparisonRanks" = compRanksGSEA2)
    gsp <- ggscatter(data = df, x = "TEcountRanks", 
                     y = "ComparisonRanks", cor.coef= T, 
                     add = "reg.line", conf.int = T, title = "TEcount vs comparison method DGE",
                     add.params = list(color = "firebrick", fill = "lightgray"))
    resList[["corPlot"]] <- gsp
    if (! returnDataOnly) {
      ggsave(gsp, filename = file.path(outdir, "DGE_with_TEcount_vs_comparison_method.png"))
    }
    print("DONE")
  }
  
  print("Analyzing trends in TE expression...")
  
  resdf <- resdf[which(! resdf$geneName %in% ens2symb$geneName),] # only TE results
  titleStr <- paste0(contrast[2], " vs. ", contrast[3]) # title string for plots
  
  # Create families etc
  teStr <- resdf$geneName
  teStr <- strsplit(teStr, split = ":")
  resdf$teGene <- sapply(teStr, "[[", 1)
  resdf$teFamily <- sapply(teStr, "[[", 2)
  resdf$teClass <- sapply(teStr, "[[", 3)
  map <- data.frame(resdf[,c("geneName", "teGene", "teFamily", "teClass")])
  resdf_sig_up <- resdf[which(resdf$GSEA > 1.3),]
  resdf_sig_dn <- resdf[which(resdf$GSEA < -1.3),]
  
  # Plot differential families and classes
  df <- as.data.frame(table(resdf_sig_up$teFamily))
  df$Freq <- (df$Freq/sum(df$Freq))*100
  df$Group <- "Overexpressed"
  df2 <- as.data.frame(table(resdf_sig_dn$teFamily))
  df2$Freq <- (df2$Freq/sum(df2$Freq))*100
  df2$Group <- "Underexpressed"
  df3 <- rbind(df, df2)
  plt <- ggbarplot(data = df3, x = "Var1", y = "Freq", fill = "Group",
                   color = "Group",  palette = "Paired", position = position_dodge(0.9), 
                   xlab = "TE family", ylab = "Percentage of differential TEs\n", 
                   title = paste0(titleStr,  " Differential TE families\n")) + rotate_x_text()
  resList[["DTEE_families"]] <- plt
  df <- as.data.frame(table(resdf_sig_up$teClass))
  df$Freq <- (df$Freq/sum(df$Freq))*100
  df$Group <- "Overexpressed"
  df2 <- as.data.frame(table(resdf_sig_dn$teClass))
  df2$Freq <- (df2$Freq/sum(df2$Freq))*100
  df2$Group <- "Underexpressed"
  df3 <- rbind(df, df2)
  plt2 <- ggbarplot(data = df3, x = "Var1", y = "Freq", fill = "Group",
                    color = "Group",  palette = "Paired", position = position_dodge(0.9), 
                    xlab = "TE class", ylab = "Percentage of differential TEs\n", 
                    title = paste0(titleStr,  " Differential TE classes\n")) + rotate_x_text()
  resList[["DTEE_classes"]] <- plt2
  
  if (! returnDataOnly) {
    pltt <- ggarrange(plt, plt2)
    ggsave(pltt, filename = file.path(outdir, "DTEE_families_and_classes.png"), height = 6, width = 12)
  }
  
  print("DONE")
  print("Running additional downstream analyses...")
  
  # Volcano
  if ("Volcano" %in% toRun) {
    print("Generating volcano plot...")
    m <- length(resdf$geneName)
    pord <- resdf$padj
    pord <- -log10(pord)
    pord <- pord[order(pord, decreasing = T)]
    pval <- pord[floor(m*.1)]
    ford <- resdf$log2FoldChange
    ford <- abs(ford)
    ford <- ford[order(ford, decreasing = T)]
    fval <- ford[floor(m*.2)]
    ev <- EnhancedVolcano(resdf, lab = resdf$teGene, x = "log2FoldChange", 
                          title = paste0(titleStr, " DTEE Volcano Plot"),
                          y = "padj", pCutoff = 10^-pval, FCcutoff = fval)
    if (! returnDataOnly) {
      ggsave(ev, filename = file.path(outdir, "Volcano.png"))
    }
    resList[["Volcano"]] <- ev
    print("DONE")
  }
  
  # Heatmap
  if ("Heatmap" %in% toRun) {
    print("Creating heatmap...")
    cts <- counts(dds, normalized = T)
    colData <- colData(dds)
    choose <- resdf[order(resdf$GSEA, decreasing = T),]
    choose <- choose[which(! is.na(choose$GSEA)),]
    n <- length(choose$geneName)
    l <- n-(numTopHeat/2)
    choose <- choose[c(1:(numTopHeat/2), l:n),]
    genes <- unique(as.character(choose$geneName))
    cts <- cts[genes,]
    # colData <- as.data.frame(colData[which(as.character(colData[, "condition"]) %in% c(contrast[2], contrast[3])),])
    colData <- as.data.frame(colData[which(as.character(colData[, contrast[1]]) %in% c(contrast[2], contrast[3])),])
    # colDataPlt <- colData[,"condition", drop = F]
    
    colDataPlt <- colData[,contrast[1], drop = F]
    if (! is.null(extraAnnoRows)) {
      colDataPlt <- colData[,c(contrast[1], extraAnnoRows), drop = F]
    }
    
    cts <- cts[,which(colnames(cts) %in% rownames(colDataPlt))]
    cts2 <- log2(cts + 1)
    cols <- c("lightcoral", "lightblue")
    names(cols) <- c(contrast[2], contrast[3])
    rowDataPltraw <- data.frame(row.names = rownames(cts2), geneName = rownames(cts2))
    rowDataPltraw <- merge(x = rowDataPltraw, y = map, by = "geneName")
    rownames(rowDataPltraw) <- rowDataPltraw$geneName
    rowDataPlt <- rowDataPltraw[,"teClass", drop = F]
    rowDataPlt <- rowDataPlt[order(match(rownames(rowDataPlt), rownames(cts2))),, drop = F]
    rowDataPltraw <- rowDataPltraw[order(match(rownames(rowDataPltraw), rownames(cts2))),, drop = F]
    
    ph <- pheatmap(mat = cts2, annotation_col = colDataPlt, 
                   annotation_row = rowDataPlt, labels_row = rowDataPltraw$teGene,
                   scale = "row", color = greenred(100), legend = T,
                   annotation_colors = list(category = cols), 
                   show_colnames = F, main = paste0(titleStr, " DTEE Heat Map\n") )
    if (! returnDataOnly) {
      ggsave(ph, filename = file.path(outdir, "Heatmap.png"), height = 8)
    }
    resList[["Heatmap"]] <- ph
    print("DONE")
  }
  
  # plotCounts
  if ("plotCounts" %in% toRun) {
    print("Plotting counts for top DTEEs...")
    resdf2 <- resdf[order(resdf$GSEA, decreasing = T),]
    m <- floor(numTopPlotCount/2)
    n <- length(resdf2$geneName)
    geneToPlot <- resdf2[c(1:m, (n-m+1):n ) ,]
    k <- length(geneToPlot$geneName)
    plotCounts <- list()
    if (! returnDataOnly) {
      dir.create(file.path(outdir, "plotCounts"))
    }
    for (i in 1:k) {
      cat("\n", i, " of ", k)
      gene <- geneToPlot$geneName[i]
      geneTitle <- geneToPlot$teGene[i]
      capText <- paste0("Family: ", geneToPlot$teFamily[i], ". Class: ", geneToPlot$teClass[i], "")
      subTitle <- paste0("FDR: ", signif(geneToPlot$padj[i], 3), 
                         "\nLog2FoldChange: ", signif(geneToPlot$log2FoldChange[i], 3))
      # counts <- plotCounts(dds, gene = gene, intgroup = "condition", returnData = T)
      counts <- plotCounts(dds, gene = gene, intgroup = contrast[1], returnData = T)
      colnames(counts)[2] <- "condition"
      counts$count <- log2(counts$count + 1)
      gbp <- ggboxplot(data = counts, x = "condition", y = "count", add = "jitter", caption = capText,
                       ylab = "Normalized read counts", title = geneTitle, subtitle = subTitle,
                       xlab = FALSE, fill = "condition") + yscale("log10", .format = T) + rremove("legend")
      plotCounts[[i]] <- gbp
      names(plotCounts)[i] <- gene
      if (! returnDataOnly) {
        ggsave(gbp, filename = file.path(outdir, "plotCounts", paste0(geneTitle, ".png")))
      }
    }
    resList[["plotCounts"]] <- plotCounts
  }
  
  # dimReduce
  if ("dimReduce" %in% toRun) {
    print("Calculating dimension reduction ... ")
    cts3 <- counts(dds, normalized = T)
    cts3 <- cts3[which(rownames(cts3) %in% resdf$geneName),]
    sampleDists <- dist(t(log2(cts3+1)))
    
    sampleDistMatrix <- as.matrix(sampleDists)
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    anns <- as.data.frame(colData(dds))
    # annoCol <- anns[,"condition", drop = F]
    annoCol <- anns[,contrast[1], drop = F]
    if (! is.null(extraAnnoRows)) {
      annoCol <- anns[,c(contrast[1], extraAnnoRows), drop = F]
    }
    ph <- pheatmap(sampleDistMatrix, annotation_row = annoCol,
                   clustering_distance_rows=sampleDists, show_rownames = F,
                   main = "Dim Reduction plot of TE expression\n",
                   clustering_distance_cols=sampleDists,
                   col=colors)
    
    
    resList[["dimReducePlot"]] <- ph
    resList[["dimReduceData"]] <- sampleDistMatrix
    resList[["dimReduceAnno"]] <- annoCol
    
    if (! returnDataOnly) {
      ggsave(ph, filename = file.path(outdir, "dimReducePlot.png"), height = 7, width = 7)
    }
    print("DONE")
  }
  
  # Return final frame
  resdf <- resdf[,c(9, 10, 11, 2:8)]
  resdf <- resdf[order(resdf$padj),]
  resList[["DGE_results"]] <- resdf
  if (! returnDataOnly) {
    write.csv(resdf, file = file.path(outdir, "DTEE_results.csv"), row.names = F, quote = F)
  }
  
  return(resList)
  print("Pipeline finished!")
  timestamp()
}

# Helper for GSEA
myGSEA <- function(ranks, TERM2GENE, plotFile,
                   outDir, Condition = Condition,
                   padjustedCutoff = .05, returnDataOnly = FALSE) {
  
  resList <- list()
  
  ranks <- ranks[which(! is.na(ranks))]
  ranks <- ranks[order(ranks, decreasing = T)]
  EGMT <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                             nPerm = 1000, pvalueCutoff = padjustedCutoff)
  resGSEA <- as.data.frame(EGMT)
  
  resList[["EGMT"]] <- EGMT
  
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                               nPerm = 1000, pvalueCutoff = padjustedCutoff + .2)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                               nPerm = 1000, pvalueCutoff = padjustedCutoff + .45)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  resGSEA <- resGSEA[order(resGSEA$NES, decreasing = T),]
  topUP <- resGSEA$ID[1:10]
  resGSEA <- resGSEA[order(resGSEA$NES, decreasing = F),]
  topDOWN <- resGSEA$ID[1:10]
  resGSEA <- resGSEA[order(resGSEA$pvalue),]
  plUP <- list()
  plDOWN <- list()
  for ( i in 1:6 ) {
    pathway <- topUP[i]
    if (nchar(pathway) > 35) {
      pathTitle <- paste0(substr(pathway, 1, 30), "...")
    } else {
      pathTitle <- pathway
    }
    gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
    gp <- gp + ggplot2::labs(title = pathTitle,
                             subtitle = paste0("Enrichment score: ",
                                               round(resGSEA$NES[which(
                                                 resGSEA$ID == pathway
                                               )], 3)))
    gg <- ggplot2::theme_classic()
    gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                               plot.subtitle = gg[["plot.subtitle"]],
                               plot.margin = gg[["plot.margin"]])
    if ( i == 1) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
    } else if (i == 2) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
    } else if (i == 3) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
    } else if (i == 4) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
    } else if (i == 5) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
    } else if (i == 6) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
    }
    plUP[[i]] <- gp
    
    pathway <- topDOWN[i]
    if (nchar(pathway) > 35) {
      pathTitle <- paste0(substr(pathway, 1, 30), "...")
    } else {
      pathTitle <- pathway
    }
    gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
    gp <- gp + ggplot2::labs(title = pathTitle,
                             subtitle = paste0("Enrichment score: ",
                                               round(resGSEA$NES[which(
                                                 resGSEA$ID == pathway
                                               )], 3)))
    gg <- ggplot2::theme_classic()
    gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                               plot.subtitle = gg[["plot.subtitle"]],
                               plot.margin = gg[["plot.margin"]])
    if (i == 1) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
    } else if (i == 2) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
    } else if (i == 3) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
    } else if (i == 4) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
    } else if (i == 5) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
    } else if (i == 6) {
      gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
    }
    plDOWN[[i]] <- gp
    
  }
  
  gaUP <- ggpubr::ggarrange(plotlist = plUP, nrow = 2, ncol = 3)
  gaUP <- ggpubr::annotate_figure(gaUP,
                                  top = ggpubr::text_grob(paste0(
                                    "Top Over-Expressed Pathways in ", Condition
                                  ),
                                  size = 35)
  )
  resList[["GSEA_up"]] <- gaUP
  if (! returnDataOnly) {
    ggplot2::ggsave(plot = gaUP,
                    filename = file.path(outDir, paste0(plotFile, "_topPathwaysUP.png")),
                    height = 14, width = 20)
  }
  
  
  gaDOWN <- ggpubr::ggarrange(plotlist = plDOWN, nrow = 2, ncol = 3)
  gaDOWN <- ggpubr::annotate_figure(gaDOWN,
                                    top = ggpubr::text_grob(paste0(
                                      "Top Under-Expressed Pathways in ", Condition
                                    ),
                                    size = 35)
  )
  
  if (! returnDataOnly) {
    ggplot2::ggsave(plot = gaDOWN,
                    filename = file.path(outDir, paste0(plotFile, "_topPathwaysDOWN.png")),
                    height = 14, width = 20)
  }
  resList[["GSEA_down"]] <- gaDOWN
  
  if (! returnDataOnly) {
    data.table::fwrite(x = resGSEA, file = file.path(outDir,
                                                     paste0(plotFile,
                                                            "_GSEA.csv")))
  }
  
  resList[["eres"]] <- resGSEA
  
  cat("Returning ... ", names(resList))
  
  return(resList)
}


DTU_pipeline <- function(
  files,
  colData,
  contrast,
  outdir,
  tx2gene,
  numTopGenesForPathwayAnalysis = NULL,
  doPathwayAnalysis = F,
  TERM2GENE = NULL,
  returnDataOnly = F,
  doPlotProportions = F,
  numPlots = NULL,
  pvalueCutoffPathwayAnalysis = .1,
  DTU_obj = NULL
) {
  
  # # Bug testing
  # files = files
  # tx2gene = tx2gene
  # colData = colData
  # contrast = c("category", "EWS", "MSC")
  # outdir = "Results/DTU/EWS_vs_Healthy/EWS_vs_MSC3"
  
  # Libraries required
  require(DESeq2)
  require(DRIMSeq)
  require(stageR)
  require(ggplot2)
  require(clusterProfiler)
  require(data.table)
  
  if (! dir.exists(outdir) & ! returnDataOnly) {
    dir.create(outdir)
  }
  
  
  
  if (is.null(DTU_obj)) {
    print("Loading data")
    timestamp()
    
    # Fix the problem of incongruous row names
    txs <- tximport(files = files[1], 
                    type="salmon",  txOut = T,
                    ignoreTxVersion = T, ignoreAfterBar = T)
    txsd <- rownames(txs$counts)
    txsd <- strsplit(txsd, split = "\\.")
    txsd <- sapply(txsd, "[[", 1)
    tx2gene <- tx2gene[which(tx2gene$ensembl_transcript_id %in% txsd),]
    txsdf <- data.frame(ensembl_transcript_id = txsd, external_gene_name = txsd)
    tx2gene2 <- rbind(tx2gene, txsdf)
    tx2gene2 <- unique(tx2gene2)
    
    # Read in data
    txi <- tximport(files = files, countsFromAbundance = "dtuScaledTPM",
                    type="salmon", tx2gene = tx2gene2, txOut = T,
                    ignoreTxVersion = T)
    
    # Wrangle count matrix
    cts <- txi$counts
    cts <- cts[rowSums(cts) > 0,]
    cts <- as.data.frame(cts)
    cts$feature_id <- rownames(cts)
    rownames(cts) <- NULL
    fid <- strsplit(cts$feature_id, split = "\\.")
    cts$feature_id <- sapply(fid, "[[", 1)
    tx2gene3 <- tx2gene[,c(2,1)]
    colnames(tx2gene3) <- c("gene_id", "feature_id")
    cts <- merge(cts, tx2gene3, by = "feature_id")
    n <- length(colnames(cts))
    cts <- cts[,c(n, 1:(n-1))]
    
    # Wrangle colData
    samples <- data.frame(sample_id = rownames(colData), condition = colData[,which(colnames(colData) == contrast[1])])
    samples <- samples[which(samples$condition %in% c(contrast[2], contrast[3])),]
    samples$condition <- factor(samples$condition, levels = c(contrast[3], contrast[2]))
    
    print("Filtering counts")
    
    # Build dmDSdata object and filter
    d <- dmDSdata(counts = cts, samples = samples)
    n <- length(samples$sample_id)
    n.small <- min(table(samples$condition))
    d <- dmFilter(d,
                  min_samps_feature_expr=n.small,
                  min_samps_feature_prop=n.small, min_feature_prop=0.05,
                  min_samps_gene_expr=n, run_gene_twice = T)
    
    design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
    coef <- paste0("condition", contrast[2])
    
    print("Starting DRIMseq")
    timestamp()
    
    # Run DRIMSeq
    set.seed(1)
    d <- dmPrecision(d, design=design_full)
    d <- dmFit(d, design=design_full)
    d <- dmTest(d, coef=coef)
    if ( ! returnDataOnly ) {
      save(d, file = file.path(outdir, "DRIMSeq_out.RData"))
      
    }
    print("DONE")
    timestamp()
  } else {
    print("Loading DTU object")
    d <- DTU_obj
    if ( ! returnDataOnly ) {
      save(d, file = file.path(outdir, "DRIMSeq_out.RData"))
      
    }
  }
  
  
  print("Gathering results -- post-hoc testing...")
  # Get results, clean, and pos-hoc filter by SD
  res <- DRIMSeq::results(d)
  res.txp <- DRIMSeq::results(d, level="feature")
  no.na <- function(x) ifelse(is.na(x), 1, x)
  res$pvalue <- no.na(res$pvalue)
  res.txp$pvalue <- no.na(res.txp$pvalue)
  res.txp$pvalue[which(res.txp$pvalue == 0)] <- .Machine$double.xmin
  res.txp$adj_pvalue[which(res.txp$adj_pvalue == 0)] <- .Machine$double.xmin
  res$pvalue[which(res$pvalue == 0)] <- .Machine$double.xmin
  res$adj_pvalue[which(res$adj_pvalue == 0)] <- .Machine$double.xmin
  smallProportionSD <- function(d, filter=0.1) {
    cts <- as.matrix(subset(DRIMSeq::counts(d), select=-c(gene_id, feature_id)))
    gene.cts <- rowsum(cts, DRIMSeq::counts(d)$gene_id)
    total.cts <- gene.cts[match(DRIMSeq::counts(d)$gene_id, rownames(gene.cts)),]
    props <- cts/total.cts
    propSD <- sqrt(rowVars(props))
    propSD < filter
  }
  filt <- smallProportionSD(d)
  res.txp$pvalue[filt] <- 1 
  res.txp$adj_pvalue[filt] <- 1
  
  # StageR
  pScreen <- res$pvalue
  names(pScreen) <- res$gene_id
  pConfirmation <- matrix(res.txp$pvalue, ncol=1)
  rownames(pConfirmation) <- res.txp$feature_id
  tx2geneNew <- res.txp[,c("feature_id", "gene_id")]
  tx2geneNew <- unique(tx2geneNew)
  stageRObj <- stageR::stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                                pScreenAdjusted=FALSE, tx2gene=tx2geneNew)
  object <- stageRObj
  pScreen <- getPScreen(object)
  pConfirmation <- getPConfirmation(object)
  pScreenAdjusted <- isPScreenAdjusted(object)
  tx2gene <- getTx2gene(object)
  stageAdjPValues <-
    .stageWiseTest(
      pScreen = pScreen,
      pConfirmation = pConfirmation,
      alpha = .05,
      method = "dtu",
      pScreenAdjusted = pScreenAdjusted,
      tx2gene = tx2gene
    )
  object@adjustedP <- stageAdjPValues[["pAdjStage"]]
  object@alphaAdjusted <-
    stageAdjPValues[["alphaAdjusted"]]
  object@method <- "dtu"
  object@alpha <- .05
  object@adjusted <- TRUE
  stageRObj <- object
  
  drim.padj <- stageRObj@adjustedP
  rws <- rownames(drim.padj)
  rws <- strsplit(rws, ":")
  genes <- sapply(rws, "[[", 1)
  transcripts <- sapply(rws, "[[", 2)
  drim.padj <- data.frame(geneName = genes, transcriptID = transcripts, stringsAsFactors = F,
                          gene_pAdj = drim.padj[,1], tx_pAdj = drim.padj[,2])
  rownames(drim.padj) <- NULL
  
  if (! returnDataOnly) {
    write.csv(drim.padj, file = file.path(outdir, "DRIM_padj_results.csv"))
  }
  
  print("DONE -- downstream plotting and analysis...")
  
  resList <- list("DTU_obj" = d, "DTU_results" = drim.padj)
  
  if (doPlotProportions) {
    print("Plotting proportions for top genes")
    if (is.null(numPlots)) {
      stop("Please specify number of plots if doPlotProportions == TRUE")
    } else {
      drim <- unique(drim.padj$geneName[order(drim.padj$gene_pAdj)])[1:numPlots]
      plotList <- list()
      for (i in 1:length(drim)) {
        gene <- drim[i]
        print(gene)
        pp <- plotProportions(x = d, gene_id = gene, group_variable = "condition")
        plotList[[i]] <- pp
        names(plotList)[i] <- gene
        if (! returnDataOnly) {
          dir.create(file.path(outdir, "plots"))
          ggsave(filename = file.path(outdir, "plots", paste0(gene, ".png")), plot = pp)
        }
      }
      resList[["plots"]] <- plotList
    }
    print("DONE")
  }
  
  if (doPathwayAnalysis) {
    print("Pathway enrichment of top genes")
    if (is.null(TERM2GENE)) {
      stop("Please specify TERM2GENE")
    } else {
      drim.padj2 <- drim.padj[which(drim.padj$gene_pAdj < .01),]
      if (! is.null(numTopGenesForPathwayAnalysis)) {
        print(paste0("Using ", numTopGenesForPathwayAnalysis, " top genes for pathway analysis."))
        drim <- unique(drim.padj2$geneName[order(drim.padj2$gene_pAdj)])[1:numTopGenesForPathwayAnalysis]
      } else {
        drim <- unique(drim.padj2$geneName)
        print("Using all significant genes for pathway analysis")
      }
      EGMT <- enricher(gene = drim, TERM2GENE = TERM2GENE, pvalueCutoff = pvalueCutoffPathwayAnalysis)
      eres <- as.data.frame(EGMT)
      if (! length(eres$ID)) {
        warning("NO SIGNIFICANT RESULTS RETURNED FROM PATHWAY ANALYSIS -- adjust p value cutoff and run again.")
      } else {
        dp <- dotplot(EGMT)
        if (! returnDataOnly) {
          ggsave(dp, filename = file.path(outdir, "pathwayAnalysis.png"), height = 5, width = 13)
          fwrite(x = eres, file = file.path(outdir, "pathwayAnalysis.csv"), quote = F, row.names = F, sep = ",")
        }
        resList[["pathwayAnalysis"]] <- EGMT
      }
    }
    print("DONE")
  }
  if (! returnDataOnly) {
    save(resList, file = file.path(outdir, "DTU_results_list.RData"))
  }
  print("PIPELINE FINISHED: ")
  timestamp()
  return(resList)
}


# This has to be in the globalEnv for some reason
.stageWiseTest <- function(pScreen, pConfirmation, alpha,
                           method=c("none","holm","dte","dtu","user"),
                           adjustment=NULL, tx2gene=NULL, pScreenAdjusted,
                           allowNA=FALSE){
  
  if(allowNA){
    if(any(is.na(pScreen))){
      naFeatures <- which(is.na(pScreen))
      message(paste0("Removing ",length(naFeatures),
                     " features with NA screening hypothesis p-values. \n"))
      pScreen <- pScreen[-naFeatures]
      pConfirmation <- pConfirmation[-naFeatures,]
    }
  }
  
  ## check for NA values
  if(!allowNA){
    if(any(is.na(pScreen)) | any(is.na(pConfirmation))){
      stop("NA p-values found in either the screening or confirmation tests.
           If you want to allow for NA p-values, set allowNA=TRUE.")
    }
    }
  method <- match.arg(method,c("none","holm","dte","dtu","user"))
  
  #screening stage
  if(!pScreenAdjusted)
    padjScreen <- p.adjust(pScreen,"BH") else
      padjScreen <- pScreen
  significanceOrdering <- order(padjScreen)
  genesStageI <- padjScreen<=alpha
  
  if(method=="none"){
    
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),
                               ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),
                                             colnames(pConfirmation)))
    pAdjConfirmation[genesStageI,] <- pConfirmation[genesStageI,]
    padjScreenReturn <- padjScreen
    
  } else if(method=="holm"){
    
    padjScreenReturn <- padjScreen
    ## only do correction for genes that passed the screening stage
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),
                               ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),
                                             colnames(pConfirmation)))
    
    for(k in seq_len(sum(genesStageI))){
      row <- pConfirmation[which(genesStageI)[k],]
      # Holm correction conditional on passing the screening stage.
      o <- order(row)
      if(all(!is.na(row))){ #if no NA's, standard Holm with screening correct.
        n <- length(row)
      } else { #if NA's present, only correct for non NA p-values
        n <- length(row[!is.na(row)])
      }
      # Holm adjustment: passing screening stage implies 1 false hypothesis
      adjustment <- c(n-1,(n-1):1)
      if(length(adjustment)!=length(row)){
        adjustment <- c(adjustment, rep(1,length(row)-length(adjustment)))
      }
      rowAdjusted <- row[o]*adjustment
      rowAdjusted <- pmin(rowAdjusted,1)
      rowAdjusted <- cummax(rowAdjusted)
      rowBack <- vector(length=length(row))
      rowBack[o] <- rowAdjusted
      pAdjConfirmation[which(genesStageI)[k],] <- rowBack
    }
    
  } else if(method=="user"){
    if(length(adjustment)!=ncol(pConfirmation)){
      stop("the length of the adjustment vector is not equal to the number of
           confirmation hypotheses as defined by the number of
           columns in pConfirmation.")
    }
    padjScreenReturn=padjScreen
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),
                               ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),
                                             colnames(pConfirmation)))
    for(k in seq_len(sum(genesStageI))){
      row <- pConfirmation[which(genesStageI)[k],]
      o <- order(row)
      rowAdjusted <- row[o]*adjustment
      rowAdjusted <- pmin(rowAdjusted,1)
      # check monotone increase of adjusted p-values
      rowAdjusted <- cummax(rowAdjusted)
      rowBack <- vector(length=length(row))
      rowBack[o] <- rowAdjusted
      rowBack
      pAdjConfirmation[which(genesStageI)[k],] <- rowBack
    }
    
    } else if(method=="dte"){
      
      if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))){
        stop("not all transcript names in pConfirmation match with
             a transcript ID from the tx2gene object.")
      }
      if(any(is.na(match(names(pScreen),tx2gene[,2])))){
        stop("not all gene names in pScreen match with
             a gene ID from the tx2gene object.")
      }
      significantGenes <- names(padjScreen)[genesStageI]
      geneForEachTx <- tx2gene[match(rownames(pConfirmation),tx2gene[,1]),2]
      txLevelAdjustments <- sapply(significantGenes,function(gene){
        id <- which(geneForEachTx %in% gene)
        row <- pConfirmation[id,]
        #make sure names are passed along if only one tx
        if(length(id)==1) names(row)=rownames(pConfirmation)[id]
        o <- order(row)
        n <- length(row)
        # DTE adjustment: passing screening stage implies 1 false hypothesis
        if(n==1) adjustment=0 else adjustment=c(n-1,(n-1):1)
        rowAdjusted <- row[o]*adjustment
        rowAdjusted <- pmin(rowAdjusted,1)
        rowAdjusted <- cummax(rowAdjusted)
        rowBack <- vector(length=length(row))
        rowBack[o] <- rowAdjusted
        names(rowBack) <- names(row)
        rowBack
      }, simplify=FALSE)
      pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=1)
      rownames(pAdjConfirmation) <- paste0(geneForEachTx,":",rownames(pConfirmation))
      # adjusted p-values for screening hypothesis
      padjScreenReturn <- padjScreen[geneForEachTx]
      # adjusted p-values for confirmation hypothesis
      #idCon <- names(unlist(txLevelAdjustments))
      # replace '.' by ':' in names to avoid confusion with ENSEMBL version names
      #idCon <- gsub(x=idCon,pattern=".",replacement=":",fixed=TRUE)
      namesList=names(txLevelAdjustments)
      namesListElements=lapply(txLevelAdjustments,names)
      idCon <- unlist(sapply(seq_len(length(namesList)), function(ii){
        gsub(x=paste(namesList[ii],namesListElements[[ii]]),
             pattern=" ",replace=":")
      }))
      pAdjConfirmation[idCon,1] <- unlist(txLevelAdjustments)
      
      } else if(method=="dtu"){
        
        if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))){
          stop("not all transcript names in pConfirmation match with
               a transcript ID from the tx2gene object.")
        }
        if(any(is.na(match(names(pScreen),tx2gene[,2])))){
          stop("not all gene names in pScreen match with
               a gene ID from the tx2gene object.")
        }
        # adjust screening
        significantGenes <- names(padjScreen)[genesStageI]
        geneForEachTx <- as.character(tx2gene[match(rownames(pConfirmation),
                                                    tx2gene[,1]),2])
        txLevelAdjustments <- sapply(significantGenes,function(gene){
          id <- which(geneForEachTx %in% gene)
          row <- pConfirmation[id,]
          o <- order(row)
          n <- length(row)
          # DTU adjustment: passing screening stage implies 2 false hypotheses
          if(n==1) stop("Some genes have only one transcript; this is incompatible with DTU correction. Remove these transcripts.")
          if(n==2) adjustment=c(0,0) else adjustment=c(n-2,n-2,(n-2):1)
          rowAdjusted <- row[o]*adjustment
          rowAdjusted <- pmin(rowAdjusted,1)
          rowAdjusted <- cummax(rowAdjusted)
          rowBack <- vector(length=length(row))
          rowBack[o] <- rowAdjusted
          names(rowBack) <- names(row)
          rowBack
        }, simplify=FALSE)
        pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=1)
        rownames(pAdjConfirmation) <- paste0(geneForEachTx,":",
                                             rownames(pConfirmation))
        # adjusted p-values for screening hypothesis
        padjScreenReturn <- padjScreen[as.character(geneForEachTx)]
        # adjusted p-values for confirmation hypothesis
        #idCon <- names(unlist(txLevelAdjustments))
        # replace '.' by ':' in names to avoid confusion with ENSEMBL version names
        #idCon <- gsub(x=idCon,pattern=".",replacement=":",fixed=TRUE)
        namesList=names(txLevelAdjustments)
        namesListElements=lapply(txLevelAdjustments,names)
        idCon <- unlist(sapply(seq_len(length(namesList)), function(ii){
          gsub(x=paste(namesList[ii],namesListElements[[ii]]),
               pattern=" ",replace=":")
        }))
        pAdjConfirmation[idCon,1] <- unlist(txLevelAdjustments)
        
        } else stop("method must be either one of 'holm' or ... ")
  
  #BH-adjusted s.l.
  alphaAdjusted <- sum(padjScreen<=alpha)/length(padjScreen)*alpha
  #Correct FWER-adjusted p-values acc. to BH-adjusted s.l.
  naPAdj <- is.na(pAdjConfirmation)
  pAdjBH <- pAdjConfirmation[!naPAdj]*length(padjScreen)/sum(padjScreen<=alpha)
  pAdjConfirmation[!naPAdj] <- pmin(pAdjBH,1)
  if(!(method %in% c("dte","dtu"))){
    pAdjStage <- cbind(padjScreenReturn,pAdjConfirmation)
    colnames(pAdjStage)[1] <- "padjScreen"
  }
  if(method %in% c("dte","dtu")){
    pAdjStage <- cbind(pAdjConfirmation,padjScreenReturn)[,2:1]
    colnames(pAdjStage) <- c("gene","transcript")
  }
  return(list(pAdjStage=pAdjStage, alphaAdjusted=alphaAdjusted))
}



DS_pipeline <- function(
  colData,
  bamFolder,
  outdir,
  contrast,
  gtf,
  toRun = c("plotFeatures", "pathwayAnalysis"),
  TERM2GENE = NULL,
  ens2symb = NULL,
  ens2symbRaw = NULL,
  cores = 12,
  txf = NULL,
  sgfc = NULL,
  returnDataOnly = FALSE,
  numPlotFeatures = 20,
  AllEvents_RNASeq = NULL,
  EventsFinal = NULL
) {
  # START PIPELINE
  require(SGSeq)
  require(biomaRt)
  require(DOSE)
  require(ggpubr)
  require(ChIPpeakAnno)
  require(data.table)
  require(rtracklayer)
  require(GenomicFeatures)
  require(gplots)
  require(pheatmap)
  require(clusterProfiler)
  require(limma)
  require(EventPointer)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  # # Bug test
  # colData = colData
  # sgfc = sgfc
  # bamFolder = bamFolder
  # AllEvents_RNASeq = AllEvents_RNASeq
  # outdir = "Results/DS/EWS_vs_Healthy/EWS_vs_NCSC"
  # contrast = c("category", "EWS", "NCSC")
  # gtf = gtf
  # TERM2GENE = TERM2GENE
  # ens2symb = NULL
  # returnDataOnly = T
  # EventsFinal = NULL
  
  # Initialize results object
  resList <- list()
  
  # Make outdir
  if (! dir.exists(outdir) & ! returnDataOnly) {
    dir.create(outdir)
  }
  
  # Get ens2symb object if does not exist
  if (is.null(ens2symb)) {
    print("No ens2symb object -- gathering from biomaRt...")
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens2symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
    if (is.null(ens2symb)) {
      ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
      ens2symb <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
      if (is.null(ens2symb)) {
        stop("Cannot contact biomaRt at this time... please supply ens2symb manually or try again later")
      }
    } 
    print("DONE")
  } 
  
  
  # Wrangle colData
  contrastCol <- colnames(colData)[which(colnames(colData) == contrast[1])]
  # Wrangle colData
  colData <- colData[which(colData[, contrastCol] %in% c(contrast[2], contrast[3])),]
  colData[, contrastCol] <- factor(colData[, contrastCol], levels = c(contrast[2], contrast[3]))
  
  if (is.null(sgfc)) {
    folders_bam <- file.path(bamFolder, rownames(colData))
    names(folders_bam) <- rownames(colData)
    
    # Find bam files -- verify existence
    if (! all(file.exists(file.path(folders_bam)))) {
      warning("Not all samples have a corresponding bam file folder... continuing anyways.")
    } 
    
    
    file_all <- list.files(file.path(folders_bam), full.names = T)
    
    # Check for bam file
    file_bam <- file_all[grep(file_all, pattern = "^.+bam$")]
    
    if (! length(file_bam) == length(folders_bam)) {
      stop("Not all bam folders contain a bam file or some contain multiple files...")
    }
    
    # Check for ResultLog file
    file_result <- file_all[grep(file_all, pattern = "ResultLog.out")]
    
    if (! length(file_result) == length(folders_bam)) {
      resdirs <- dirname(file_result)
      bamdirs <- dirname(file_bam)
      msg <- paste0("\nNot all bam folders contain a result log file. Missing: ", bamdirs[which(! bamdirs %in% resdirs)])
      stop(msg)
    }
    
    # Check for final result file
    file_result2 <- file_all[grep(file_all, pattern = "ResultLog.final.out")]
    
    if (! length(file_result2) == length(folders_bam)) {
      resdirs <- dirname(file_result2)
      bamdirs <- dirname(file_bam)
      msg <- paste0("\nNot all bam folders contain a final results file. Missing: ", bamdirs[which(! bamdirs %in% resdirs)])
      stop(msg)
    }
    
    print("All necessary files exist -- continuing with differential splicing analysis...")
    
    print("Getting bam file info! ... ")
    
    si <- determineLibraryInfoFromSTAR(sampleNames = rownames(colData), resultLogs = file_result, 
                                       resultLogFinals = file_result2, bam_files = file_bam)
    resList[["si"]] <- si
    print("DONE")
    
    if (is.null(txf)) {
      print("Loading TXF data from gtf file")
      txdb <- makeTxDbFromGFF(file = gtf)
      txf <- convertToTxFeatures(txdb)
      txf <- keepStandardChromosomes(x = txf, pruning.mode = "coarse")
      if (! returnDataOnly) {
        save(txf, file = file.path(outdir, "txf.RData"))
      }
    } else {
      # Just to be safe...
      txf <- keepStandardChromosomes(x = txf, pruning.mode = "coarse")
      if (! returnDataOnly) {
        save(txf, file = file.path(outdir, "txf.RData"))
      }
    }
    
    
    
    colnames(ens2symb) <- c("geneName", "geneSymb")
    
    resList[["txf"]] <- txf
    
    print("Starting SGSeq run...")
    timestamp()
    
    sgfc <- analyzeFeatures(si, features = txf, cores = cores)
    
    if (! returnDataOnly) {
      save(sgfc, file = file.path(outdir, "sgfc.RData"))
    }
    resList[["sgfc"]] <- sgfc
    print("SGSeq DONE")
    timestamp()
    
  }
  
  print("Starting event pointer for differential analysis of feature counts...")
  if (is.null(EventsFinal)) {
    if (returnDataOnly) {
      TxtPath <- tempdir()
    } else {
      TxtPath <- file.path(outdir, "eventPointerOut")
      dir.create(TxtPath)
    }
    if (is.null(AllEvents_RNASeq)) {
      timestamp()
      print("Detecting events ... ")
      AllEvents_RNASeq <- EventDetection(sgfc, cores=1, Path=TxtPath)
    }
    if (! returnDataOnly) {
      save(AllEvents_RNASeq, file = file.path(outdir, "AllEvents_RNASeq.RData"))
    }
    print("DONE")
    timestamp()
    print("Getting final events info")
    #load(file.path(outdir, "AllEvents_RNASeq.RData"))
    group <- colData[,contrastCol]
    Dmatrix<-model.matrix(~0 + group)
    colnames( Dmatrix ) <- levels( group )
    contrastStrin <- paste0(contrast[2], "-", contrast[3])
    contrasts <- contrastStrin
    levels <- c(contrast[2], contrast[3])
    # mod of makeContrasts... 
    # The package version didn't work in a subshell for some reason
    n <- length(levels)
    if (n < 1)
      stop("No levels to construct contrasts from.")
    if (is.factor(levels))
      levels <- levels(levels)
    if (!is.character(levels))
      levels <- colnames(levels)
    
    indicator <- function(i,n) {
      out <- rep(0,n)
      out[i] <- 1
      out
    }
    #In order to remove invalid level values which are not syntactically valid variable names in R and, for example, do not begin with a letter.
    from <- levels
    to <- make.names(levels)
    
    gsub2 <- function(pattern, replacement, x, ...) {
      for(i in 1:length(pattern))
        x <- gsub(pattern[i], replacement[i], x, ...)
      x
    }
    
    e <- gsub2(from, to, contrasts)
    levelsenv <- new.env()
    for (i in 1:n) assign(to[i], indicator(i, n), pos = levelsenv)
    ne <- length(contrasts)
    L <- matrix(0, nrow=length(levels), ncol=length(contrasts),
                dimnames=list(Levels = levels, Contrasts = contrasts))
    for (j in 1:ne) {
      ej <- parse(text = e[j])
      L[, j] <- eval(ej, envir = levelsenv)
    }
    Cmatrix <- L
    # Final input for EventPointer
    EventsFinal <- EventPointer_RNASeq(Events = AllEvents_RNASeq, 
                                       Design = Dmatrix, 
                                       Contrast = Cmatrix, 
                                       Statistic = "LogFC", 
                                       PSI = TRUE)
    if (! returnDataOnly) {
      save(EventsFinal, file = file.path(outdir, "EventsFinal.RData"))
    }
    resList[["EventsFinal"]] <- EventsFinal
    print("DONE")
  } else {
    print("Loading events info...")
    if (! returnDataOnly) {
      save(EventsFinal, file = file.path(outdir, "EventsFinal.RData"))
    }
    
    print("DONE")
  }
  
  # Wrangling EventsFinals
  print("Analyzing events data ...")
  EventsFinalRaw <- EventsFinal
  # EventsFinal <- EventsFinalRaw
  EventsFinal$Gene <- substr(EventsFinal$Gene, 1, 15)
  EventsFinal <- merge(x = ens2symb, y = EventsFinal, by.x = "ensembl_gene_id", by.y = "Gene")
  colnames(EventsFinal)[c(1,2)] <- c("geneID", "geneName")
  EventsFinal$GSEA <- -log10(EventsFinal$Pvalue) * sign(EventsFinal$`Delta PSI`)
  EventsFinal <- EventsFinal[order(EventsFinal$GSEA, decreasing = T),]
  resList[["EventsFinal"]] <- EventsFinal
  if (! returnDataOnly) {
    fwrite(EventsFinal, file = file.path(outdir, "EventsFinal.csv"), sep = ",", row.names = F, quote = F)
  }
  titleStr <- paste0(contrast[2], " vs. ", contrast[3])
  # Get event type plots
  EventsFinal_sig_up <- EventsFinal[which(EventsFinal$GSEA > 1.3),]
  EventsFinal_sig_dn <- EventsFinal[which(EventsFinal$GSEA < -1.3),]
  df <- as.data.frame(table(EventsFinal_sig_up$Event_Type))
  df$Freq <- (df$Freq/sum(df$Freq))*100
  df$Group <- "Overexpressed"
  df2 <- as.data.frame(table(EventsFinal_sig_dn$Event_Type))
  df2$Freq <- (df2$Freq/sum(df2$Freq))*100
  df2$Group <- "Underexpressed"
  df3 <- rbind(df, df2)
  capText <- paste0("# Increasing events: ", length(EventsFinal_sig_up$geneID),
                    "\n# Decreasing events: ", length(EventsFinal_sig_dn$geneID))
  plt <- ggbarplot(data = df3, x = "Var1", y = "Freq", fill = "Group", caption = capText,
                   color = "Group",  palette = "Paired", position = position_dodge(0.9), 
                   xlab = "Event Type", ylab = "Percentage of significant splice events\n", 
                   title = paste0(titleStr,  " Differential Splicing Events\n")) + rotate_x_text()
  resList[["EventType_plot"]] <- plt
  
  if (! returnDataOnly) {
    ggsave(plt, filename = file.path(outdir, "DS_events.png"), height = 8, width = 6)
  }
  
  print("DONE")
  
  # IGV visualization
  print("Returning IGV data...")
  bed <- EventsFinal
  trackstr <- strsplit(EventsFinal$Position, ":")
  trackstr2 <- sapply(trackstr, "[[", 2)
  trackstr <- sapply(trackstr, "[[", 1)
  trackstr3 <- strsplit(trackstr2, "-")
  start <- sapply(trackstr3, "[[", 1)
  end <- sapply(trackstr3, "[[", 2)
  badInd <- which(start > end)
  length(start) == length(end) & length(end) == length(trackstr) & length(trackstr) == length(bed$geneID)
  start <- start[-badInd]
  end <- end[-badInd]
  trackstr <- trackstr[-badInd]
  bed <- bed[-badInd, ]
  length(start) == length(end) & length(end) == length(trackstr) & length(trackstr) == length(bed$geneID)
  bed$seqnames <- trackstr
  bed$start <- start
  bed$end <- end
  bed <- bed[,c(-1)]
  bed <- unique(bed)
  bed <- bed[which(bed$Pvalue < .05),]
  bed <- bed[which(abs(bed$`Delta PSI`) > .1),]
  bed <- bed[,c(8,9,10,1:7)]
  colnames(bed)[10] <- "score"
  bed <- bed[,c(-6)]
  bed <- toGRanges(bed)
  names(bed) <- paste0(bed$geneName, "_", gsub(x = bed$Event_Type, pattern = " ", replacement = "_"),
                       "_", signif(bed$score, 3))
  if (! returnDataOnly) {
    export(bed, format = "bed", con = file.path(outdir, paste0(contrast[2], "_vs_", contrast[3], "_DS.bed")))
    export(bed, format = "bedGraph", con = file.path(outdir, paste0(contrast[2], "_vs_", contrast[3], "_DS.bedGraph")))
  }
  resList[["dsGR"]] <- bed
  if ("plotFeatures" %in% toRun) {
    print("Plotting splice features ... ")
    names(bed) <- NULL
    bed <- as.data.frame(bed)
    bed <- bed[,c(6:11)]
    bed <- unique(bed)
    bed <- bed[order(bed$score),]
    genes <- bed$geneName
    genes <- unique(genes)
    n <- length(genes)
    m <- floor(numPlotFeatures/2)
    genesToPlot <- genes[c(1:m, (n-m+1):n)]
    map <- ens2symb[which(ens2symb$external_gene_name %in% genesToPlot),]
    p <- length(map$ensembl_gene_id)
    group <- colData[,contrastCol, drop = F]
    cd <- as.data.frame(colData(sgfc))
    group$sample_name <- rownames(group)
    group2 <- merge(x = group, y = cd, by = "sample_name")
    group2 <- group2[,c(1, 2, 4)]
    rownames(group2) <- group2$sample_name
    group2 <- group2[,c(-1)]
    annoCol <- group2[order(match(rownames(group2), rownames(colData(sgfc)))),, drop = F]
    annoCol$paired_end <- as.character(annoCol$paired_end)
    if (! returnDataOnly) {
      dir.create(file.path(outdir, "featurePlots"))
      
    }
    featurePlots <- list()
    features <- as.data.frame(rowRanges(sgfc))
    features$geneName <- substr(as.character(features$geneName), 1, 15)
    for (i in 1:p) {
      cat("\n", i, " of ", p)
      id <- map$ensembl_gene_id[i]
      geneName <- map$external_gene_name[i]
      

      features2 <- features[which(features$geneName == id & features$type %in% c("E", "J")),]
      if (! length(features2$seqnames)) {
        next
      }
      annoRow <- data.frame(row.names = rownames(features2), "featureType" = features2$type)
      ind <- features2$featureID
      cts <- log2(FPKM(sgfc[ind,]) +1)
      rownames(cts) <- ind
      
      ph <- pheatmap(cts, show_rownames = F, show_colnames = F, cluster_rows = F, 
                     main = paste0(geneName, " diffSplice plot\n"), 
                     scale = "row", annotation_row = annoRow, annotation_col = annoCol,
                     color = greenred(100))
      if (! returnDataOnly) {
        ggsave(ph, filename = file.path(outdir, "featurePlots", paste0(geneName, ".png")))
      }
      featurePlots[[i]] <- ph
      names(featurePlots)[i] <- geneName
    }
    resList[["featurePlots"]] <- featurePlots
    print("DONE")
  }
  if ("pathwayAnalysis" %in% toRun) {
    print("Pathway analysis ... ")
    print("Analysis of event direction ... ")
    EventsFinal_sig_up <- EventsFinal_sig_up[,"geneName", drop = F]
    EventsFinal_sig_dn <- EventsFinal_sig_dn[,"geneName", drop = F]
    EventsFinal_sig_up$Group <- "Upregulated"
    EventsFinal_sig_dn$Group <- "Downregulated"
    EventsFinal2 <- rbind(EventsFinal_sig_up, EventsFinal_sig_dn)
    pathUp <- enricher(EventsFinal_sig_up$geneName,  pvalueCutoff = 0.05,
                       TERM2GENE = hsapiens_complex_TERM2GENE)
    pathUperes <- as.data.frame(pathUp)
    pathUperes$group <- "Upregulated"
    pathDn <- enricher(EventsFinal_sig_dn$geneName,  pvalueCutoff = 0.05,
                       TERM2GENE = hsapiens_complex_TERM2GENE)
    pathDneres <- as.data.frame(pathDn)
    pathDneres$group <- "Downregulated"
    
    
    path <- enricher(unique(c(EventsFinal_sig_up$geneName, EventsFinal_sig_dn$geneName)), pvalueCutoff = 0.05,
                     TERM2GENE = hsapiens_complex_TERM2GENE)
    patheresraw <- as.data.frame(path)
    patheresraw$group <- "Either"
    patheres2 <- rbind(pathUperes, pathDneres, patheresraw)
    patheres <- rbind(pathUperes[c(1:20),], pathDneres[c(1:20),], patheresraw[c(1:20),])
    patheres <- patheres[order(patheres$Count, decreasing = T),]
    idx <- order(patheres$Count, decreasing = T)
    patheres$Description <- factor(patheres$Description, levels=rev(unique(patheres$Description[idx])))
    patheres$group <- factor(patheres$group, levels = c("Upregulated", "Downregulated", "Either"))
    gp <- ggplot(patheres, aes_string(x="Count", y="Description", color="p.adjust")) +
      geom_point() +
      scale_color_continuous(low="red", high="blue", name = "p.adjust", 
                             guide=guide_colorbar(reverse=TRUE, title = "FDR")) + 
      ylab(NULL) + scale_size(range=c(3, 8)) + facet_grid(~group) + theme_dose(font.size = 10.5)
    
    resList[["enrichmentPlot_up_dn_either"]] <- gp
    resList[["enrichmentPlot_up_dn_either"]] <- gp
    
    if (! returnDataOnly) {
      ggsave(gp, filename = file.path(outdir, "pathwayEnrichment_by_event_direction.png"), height = 6, width = 12)
    }
    
    EventsFinal_sig <- EventsFinal[which(EventsFinal$Pvalue < .05),]
    events <- as.data.frame(table(EventsFinal_sig$Event_Type))
    events <- events[which(events$Freq > 100),]
    eventTypes <- as.character(events$Var1)
    print("Pathway analysis per event type ...")
    
    for (i in 1:length(eventTypes)) {
      eventType <- eventTypes[i]
      print(eventType)
      typeFinal <- EventsFinal_sig[which(EventsFinal_sig$Event_Type == eventType),]
      geneList <- unique(typeFinal$geneName)
      EGMT <- enricher(gene = geneList, TERM2GENE = hsapiens_complex_TERM2GENE)
      eres <- as.data.frame(EGMT)
      eres$group <- eventType
      eres2 <- eres[c(1:10),]
      if (i == 1) {
        typeFrame <- eres2
        typeFrame2 <- eres
      } else {
        typeFrame <- rbind(typeFrame, eres2)
        typeFrame2 <- rbind(typeFrame2, eres)
      }
    }
    
    typeFrame <- typeFrame[order(typeFrame$Count, decreasing = T),]
    idx <- order(typeFrame$Count, decreasing = T)
    typeFrame$Description <- factor(typeFrame$Description, levels=rev(unique(typeFrame$Description[idx])))
    gp2 <- ggplot(typeFrame, aes_string(x="Count", y="Description", color="p.adjust")) +
      geom_point() +
      scale_color_continuous(low="red", high="blue", name = "p.adjust", 
                             guide=guide_colorbar(reverse=TRUE, title = "FDR")) + 
      ylab(NULL) + facet_grid(~group) + theme_dose(font.size = 8.5)
    
    resList[["enrichmentPlot_splice_types"]] <- gp2
    typeFrame2 <- rbind(patheres2, typeFrame2)
    resList[["enrichmentData_splice_types"]] <- typeFrame2
    if (! returnDataOnly) {
      ggsave(gp2, filename = file.path(outdir, "pathwayEnrichment_by_event_type.png"), height = 10, width = 20)
      fwrite(typeFrame2, file = file.path(outdir, "pathwayEnrichment.csv"), row.names = F, quote = F)
    }
    print("DONE")
  }
  print("PIPELINE DONE ... RETURNING RESULTS ...")
  timestamp()
  return(resList)
}


determineLibraryInfoFromSTAR <- function(sampleNames, resultLogs, resultLogFinals, bam_files) {
  # # Bug testing
  # sampleNames <- samples$sample_id
  # resultLogs <- file_result
  # resultLogFinals <- file_result2
  # bam_files <- file_bam
  
  n <- length(sampleNames)
  paired_end <- rep(F, n)
  read_length <- rep(NA, n)
  frag_length <- rep(NA, n)
  lib_size <- rep(NA, n)
  
  for (i in 1:n) {
    cat("\n", i, " of ", n)
    # Read log files and get data
    lines <- readLines(resultLogs[i])
    line <- lines[grep(x = lines, pattern = "--readFilesIn.+")][1] #CMD line
    line <- strsplit(line, split = " ")
    line <- unlist(line, use.names = F)
    res <- grep(x = line, pattern = "fastq")
    if (length(res) == 2) {
      PE <- T
      paired_end[i] <- PE
    } else {
      PE <- F
    }
    # From final log -- read and frag length
    lines <- readLines(resultLogFinals[i])
    line <- lines[grep(x = lines, pattern = "Average input read length")][1] #CMD line
    line <- strsplit(line, split = "\t")
    line <- sapply(line, "[[", 2)
    lent <- as.numeric(line)
    if (PE) {
      lent <- floor(lent/2)
      # Get median frag length
      fl <- system(command = paste0("samtools view ", bam_files[i], " 2>/dev/null | head -n 100000 | awk '{print $9}' "),
                   intern = T, ignore.stderr = T, ignore.stdout = F)
      fl <- as.numeric(fl)
      fl <- median(abs(fl), na.rm = T)
      frag_length[i] <- fl
    } else {
      frag_length[i] <- NA_real_
    }
    read_length[i] <- lent
    # From final log -- library size 
    lines <- readLines(resultLogFinals[i])
    line <- lines[grep(x = lines, pattern = "Uniquely mapped reads number")][1] #CMD line
    line <- strsplit(line, split = "\t")
    line <- sapply(line, "[[", 2)
    libSize <- as.numeric(line)
    lib_size[i] <- libSize
  }
  
  si <- data.frame(sample_name = sampleNames, file_bam = bam_files,
                   paired_end = paired_end, read_length = read_length,
                   frag_length = frag_length, lib_size = lib_size, stringsAsFactors = F)
  
  return(si)
  
}









