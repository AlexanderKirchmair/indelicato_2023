
## R analysis functions

nf_import <- function (nfdir, ...){
  
  cat(crayon::red("Use '$counts' for 3'-tagged sequencing data\n"))
  files <- list.files(path = nfdir, pattern = "quant.sf", recursive = TRUE, 
                      full.names = TRUE)
  names(files) <- sapply(strsplit(files, "/", fixed = TRUE), 
                         function(x) x[[length(x) - 1]])
  tx2gene <- read.delim(file.path(nfdir, "star_salmon", "salmon_tx2gene.tsv"), 
                        header = FALSE)
  colnames(tx2gene) = c("tx", "gene_id", "gene_name")
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene[, 
                                                                      c(1, 3)], ...)
  txi$counts <- matrix(as.integer(round(txi$counts)), nrow = nrow(txi$counts), 
                       dimnames = dimnames(txi$counts))
  txi
}



nf_summary <- function (nfdir, design = NULL, ignore = FALSE){
  
  multiqc_dir <- file.path(nfdir, "multiqc", "star_salmon", "multiqc_data")
  
  multiqc_files <- list.files(multiqc_dir, pattern = "multiqc_", 
                              full.names = TRUE)
  multiqc_files <- multiqc_files[grep("txt", multiqc_files)]
  names(multiqc_files) <- gsub("multiqc_|\\.txt", "", basename(multiqc_files))
  multiqc_data <- lapply(multiqc_files, function(tmpfile) read.delim(file = tmpfile))
  multiqc <- data.frame(row.names = multiqc_data$general_stats$Sample, 
                        sample = multiqc_data$general_stats$Sample)
  fastqc_raw <- multiqc_data$fastqc[match(multiqc$sample, multiqc_data$fastqc$Sample), 
  ]
  multiqc$raw <- fastqc_raw$Total.Sequences
  
  if (!is.null(multiqc_data$fastqc_1)){
    fastqc_trimmed <- multiqc_data$fastqc_1[match(multiqc$sample, 
                                                multiqc_data$fastqc_1$Sample), ]
    multiqc$trimmed <- multiqc$raw - fastqc_trimmed$Total.Sequences
  }
  
  star <- multiqc_data$star[match(multiqc$sample, multiqc_data$star$Sample), ]
  
  if (!ignore) 
    stopifnot(all(star$total_reads == multiqc$raw - multiqc$trimmed))
  multiqc$unmapped <- rowSums(star[, grepl("unmapped_", colnames(star)) & 
                                     !grepl("percent", colnames(star))])
  multiqc$multimapped <- star$multimapped + star$multimapped_toomany
  multiqc$mapped <- star$uniquely_mapped
  
  if (!ignore) 
    stopifnot(all(star$total_reads == multiqc$unmapped + 
                    multiqc$multimapped + multiqc$mapped))
  rseqc <- multiqc_data$rseqc_read_distribution[match(multiqc$sample, 
                                                      multiqc_data$rseqc_read_distribution$Sample), ]
  exon_fraction <- (rseqc$cds_exons_tag_count + rseqc$X5_utr_exons_tag_count + 
                      rseqc$X3_utr_exons_tag_count)/rseqc$total_tags
  multiqc$nonexon <- star$uniquely_mapped * (1 - exon_fraction)
  multiqc$exon <- star$uniquely_mapped * exon_fraction
  
  if (!ignore & !is.null(multiqc$trimmed)) {
    with(multiqc, stopifnot(all(nonexon + exon + multimapped + unmapped + trimmed == raw)))
  }
  
  if (!ignore & is.null(multiqc$trimmed)) {
    with(multiqc, stopifnot(all(nonexon + exon + multimapped + unmapped == raw)))
  }
  
  ggdf <- multiqc %>% dplyr::select(-c(raw, mapped)) %>% tidyr::pivot_longer(cols = -sample, 
                                                                             names_to = "fate", values_to = "reads") %>% as.data.frame()
  ggdf$fate <- factor(ggdf$fate, ordered = TRUE, levels = c("unmapped", 
                                                            "multimapped", "trimmed", "nonexon", "exon"))
  if (!is.null(design)) {
    ggdf <- subset(ggdf, sample %in% rownames(design))
    ggdf$sample <- factor(ggdf$sample, ordered = TRUE, levels = rownames(design))
  }
  gg <- ggplot2::ggplot(data = ggdf, ggplot2::aes(x = sample, 
                                                  y = reads, fill = fate)) + ggplot2::theme_bw(base_size = 10) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   axis.line = ggplot2::element_line(colour = "black"), 
                   axis.text = ggplot2::element_text(colour = "black"), 
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, 
                                                       vjust = 0.5)) + ggplot2::geom_bar(stat = "identity", 
                                                                                         width = 0.65) + ggplot2::scale_fill_manual(values = c("#e34c39", 
                                                                                                                                               "#d99b27", "#e3e649", "#1ea9d4", "#1136ed")) + ggplot2::labs(x = "", 
                                                                                                                                                                                                            y = "million reads", fill = "") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
                                                                                                                                                                                                                                                                                                               0.02)))
  list(data = multiqc, plot = gg)
}









runDESeq2 <- function(data, design = NULL, formula = ~ 1, contrasts = NULL, ctrlgenes = NULL, sizefactors = NULL,
                      alpha = 0.05, ordered = TRUE, df = TRUE, ncores = NULL,
                      shrink = TRUE, ihw = TRUE, vst = TRUE, rlog = FALSE, ...){
  
  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  if (ifelse(is.null(ncores), TRUE, ncores > 1)) stopifnot(requireNamespace("BiocParallel", quietly = TRUE))
  if (ihw == TRUE) stopifnot(requireNamespace("IHW", quietly = TRUE))
  if (shrink == TRUE) stopifnot(requireNamespace("ashr", quietly = TRUE))
  
  results <- list()
  
  
  # Parallel setup ----
  if (is.null(ncores)) ncores <- min(c(4, max(c(1, length(contrasts)))))
  bppar <- NULL
  if (ncores > 1){
    BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores))
    bppar <- BiocParallel::bpparam()
    message(paste0("Parallel on ", ncores, " cores..."))
  }
  
  
  # Input data ----
  if (class(data) %in% "SummarizedExperiment"){
    if (is.null(design)) design <- data.frame(row.names = colnames(assays(data)[[1]]))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSet(data, colData = design, design = formula)
    
  } else if (class(data) %in% "matrix"){
    if (is.null(design)) design <- data.frame(row.names = colnames(data))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromMatrix(data, colData = design, design = formula)
    
  } else if (class(data) %in% "list") {
    if (is.null(design)) design <- data.frame(row.names = colnames(data$abundance))
    for (i in seq_along(data)){
      if (!is.null(ncol(data[[i]]))){
        data[[i]] <- data[[i]][,rownames(design)]
      }
    }
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromTximport(data, colData = design, design = formula)
    
  } else {
    stop("Wrong input data format!")
  }
  
  
  # Pre-filtering ----
  dds <- dds[rowSums(counts(dds) > 0, na.rm = TRUE) != 0,]
  
  ### filtering -------------
  
  
  # Model fitting ----
  if (!is.null(sizefactors)){
    dds$sizeFactor <- sizefactors
  } else if (!is.null(ctrlgenes)){
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = rownames(dds) %in% ctrlgenes)
  }
  dds <- DESeq2::DESeq(dds, parallel = (ncores > 1), BPPARAM = bppar, ...)
  
  
  # Normalized counts ----
  results$normcounts <- DESeq2::counts(dds, normalized = TRUE)
  if (vst == TRUE) results$vst <- assays(vst(dds, blind = TRUE))[[1]]
  if (rlog == TRUE) results$rlog <- assays(rlog(dds, blind = TRUE))[[1]]
  
  
  # Contrasts ----
  if (!is.null(contrasts)){
    
    if (ihw == TRUE) filterFun <- IHW::ihw else filterFun <- NULL
    
    results$results <- lapply(contrasts, function(tmp){
      
      if (is.null(filterFun)){
        mle <- DESeq2::results(dds, contrast = tmp, alpha = alpha)
      } else {
        mle <- DESeq2::results(dds, filterFun = filterFun, contrast = tmp, alpha = alpha)
      }
      
      
      if (shrink == TRUE){
        mmse <- DESeq2::lfcShrink(dds, contrast = tmp, res = mle, type = "ashr", svalue = TRUE)
        mmse <- mmse[rownames(mle),]
        mle$svalue <- mmse$svalue
        mle$lfcShrink <- mmse$log2FoldChange
      }
      
      mle$gene <- rownames(mle)
      mle <- as.data.frame(mle)
      mle <- dplyr::relocate(mle, gene)
      if (ordered == TRUE) mle <- dplyr::arrange(mle, padj)
      mle
    })
    
  }
  
  
  # Results ----
  results$dds <- dds
  results$design <- as.data.frame(colData(dds))
  results
  
}









runLIMMA <- function(data, design, formula = ~ 1, contrasts = NULL, trend = TRUE, robust = FALSE, p.adj.method = "fdr", ...){
  
  stopifnot(requireNamespace("limma", quietly = TRUE))
  
  formula <- update(formula,  ~  0 + .)
  
  mm <- model.matrix(formula, design)
  
  contrasts_named <- lapply(contrasts, function(contr){  paste(c(paste0(contr[1], contr[2]), paste0(contr[1], contr[3])), collapse = " - ")})
  contrasts_limma <- do.call(limma::makeContrasts, c(contrasts_named, list(levels = mm)))
  
  fit <- limma::lmFit(data, mm)
  cfit <- limma::contrasts.fit(fit = fit, contrasts = contrasts_limma)
  efit <- limma::eBayes(fit = cfit, trend = trend, robust = robust, ...)
  
  coefnames <- colnames(efit$coefficients)
  results <- lapply(setNames(coefnames, coefnames), function(tmpcoef) limma::topTable(efit, coef = tmpcoef, number = Inf, adjust.method = p.adj.method) )
  results <- lapply(results, function(tmpres) data.frame(row.names = rownames(tmpres),
                                                         "id" = rownames(tmpres),
                                                         "mean" = tmpres$AveExpr,
                                                         "stat" = tmpres$t,
                                                         "lfc" = tmpres$logFC,
                                                         "pvalue" = tmpres$P.Value,
                                                         "padj" = tmpres$adj.P.Val) )
  
  
  results
}


runORA_GO <- function(genes, universe = NULL){
  if (is.null(universe)){
    ora <- clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "fdr", keyType = "SYMBOL")
  } else {
    ora <- clusterProfiler::enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "fdr", keyType = "SYMBOL")
  }
  
  as.data.frame(ora)
}





getGeneSets <- function(collections = c("H", "C2|CP:KEGG|CP:REACTOME"), species = "human", id_type = "gene_symbol", format = c("list", "dataframe", "GeneSetCollection")){
  
  stopifnot(requireNamespace("msigdbr"))
  # msigdbr::msigdbr_species()
  # msigdbr::msigdbr_collections()
  
  if (!is.null(collections)){
    
    collections <- strsplit(collections, split = "|", fixed = TRUE)
    categories <- sapply(collections, function(x) x[1] )
    subcategories <- lapply(setNames(collections, categories), function(x) unique(x[-1]) )
    subcategories[sapply(subcategories, length) == 0] <- NA
    df <- stack(subcategories)[,c(2,1)]
    colnames(df) <- c("category", "subcategory")
    
    
    gslist <- lapply(1:nrow(df), function(i){
      categ <- as.character(df[i,1])
      subcateg <- df[i,2]
      if (is.na(subcateg)) subcateg <- NULL
      msigdbr::msigdbr(species = species, category = categ, subcategory = subcateg)
    })
    
    genesets <- Reduce(x = gslist, f = rbind)
    
  } else {
    
    genesets <- msigdbr::msigdbr(species = species)
  }

  
  
  # format
  format <- tolower(format[1])
  
  if (format == "list"){
    GS <- unstack(unique(data.frame(genesets[,c(id_type, "gs_name")])))
    
  } else if (format == "GeneSetCollection"){
    GS <- unstack(unique(data.frame(genesets[,c(id_type, "gs_name")])))
    GS <- lapply(setNames(seq_along(genesets), names(genesets)), function(i){
      GSEABase::GeneSet(genesets[[i]], setName = names(genesets[i]))
    })
    GS <- GSEABase::GeneSetCollection(GS)
    
  } else if (format %in% c("dataframe", "data.frame", "df", "table")){
    GS <- unique(as.data.frame(genesets)[,c("gs_name", id_type)])
    colnames(GS) <- c("term", "gene")
    
  } else {
    GS <- genesets
    
  }
  
  
  GS
}


saveplot <- function (p, file = NULL, dev = "png", width = 3000, height = 2500, dpi = 300, units = "px", bg = "white", ggsave = TRUE, PRINTFUN = print, ...){
  
  stopifnot(dev %in% c("png", "pdf", "tiff", "svg", "jpeg", 
                       "jpg", "cairo_ps", "ps"))
  
  if (is.null(file)) 
    file <- rlang::as_name(rlang::enquo(p))
  if (length(file) > 1) {
    file <- do.call(file.path, as.list(file))
  }
  if (nat(baseext(file) != dev)) 
    file <- paste0(file, ".", dev)
  get_inches <- function(x, dpi, units) {
    if (units == "px") {
      x <- x/dpi
    }
    if (units == "cm") {
      x <- x/dpi/2.54
    }
    x
  }
  if ("gg" %in% class(p) & ggsave == TRUE & !(dev %in% c("pdf", 
                                                         "ps", "cairo_ps"))) {
    ggplot2::ggsave(filename = file, plot = p, device = dev, 
                    width = width, height = height, units = units, dpi = dpi, 
                    bg = bg, ...)
  }
  else {
    if (dev == "png") {
      if (is.null(bg)) 
        bg <- "white"
      grDevices::png(filename = file, width = width, height = height, 
                     units = units, bg = bg, res = dpi, type = "cairo", 
                     ...)
      PRINTFUN(p)
      dev.off()
    }
    if (dev == "pdf") {
      if (!"list" %in% class(p)) 
        p <- list(p)
      grDevices::pdf(file = file, width = get_inches(width, 
                                                     dpi = dpi, units = units), height = get_inches(height, 
                                                                                                    dpi = dpi, units = units), bg = bg, onefile = TRUE, 
                     ...)
      tmp <- lapply(p, PRINTFUN)
      dev.off()
    }
    if (dev == "svg") {
      grDevices::svg(filename = file, width = get_inches(width, 
                                                         dpi = dpi, units = units), height = get_inches(height, 
                                                                                                        dpi = dpi, units = units), onefile = TRUE, ...)
      PRINTFUN(p)
      dev.off()
    }
    if (dev == "tiff") {
      grDevices::tiff(filename = file, width = width, height = height, 
                      res = dpi, bg = bg, type = "cairo", units = units, 
                      ...)
      PRINTFUN(p)
      dev.off()
    }
    if (dev %in% c("jpeg", "jpg")) {
      grDevices::jpeg(filename = file, width = width, height = height, 
                      res = dpi, bg = bg, type = "cairo", units = units, 
                      ...)
      PRINTFUN(p)
      dev.off()
    }
    if (dev %in% c("ps", "cairo_ps")) {
      grDevices::cairo_ps(filename = file, width = get_inches(width, 
                                                              dpi = dpi, units = units), height = get_inches(height, 
                                                                                                             dpi = dpi, units = units), fallback_resolution = dpi, 
                          bg = bg, ...)
      PRINTFUN(p)
      dev.off()
    }
  }
  invisible(p)
}




getRanks <- function(data, rank_by = 1, type = NULL){
  
  # Column 1 is always the main ranking, should be positive and negative
  
  rank_by <- rlang::enquo(rank_by)
  
  
  rankdf <- dplyr::select(data, !!rank_by)
  
  if (!is.null(type)){
    stopifnot(length(type) == ncol(rankdf))
    for (i in seq(type)){
      if (type[i] %in% c("+", "positive")) rankdf[[i]] <- rankdf[[i]] * sign(rankdf[[1]])
      if (type[i] %in% c("p", "probability")) rankdf[[i]] <- -log10(rankdf[[i]]) * sign(rankdf[[1]])
    }
  }
  
  df <- dplyr::arrange_all(rankdf)
  df <- df[rev(1:nrow(df)),,drop = FALSE]
  
  df$id <- as.character(df[[1]])
  df$id[is.na(df$id)] <- "na"
  df$dup <- df$id %in% df$id[duplicated(df$id)]
  df$newstat <- df[[1]]
  
  for (id in unique(subset(df, dup == TRUE)$id)){
    
    tmp <- df[df$id == id,,drop = FALSE]
    ix <- range(which(id == df$id))
    
    tmp2 <- df[c(ix[1]-1, ix[2]+1),,drop = FALSE]
    maxdiff <- min(abs(tmp2[[1]] - mean(tmp[[1]]))) * 0.4
    
    tmpranks <- rev(1:nrow(tmp))
    tmpranks <- (tmpranks - median(tmpranks)) / max(abs(tmpranks))
    
    tmp$newstat <- tmp[[1]] + maxdiff * tmpranks
    df[df$id == id,]$newstat <- tmp$newstat
  }
  
  ranks <- setNames(df$newstat, rownames(df))
  ranks
}



runGSEA <- function(data, genesets = NULL, rank_by = c(stat, baseMean, svalue), type = c("", "+", "p"), ...){
  
  rank_by <- rlang::enquo(rank_by)
  
  stopifnot(requireNamespace("clusterProfiler"))
  stopifnot(requireNamespace("fgsea"))
  stopifnot(requireNamespace("msigdbr"))
  
  if (is.null(genesets)){
    genesets <- NA
  }
  
  if (class(genesets) == "list"){
    
  }
  
  ranks <- getRanks(data, rank_by = !!rank_by, type = type)
  ranks <- sort(ranks, decreasing = TRUE)
  
  results <- clusterProfiler::GSEA(ranks,
                                   seed = 123,
                                   eps = 0,
                                   minGSSize = 1,
                                   maxGSSize = 2000,
                                   TERM2GENE = genesets,
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "fdr",
                                   ...)
  
  results
}





runGSVA <- function(data, genesets = NULL, method = "gsva", kcdf = "Gaussian", ncores = 1, ...){
  
  stopifnot(requireNamespace("GSVA"))
  
  ncores <- min(8, ncores, parallel::detectCores())
  data <- data.matrix(data)
  
  results <- GSVA::gsva(data,
                        genesets,
                        method = "gsva",
                        kcdf = "Gaussian",
                        parallel.sz = ncores,
                        ...)
  
  results <- data.frame(results)
  results
}




getHUMAN1sets <- function(type = "genes"){
  
  stopifnot(type %in% c("genes", "metabolites"))
  
  json_subsystems <- rjson::fromJSON(paste(readLines("https://metabolicatlas.org/api/v2/maps/listing?model=HumanGem"), collapse=""))
  subsystems <- sapply(json_subsystems$subsystems, function(tmp) tmp$id )
  subsystem_data <- lapply(setNames(subsystems, subsystems), function(tmp){
    url <- paste0("https://metabolicatlas.org/api/v2/subsystems/", tmp,"?model=HumanGem")
    rjson::fromJSON(paste(readLines(url, warn = FALSE), collapse=""))
  })
  
  
  if (type == "genes"){
    HUMAN1 <- sapply(subsystem_data, function(tmp){
      unique(unlist( sapply(tmp$genes, function(tmp2){ tmp2$name }) ))
    })
  }
  
  if (type == "metabolites"){
    HUMAN1 <- sapply(subsystem_data, function(tmp){
      unique(unlist( sapply(tmp$metabolites, function(tmp2){ tmp2$id }) ))
    })
  }
  
  names(HUMAN1) <- paste0("HUMAN1_", names(HUMAN1))
  HUMAN1
}



readGTF <- function(file, columns = NULL, ...){
  
  gtf <- rtracklayer::import(file)
  
  gtf <- subset(gtf, ...)
  
  if (is.null(columns)){
    columns <- c("gene_name", "seqnames", "gene_type")
    cat(paste0("Returning columns: ", paste(columns, collapse = ", ")))
    cat(paste0("Available columns: ", paste(colnames(as.data.frame(head(gtf))), collapse = ", ")))
  } else if ("all" %in% tolower(columns)){
    columns <- colnames(as.data.frame(head(gtf)))
  }
  
  df <- as.data.frame(gtf)[,columns, drop = FALSE]
  df <- dplyr::distinct(df)
  df
}



readTables <- function (file, rowNames = TRUE, ...) 
{
  sheets <- openxlsx::getSheetNames(file)
  res <- lapply(setNames(sheets, sheets), function(tmp) openxlsx::read.xlsx(sheet = tmp, 
                                                                            xlsxFile = file, rowNames = rowNames, ...))
  if (length(res) == 1) 
    res <- res[[1]]
  res
}



pval_format <- function(p, min = 0.001, scientific = NULL, add = "p ", stars = FALSE){
  
  if (!is.null(min)){
    digits <- nchar(sub("^-?\\d*\\.?","", min))
    ptext <- round(p, digits = digits)
    poob <- p < min
    ptext[poob] <- min
    ptext <- paste0(add, ifelse(poob, "<", "="), " ", ptext)
    
    psc <- signif(p, digits) |> format(scientific=TRUE)
    if (!is.null(scientific)) ptext[poob] <- psc[poob]
    
  }
  
  # format(p, scientific=TRUE)
  
  ptext
}


convertGeneIDs <- function(ids, from = "ENTREZID", to = "SYMBOL", annotation = org.Hs.eg.db::org.Hs.eg.db, multiVals = "collapse", ...){
  
  stopifnot(requireNamespace("AnnotationDbi"))
  
  ktys <- AnnotationDbi::keytypes(annotation)
  if (!from %in% ktys | !to %in% ktys) {
    print(paste0("Available ID types are: ", paste(ktys, 
                                                   collapse = ", ")))
    stop("Error: Keytype not found!")
  }
  if (multiVals == "collapse") {
    multiVals <- function(x) {
      paste0(x, collapse = "|")
    }
  }
  res <- AnnotationDbi::mapIds(x = annotation, keys = ids, 
                               column = to, keytype = from, multiVals = multiVals, ...)
  res[res == "NA"] <- NA
  res
}



cxheatmap <- function (data, rowdf = NULL, coldf = NULL, scale = FALSE, cluster_rows = NULL, 
          cluster_cols = NULL, rowdf_side = "left", coldf_side = "top", 
          rowdf_legend = TRUE, coldf_legend = TRUE, legend_border = "black", 
          anno_border = "black", fontsize = 12, rowcex = NULL, colcex = 1, 
          rownames_width = 0.3, colnames_width = 0.3, heatpal = NULL, 
          border = NULL, title = NULL, colors = NULL, inf = F, na = 0, 
          mat = NULL, markoob = FALSE, markshape = 4, marksize = NULL, 
          na_col = "grey", maxchar = 35, ...){
  
  datacall <- substitute(data)
  if (is.null(title)) {
    if (grepl("row|col", scale, ignore.case = TRUE)) {
      title <- "z-score"
    }
    else {
      title <- deparse1(datacall)
    }
  }
  if (title == FALSE) 
    title <- " "
  heatdata <- eval(datacall, envir = parent.frame())
  heatdata <- data.matrix(heatdata)
  heatdata <- matScale(heatdata, rows = grepl("row", scale, 
                                              ignore.case = TRUE), cols = grepl("col", scale, ignore.case = TRUE))
  clust <- clusterData(heatdata, rows = cluster_rows, cols = cluster_cols, 
                       inf = inf, na = na)
  if (is.null(clust$rows)) 
    clust$rows <- FALSE
  if (is.null(clust$cols)) 
    clust$cols <- FALSE
  if (is.null(heatpal)) {
    heatpal_colors <- getColorScale(heatdata)
    heatpal <- circlize::colorRamp2(breaks = heatpal_colors, 
                                    colors = names(heatpal_colors))
  }
  docol <- setdiff(unlist(lapply(list(coldf, rowdf), colnames)), 
                   names(colors))
  addcol <- NULL
  if (length(docol) > 0) {
    if (!is.null(coldf)) 
      all <- coldf
    if (!is.null(rowdf)) 
      all <- rowdf
    if (!is.null(coldf) & !is.null(rowdf)) 
      all <- dplyr::full_join(coldf, rowdf, by = character())
    addcol <- getColors(all[, docol, drop = FALSE])
  }
  colors <- c(colors, addcol)
  colors <- lapply(colors, function(tmp) tmp[!is.na(tmp) & 
                                               !is.na(names(tmp))])
  if (is.null(border)) {
    if (nrow(heatdata) < 100 & ncol(heatdata) < 100) {
      border <- grid::gpar(col = rgb(1, 1, 1), lwd = grid::unit(1, 
                                                                "pt"))
    }
    else {
      border <- grid::gpar(col = NA)
    }
  }
  else {
    if (any(border == TRUE)) {
      border <- grid::gpar(col = rgb(1, 1, 1), lwd = grid::unit(1, 
                                                                "pt"))
    }
    else if (length(border) > 1) {
      border <- grid::gpar(col = border[is.na(as.numeric(border))], 
                           lwd = grid::unit(as.numeric(border)[!is.na(as.numeric(border))], 
                                            "pt"))
    }
    else {
      border <- grid::gpar(col = NA)
    }
  }
  legend_params <- list(title_gp = grid::gpar(fontsize = fontsize, 
                                              fontface = "bold"), legend_height = grid::unit(0.2, "npc"), 
                        border = legend_border, labels_gp = grid::gpar(fontsize = fontsize))
  rowAnn <- NULL
  if (!is.null(rowdf)) 
    rowAnn <- getCXanno(df = rowdf[rownames(heatdata), , 
                                   drop = FALSE], colors = colors, anno_border = anno_border, 
                        side = rowdf_side, legend = rowdf_legend, legend_params = legend_params)
  colAnn <- NULL
  if (!is.null(coldf)) 
    colAnn <- getCXanno(coldf[colnames(heatdata), , drop = FALSE], 
                        colors = colors, anno_border = anno_border, side = coldf_side, 
                        legend = coldf_legend, legend_params = legend_params)
  if (markoob == TRUE & is.null(mat)) {
    mat <- matrix(data = FALSE, nrow = nrow(heatdata), ncol = ncol(heatdata), 
                  dimnames = dimnames(heatdata))
    mat[heatdata < min(heatpal_colors)] <- TRUE
    mat[heatdata > max(heatpal_colors)] <- TRUE
  }
  if (is.null(marksize)) 
    marksize <- fontsize * 0.6
  cellFUN <- NULL
  if (!is.null(mat)) {
    if (!is.logical(mat)) 
      stop("'Mat' must be a logical indicator of whether cells should be marked!")
    cellmat <- mat[rownames(heatdata), colnames(heatdata)]
    cellFUN <- function(j, i, x, y, width, height, fill) {
      if (naf(cellmat[i, j] == TRUE)) {
        grid::grid.points(x, y, pch = markshape, size = unit(marksize, 
                                                             "pt"))
      }
    }
  }
  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, 
                                                                      maxchar = maxchar))
  if (is.null(rowcex) & nrow(heatdata) > 10 * ncol(heatdata)) 
    rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) 
    rowcex <- 1
  if (is.null(colcex)) 
    colcex <- 1
  hm <- ComplexHeatmap::Heatmap(name = title, matrix = heatdata, 
                                row_names_max_width = grid::unit(rownames_width, "npc"), 
                                column_names_max_height = grid::unit(colnames_width, 
                                                                     "npc"), column_title_gp = grid::gpar(fontsize = fontsize, 
                                                                                                          fontface = "bold"), rect_gp = border, na_col = na_col, 
                                left_annotation = rowAnn, top_annotation = colAnn, row_names_gp = grid::gpar(fontsize = fontsize * 
                                                                                                               rowcex), column_names_gp = grid::gpar(fontsize = fontsize * 
                                                                                                                                                       colcex), col = heatpal, heatmap_legend_param = legend_params, 
                                cluster_rows = clust$rows, cluster_columns = clust$cols, 
                                cell_fun = cellFUN, ...)
  hm
}


