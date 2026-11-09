#!/usr/bin/env Rscript
# Read-only expression visualization. Never resave or relabel a source object.
suppressPackageStartupMessages({library(qs); library(Seurat); library(Matrix); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
out <- args[[1]]; mode <- args[[2]]
for (d in c("tables", "audit", "pdf", "png", "runtime")) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
genes <- c("CSF1R","P2RY12","AIF1","CX3CR1","TMEM119","OLFML3","GPR34","HEXB","SALL1","GPNMB","SPP1","APOE","TREM2","CST3")
studies <- c("Su","Chia","Fullard","Jakel","Grubman","Yao")
paths <- c("Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs", "Chia_et_al./data/processed_data/Chia.qs", "Fullard_et_al./data/processed_data/Fullard.qs", "Jakel_et_al./data/processed_data/Jakel.qs", "Grubman_et_al./data/processed_data/Grubman.qs", "Yao_et_al./data/processed_data/Yao_v2.qs")
columns <- c("CellType","cell.type","celltype","Celltypes","cellType","CellType")
expected_n <- c(224464L,87470L,68557L,17644L,13214L,73721L)
notes <- c("Existing Su CellType; source annotation retained", "Source cell.type, not project CellType", "Source celltype", "Source Celltypes (23 labels); 17,644 cells, not the 17,799-cell v2 object", "Source cellType, including unID and doublet; not project CellType", "Existing project CellType; retains Unknown1/Unknown2; not claimed article-native")
manifest <- data.frame(study=studies,source_path=file.path(root,paths),annotation_column=columns,expected_cells=expected_n,annotation_note=notes,stringsAsFactors=FALSE)
write.csv(manifest,file.path(out,"audit/source_manifest.csv"),row.names=FALSE)
if (mode=="extract") {
  inventory <- list(); gene_inventory <- list()
  for (i in seq_along(studies)) {
    study <- studies[[i]]; path <- manifest$source_path[[i]]; field <- columns[[i]]
    message(format(Sys.time(),"%F %T")," Reading ",study,": ",path)
    obj <- qread(path,nthreads=1)
    stopifnot(ncol(obj)==expected_n[[i]], field %in% colnames(obj[[]]), "RNA" %in% Assays(obj))
    cnt <- GetAssayData(obj,assay="RNA",layer="counts")
    stopifnot(identical(colnames(cnt),colnames(obj)), nrow(cnt)>1000L)
    if (any(!is.finite(cnt@x)) || any(cnt@x<0)) stop(study,": invalid RNA counts")
    lib <- Matrix::colSums(cnt)
    if (any(!is.finite(lib)) || any(lib<=0)) stop(study,": invalid RNA library totals")
    labels <- as.character(obj[[]][[field]])
    labels[is.na(labels)|!nzchar(labels)] <- "<missing>"
    original <- obj[[]][[field]]
    lev <- if (is.factor(original)) levels(droplevels(original)) else sort(unique(labels))
    lev <- unique(c(lev[lev %in% labels],setdiff(unique(labels),lev)))
    grp <- factor(labels,levels=lev)
    features <- vapply(genes,function(g) {hit <- rownames(cnt)[toupper(rownames(cnt))==g]; if(length(hit)>1L) stop(study,": duplicated feature ",g); if(length(hit)) hit else NA_character_},character(1))
    present <- genes[!is.na(features)]
    selected <- cnt[features[present],,drop=FALSE]
    rownames(selected) <- present
    cp10k <- selected %*% Diagonal(x=10000/lib)
    dimnames(cp10k) <- dimnames(selected)
    lognorm <- cp10k; lognorm@x <- log1p(lognorm@x)
    # Seurat DotPlot on a plotting-only 14-gene object, normalized using all RNA genes.
    thin <- CreateSeuratObject(counts=selected,meta.data=data.frame(native_celltype=grp,row.names=colnames(cnt)),min.cells=0,min.features=0)
    thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=lognorm)
    dp <- DotPlot(thin,features=present,assay="RNA",group.by="native_celltype",scale=TRUE,col.min=-2.5,col.max=2.5,dot.min=0)$data
    rows <- list()
    for (g in genes) {
      vals <- if(g %in% present) as.numeric(cp10k[g,]) else NULL
      for (label in lev) {
        idx <- which(labels==label); available <- g %in% present
        z <- if(available) dp$avg.exp.scaled[as.character(dp$features.plot)==g & as.character(dp$id)==label] else NA_real_
        pct <- if(available) 100*mean(vals[idx]>0) else NA_real_
        mu <- if(available) mean(vals[idx]) else NA_real_
        if(available) {
          q <- dp[as.character(dp$features.plot)==g & as.character(dp$id)==label,,drop=FALSE]
          stopifnot(nrow(q)==1L,abs(q$pct.exp-pct)<1e-8,abs(q$avg.exp-mu)<1e-8)
        }
        rows[[length(rows)+1L]] <- data.frame(study=study,annotation_column=field,celltype=label,gene=g,actual_feature=features[[g]],n_cells=length(idx),gene_present=available,percent_expressed=pct,mean_normalized_CP10K=mu,log1p_mean_CP10K=if(available)log1p(mu) else NA_real_,mean_log1p_CP10K=if(available)mean(log1p(vals[idx])) else NA_real_,scaled_average_expression=z,stringsAsFactors=FALSE)
      }
    }
    result <- do.call(rbind,rows)
    stopifnot(nrow(result)==length(genes)*length(lev),sum(table(labels))==ncol(obj))
    write.csv(result,file.path(out,"tables",paste0(study,"_microglia_markers_by_source_celltype.csv")),row.names=FALSE,na="NA")
    write.csv(data.frame(study=study,annotation_column=field,celltype=lev,n_cells=as.integer(table(grp))),file.path(out,"tables",paste0(study,"_source_celltype_cell_counts.csv")),row.names=FALSE)
    inventory[[study]] <- data.frame(study=study,source_path=path,source_size=file.info(path)$size,source_mtime=as.character(file.info(path)$mtime),annotation_column=field,n_cells=ncol(obj),n_celltypes=length(lev),n_RNA_features=nrow(cnt),n_markers_present=length(present),missing_markers=paste(setdiff(genes,present),collapse=";"),n_missing_annotations=sum(labels=="<missing>"),max_fractional_count_deviation=max(abs(cnt@x-round(cnt@x))),stringsAsFactors=FALSE)
    gene_inventory[[study]] <- data.frame(study=study,gene=genes,actual_feature=unname(features),gene_present=!is.na(features),stringsAsFactors=FALSE)
    write.csv(do.call(rbind,inventory),file.path(out,"audit/object_and_annotation_inventory.csv"),row.names=FALSE)
    write.csv(do.call(rbind,gene_inventory),file.path(out,"audit/feature_availability.csv"),row.names=FALSE,na="NA")
    message(study,": ",ncol(obj)," cells; ",length(lev)," source labels; ",length(present),"/14 markers; table validated against Seurat DotPlot")
    rm(obj,cnt,thin,cp10k,lognorm,selected,result,dp); gc(verbose=FALSE)
  }
  writeLines(capture.output(sessionInfo()),file.path(out,"audit/extraction_sessionInfo.txt"))
  message("EXTRACTION_COMPLETE")
} else if(mode=="plot") {
  plots <- list(); dimensions <- list()
  for (study in studies) {
    dat <- read.csv(file.path(out,"tables",paste0(study,"_microglia_markers_by_source_celltype.csv")),check.names=FALSE)
    info <- read.csv(file.path(out,"audit/object_and_annotation_inventory.csv"),check.names=FALSE)
    info <- info[info$study==study,,drop=FALSE]
    lev <- unique(dat$celltype)
    dat$gene <- factor(dat$gene,levels=genes)
    dat$celltype <- factor(dat$celltype,levels=rev(lev))
    p <- ggplot(dat,aes(x=gene,y=celltype)) +
      geom_hline(yintercept=seq_along(lev),colour="#F0F0F0",linewidth=0.3) +
      geom_vline(xintercept=seq_along(genes),colour="#F6F6F6",linewidth=0.3) +
      geom_point(data=dat[dat$gene_present,,drop=FALSE],aes(size=percent_expressed,colour=scaled_average_expression)) +
      scale_size_area(name="Percent expressed",limits=c(0,100),breaks=c(0,25,50,75,100),max_size=6.8) +
      scale_colour_gradientn(name="Scaled average\nexpression",colours=c("#F0F0F0","#FEE8A0","#FC8D59","#D7301F","#99000D"),values=c(0,0.25,0.5,0.75,1),limits=c(-2.5,2.5),breaks=c(-2,-1,0,1,2)) +
      scale_x_discrete(drop=FALSE,expand=expansion(add=0.6)) + scale_y_discrete(drop=FALSE,expand=expansion(add=0.65)) +
      labs(title=paste0(study,": microglial marker expression"),subtitle=paste0(format(info$n_cells,big.mark=",")," cells; ",info$n_celltypes," existing cell-type labels (",info$annotation_column,")"),x=NULL,y=NULL,
           caption="Colour: per-gene z-score of log1p(mean CP10K), calculated within this study.\nDot area: cells with RNA counts > 0 (%). No new filtering, clustering or annotation.") +
      guides(size=guide_legend(order=1,override.aes=list(colour="#333333")),colour=guide_colourbar(order=2,barheight=grid::unit(35,"mm"))) +
      theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10,colour="#222222"),axis.text.y=element_text(size=10,colour="#222222"),axis.ticks=element_blank(),axis.line=element_line(linewidth=0.4),legend.position="right",legend.text=element_text(size=9),legend.title=element_text(size=10),plot.title=element_text(face="bold",size=16),plot.subtitle=element_text(size=10.5,margin=margin(b=10)),plot.caption=element_text(hjust=0,size=8.5,colour="#4D4D4D",margin=margin(t=12)),plot.margin=margin(12,16,12,12),panel.border=element_rect(colour="#999999",fill=NA,linewidth=0.4))
    if(any(!dat$gene_present)) p <- p + geom_point(data=dat[!dat$gene_present,,drop=FALSE],shape=4,size=2,colour="#777777") + labs(caption="X: gene absent from this object's feature list, not zero expression.\nColour: within-study per-gene z-score of log1p(mean CP10K). Dot area: RNA count > 0 (%).")
    h <- max(5.2,3.25+0.27*length(lev)); w <- if(max(nchar(lev))>28) 13.8 else 12.4
    pdf <- file.path(out,"pdf",paste0(study,"_microglia_markers_source_celltype_DotPlot.pdf"))
    png <- file.path(out,"png",paste0(study,"_microglia_markers_source_celltype_DotPlot.png"))
    ggsave(pdf,p,device=cairo_pdf,width=w,height=h,bg="white",limitsize=FALSE)
    ggsave(png,p,device=ragg::agg_png,width=w,height=h,dpi=240,bg="white",limitsize=FALSE)
    plots[[study]] <- p
    dimensions[[study]] <- data.frame(study=study,width_inches=w,height_inches=h,n_rows=length(lev),n_genes=length(genes))
    message(study," PDF/PNG saved")
  }
  cairo_pdf(file.path(out,"Six_study_microglia_markers_source_celltype_DotPlots.pdf"),width=13.8,height=10,onefile=TRUE)
  for(p in plots) print(p)
  dev.off()
  write.csv(do.call(rbind,dimensions),file.path(out,"audit/plot_dimensions.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"audit/plotting_sessionInfo.txt"))
  message("PLOTTING_COMPLETE")
} else stop("mode must be extract or plot")
