#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(scattermore); library(scales); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
task <- normalizePath(args[[1]], mustWork=TRUE)
for (d in c("output/pdf","output/png","audit")) dir.create(file.path(task,d),recursive=TRUE,showWarnings=FALSE)
major_order <- c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS","OPC","OLG","MG","EC","Pyr","CR","Per","UN1","UN2","Unmatched")
display <- setNames(major_order,major_order)
display[c("IN","M-AS","AS/qNSC")] <- c("GABA-IN","AS1","AS2/qNSC")
palette <- setNames(c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
  "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4","#8B658B","#FFFF00",
  "#00F5FF","#000000","#7E6148B2","#7FFF00","#B8B8B8"),major_order)
specs <- list(
  list(id="Yao_current_UMAP_author_labels",bundle="Yao_current_UMAP_author_labels.rds",
       title="Yao | author cell types", subtitle="Total: 73,721 cells; 16 author cell types; 174 unmatched"),
  list(id="Yao_author_atlas_original_UMAP_labels",bundle="Yao_author_atlas_original_UMAP_labels.rds",
       title="Yao | original author atlas",subtitle="Total: 92,966 cells; 16 author cell types; original UMAP"))
all_anchors <- list()
all_counts <- list()
audit <- list()
for (spec in specs) {
  b <- readRDS(file.path(task,"plot_data",spec$bundle))
  stopifnot(identical(b$cells,rownames(b$umap)),identical(b$cells,rownames(b$metadata)))
  values <- as.character(b$metadata$AuthorMajor)
  stopifnot(!anyNA(values),all(values %in% major_order))
  types <- major_order[major_order %in% values]
  coords <- b$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  counts <- table(factor(values,levels=types))
  data <- data.frame(x=coords[,1],y=coords[,2],type=factor(values,levels=major_order),row.names=b$cells)
  anchors <- do.call(rbind,lapply(types,function(type) {
    idx <- which(values==type)
    med <- apply(coords[idx,,drop=FALSE],2,median)
    distances <- ((coords[idx,1]-med[1])/diff(xr))^2+((coords[idx,2]-med[2])/diff(yr))^2
    k <- idx[which.min(distances)]
    data.frame(type=type,display_label=display[[type]],n_cells=as.integer(counts[[type]]),
      x=coords[k,1],y=coords[k,2],anchor_cell=b$cells[k],
      label_text=paste0(display[[type]],"\n(",comma(as.integer(counts[[type]])),")"))
  }))
  # Optional explicit nudges affect text only, never cells or coordinates.
  overrides_path <- file.path(task,"code/label_position_overrides.csv")
  anchors$nudge_x <- 0; anchors$nudge_y <- 0
  if (file.exists(overrides_path)) {
    ov <- read.csv(overrides_path,stringsAsFactors=FALSE)
    ov <- ov[ov$plot_id==spec$id,,drop=FALSE]
    for (j in seq_len(nrow(ov))) {
      k <- match(ov$type[j],anchors$type)
      stopifnot(!is.na(k))
      anchors$nudge_x[k] <- ov$nudge_x[j]; anchors$nudge_y[k] <- ov$nudge_y[j]
    }
  }
  p <- ggplot(data,aes(x=x,y=y,colour=type))+
    scattermore::geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    ggrepel::geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      colour="black",size=3.7,lineheight=0.95,seed=20261016,
      nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=grid::unit(0.45,"lines"),point.padding=grid::unit(0.15,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=0.25,
      force=2,force_pull=0.3,max.overlaps=Inf,max.time=Inf,max.iter=10000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xr,ylim=yr,clip="off")+
    labs(title=spec$title,subtitle=spec$subtitle,x="UMAP_1",y="UMAP_2")+
    theme_classic(base_size=12)+
    theme(plot.title=element_text(face="bold",size=19),plot.subtitle=element_text(size=11,colour="#333333"),
      axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=11),
      legend.position="none",plot.margin=margin(12,25,12,15))
  ggsave(file.path(task,"output/pdf",paste0(spec$id,".pdf")),plot=p,device=grDevices::cairo_pdf,
         width=12.5,height=10,units="in",bg="white")
  ggsave(file.path(task,"output/png",paste0(spec$id,".png")),plot=p,device=ragg::agg_png,
         width=12.5,height=10,units="in",dpi=240,bg="white")
  anchors$plot_id <- spec$id
  all_anchors[[spec$id]] <- anchors
  all_counts[[spec$id]] <- data.frame(plot_id=spec$id,author_Major=types,display_label=unname(display[types]),
    n_cells=as.integer(counts),colour=unname(palette[types]))
  audit[[spec$id]] <- list(n_cells=nrow(data),n_author_types=sum(types!="Unmatched"),
    unmatched_n=sum(values=="Unmatched"),n_labels=nrow(anchors),coordinate_ratio=1,
    x_range=xr,y_range=yr,labels_from_author_Major=TRUE,source_bundle=spec$bundle,
    coordinates_identical_to_bundle=TRUE,umap_recomputed=FALSE,
    labels_have_background=FALSE,pointsize=2.2,width_in=12.5,height_in=10,png_dpi=240)
  message("PLOTTED ",spec$id)
}
write.csv(do.call(rbind,all_anchors),file.path(task,"audit/plot_label_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,all_counts),file.path(task,"audit/plot_counts_and_colours.csv"),row.names=FALSE)
write_json(audit,file.path(task,"audit/plot_preservation.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("AUTHOR_LABEL_PLOTTING_COMPLETE")
