#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")
