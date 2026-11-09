#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
base <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results")
current_path <- file.path(base, "source_celltype_umap_20260920/revisions/20261016_celltype_labels_on_umap_no_boxes/plot_data/Yao_labels_umap.rds")
qc_path <- file.path(base, "qc_per_cell_nCount_nFeature_by_source_celltype_20260920/per_study_csv/Yao_per_cell_nCount_nFeature.csv")
author <- readRDS(file.path(task, "inputs/author_metadata_umaps.rds"))
current <- readRDS(current_path)
dir.create(file.path(task, "plot_data"), showWarnings = FALSE)
stopifnot(length(author$cells) == 92966L, length(current$cells) == 73721L,
          identical(rownames(author$metadata), author$cells),
          identical(rownames(current$umap), current$cells),
          identical(rownames(current$metadata), current$cells))
parse_cells <- function(cells) {
  valid <- grepl("^.+_[ACGT]{16}-[0-9]+$", cells)
  stopifnot(all(valid))
  barcode <- sub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", cells)
  sample <- sub("_[ACGT]{16}-[0-9]+$", "", cells)
  data.frame(cell_id = cells, sample = sample, barcode = barcode, stringsAsFactors = FALSE)
}
source_ids <- parse_cells(author$cells)
target_ids <- parse_cells(current$cells)
stopifnot(identical(source_ids$sample, as.character(author$metadata$orig.ident)))
aliases <- data.frame(target_sample = c("l_4D", "l_32y", "l_48y", "l_50y", "l_56y", "l_60y", "l_64f", "l_64m", "l_68y"),
                      author_sample = c("4D", "32y", "48y", "50y", "56y", "60y", "64f", "64m", "68y"))
stopifnot(setequal(unique(target_ids$sample), aliases$target_sample))
target_ids$author_sample <- aliases$author_sample[match(target_ids$sample, aliases$target_sample)]
source_key <- paste(source_ids$sample, source_ids$barcode, sep = "|")
target_key <- paste(target_ids$author_sample, target_ids$barcode, sep = "|")
stopifnot(!anyDuplicated(source_key), !anyDuplicated(target_key))
idx <- match(target_key, source_key)
matched <- !is.na(idx)
stopifnot(sum(matched) == 73547L, sum(!matched) == 174L)
order <- c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS","OPC","OLG","MG","EC","Pyr","CR","Per","UN1","UN2")
expected_original <- c(GC=24671, IN=8601, NB=2607, aNSC=2140, pNSC=2798,
  "AS/qNSC"=11071, "M-AS"=1146, OPC=5396, OLG=15796, MG=11823, EC=1232, Pyr=676, CR=218, Per=981, UN1=1999, UN2=1811)
expected_current <- c(GC=17166, IN=4566, NB=1770, aNSC=2082, pNSC=2724,
  "AS/qNSC"=10130, "M-AS"=1045, OPC=4504, OLG=15270, MG=9495, EC=935, Pyr=231, CR=147, Per=715, UN1=1341, UN2=1426)
display <- setNames(order, order)
display[c("IN","M-AS","AS/qNSC")] <- c("GABA-IN","AS1","AS2/qNSC")
major <- as.character(author$metadata$Major)
stopifnot(setequal(unique(major), order), all(as.integer(table(factor(major, levels=order))) == expected_original[order]))
current_major <- rep("Unmatched", length(idx))
current_major[matched] <- major[idx[matched]]
stopifnot(all(as.integer(table(factor(current_major[matched], levels=order))) == expected_current[order]))
display <- c(display, Unmatched = "Unmatched")
definitions <- c(GC="Granule cell", IN="GABAergic interneuron", NB="Neuroblast",
  aNSC="Active neural stem cell", pNSC="Primed neural stem cell",
  "AS/qNSC"="Astrocyte2 / quiescent neural stem cell", "M-AS"="Astrocyte1",
  OPC="Oligodendrocyte progenitor", OLG="Oligodendrocyte", MG="Microglia",
  EC="Endothelial cell", Pyr="Pyramidal neuron", CR="Cajal-Retzius cell",
  Per="Pericyte", UN1="Unidentified population 1", UN2="Unidentified population 2")
definition_table <- data.frame(author_Major=order, display_label=unname(display[order]),
  definition=unname(definitions[order]), author_atlas_n=as.integer(expected_original[order]),
  current_umap_n=as.integer(expected_current[order]), source="Author Major; eLife 89507 Figure 1",
  paper_count_basis=ifelse(order %in% c("UN1","UN2"),
    "Paper reports combined UN1+UN2=3810; separate counts from author Major",
    "Individual author count agrees with paper Figure 1 legend"))
write.csv(definition_table, file.path(task, "audit/author_definitions_and_counts.csv"), row.names=FALSE)
write.csv(aliases, file.path(task, "audit/sample_aliases.csv"), row.names=FALSE)
mapping <- data.frame(current_cell_id=current$cells, target_sample=target_ids$sample,
  author_sample=target_ids$author_sample, barcode=target_ids$barcode,
  author_cell_id=ifelse(matched, author$cells[idx], NA_character_),
  old_CellType=as.character(current$metadata$CellType), author_Major=current_major,
  display_label=unname(display[current_major]), author_label_available=matched)
write.csv(mapping, file.path(task, "audit/current_cell_author_label_mapping.csv"), row.names=FALSE)
write.csv(mapping[!matched, ], file.path(task, "audit/unmatched_current_cells.csv"), row.names=FALSE)
cross <- as.data.frame(table(old_CellType=mapping$old_CellType, author_Major=mapping$author_Major))
write.csv(cross[cross$Freq>0, ], file.path(task, "audit/old_CellType_vs_author_Major.csv"), row.names=FALSE)
by_sample <- do.call(rbind, lapply(seq_len(nrow(aliases)), function(i) {
  j <- target_ids$sample == aliases$target_sample[[i]]
  data.frame(aliases[i, ], current_n=sum(j), author_n=sum(source_ids$sample==aliases$author_sample[[i]]),
             matched_n=sum(matched[j]), unmatched_n=sum(!matched[j]))
}))
write.csv(by_sample, file.path(task, "audit/barcode_match_by_sample.csv"), row.names=FALSE)
qc <- read.csv(qc_path, stringsAsFactors=FALSE, check.names=FALSE)
qi <- match(current$cells, qc$cell_id)
stopifnot(!anyNA(qi), !anyDuplicated(qc$cell_id), nrow(qc)==length(current$cells))
d_count <- qc$nCount_RNA[qi[matched]]-author$metadata$nCount_RNA[idx[matched]]
d_feature <- qc$nFeature_RNA[qi[matched]]-author$metadata$nFeature_RNA[idx[matched]]
stopifnot(all(d_count<=0), all(d_feature<=0), min(d_count)== -22, min(d_feature)== -20,
          sum(d_count==0 & d_feature==0)==60044L)
write.csv(data.frame(current_cell_id=current$cells[matched], nCount_difference=d_count, nFeature_difference=d_feature),
          file.path(task, "audit/matched_cell_QC_identity_check.csv"), row.names=FALSE)
metadata_current <- data.frame(AuthorMajor=current_major, AuthorLabel=unname(display[current_major]),
                              OldCellType=mapping$old_CellType, row.names=current$cells)
metadata_original <- data.frame(AuthorMajor=major, AuthorLabel=unname(display[major]), row.names=author$cells)
current_plot <- list(cells=current$cells, umap=current$umap, metadata=metadata_current,
                     original_bundle=current_path, annotation_source=author$source)
original_plot <- list(cells=author$cells, umap=author$reductions$umap, metadata=metadata_original,
                      annotation_source=author$source)
stopifnot(identical(rownames(original_plot$umap), author$cells),
          identical(current_plot$umap, current$umap), all(is.finite(current_plot$umap)),
          all(is.finite(original_plot$umap)))
saveRDS(current_plot, file.path(task, "plot_data/Yao_current_UMAP_author_labels.rds"))
saveRDS(original_plot, file.path(task, "plot_data/Yao_author_atlas_original_UMAP_labels.rds"))
audit <- list(author_annotation="Major", source_n=length(author$cells), target_n=length(current$cells),
  matched_n=sum(matched), unmatched_n=sum(!matched), matched_fraction=mean(matched),
  sample_and_barcode_keys_unique=TRUE, exact_original_names_match=0,
  match_method="Explicit nine sample aliases plus complete 10X barcode; no expression or coordinate inference",
  original_16_class_counts_match_paper=TRUE, current_16_class_counts_independently_verified=TRUE,
  paper_count_validation="14 individual class counts and combined UN1+UN2=3810 agree; separate UN1/UN2 counts are from author Major",
  paper_unknown_count_reported_as_combined=TRUE,
  current_coordinates_identical=identical(current_plot$umap,current$umap),
  original_coordinates_identical=identical(original_plot$umap,author$reductions$umap),
  author_sample_missing_from_current="31y", missing_sample_original_n=12406,
  qc_source=qc_path, matched_QC_both_exact=sum(d_count==0 & d_feature==0),
  target_minus_author_nCount_range=range(d_count), target_minus_author_nFeature_range=range(d_feature),
  umap_recomputed=FALSE, source_objects_modified=FALSE, unmatched_labels_inferred=FALSE)
write_json(audit, file.path(task, "audit/annotation_mapping_validation.json"), pretty=TRUE, auto_unbox=TRUE)
print(audit)
print(definition_table)
message("AUTHOR_LABEL_PREPARATION_COMPLETE")
