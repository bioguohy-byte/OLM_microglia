#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Supplementary Figure 8

# S8A  Three representative TEM fields each for CC-saline, CC-HMC3 and CC-HMC3-mye
# Kind: experimental_micrographs
mg_register_panel("S8A", NULL)
