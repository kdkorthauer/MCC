library(XenofilteR)
library(BiocParallel)

bamdir <- "../../../PREPROCESS/DNA/pdx_filter"
dir.create(file.path(bamdir, "xenofilteR"), showWarnings = FALSE)
bams_human <- list.files(path = bamdir, pattern = "human", full.names = TRUE)
bams_mouse <- list.files(path = bamdir, pattern = "mouse", full.names = TRUE)

bams_mouse <- bams_mouse[!grepl("tmp", bams_mouse)]

x <- match(gsub("_human", "", bams_human),
           gsub("_mouse", "", bams_mouse))

sample.tab <- cbind(bams_human[!is.na(x)],
    bams_mouse[x[!is.na(x)]])

bp.param <- MulticoreParam(1)

# loop through sample pairs

for (i in seq_len(nrow(sample.tab))) {
  if ( ! file.exists()) {
    XenofilteR(sample.tab[2, , drop = FALSE],
        destination.folder = file.path(bamdir, "xenofilteR"),
        bp.param = bp.param)
  }
}