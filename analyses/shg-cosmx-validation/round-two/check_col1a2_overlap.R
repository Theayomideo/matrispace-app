library(data.table)
library(hexbin)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript check_col1a2_overlap.R analysis_root [output_dir]")
root <- args[1]
out <- if (length(args) > 1L) args[2] else "results"
obj <- readRDS(file.path(root, "derived/section04_matrispace_collagen_scored.rds"))
expr <- SeuratObject::LayerData(obj, assay="SCT", layer="data")
cells <- fread(file.path(root,"derived/section04_cell_collagen_scores.csv.gz"))
hex <- fread(file.path(root,"results/primary_hex_data_area_matched.csv.gz"))
baselines <- fread(file.path(out,"shg_comparator_bins.csv"))
genes <- strsplit(fread(file.path(root,"derived/score_signature_coverage.csv"))[
  signature=="collagens", measured],";",fixed=TRUE)[[1]]
stopifnot(length(genes)==25,"COL1A2" %in% genes,nrow(expr)==960,
          as.character(packageVersion("UCell"))=="2.14.0")
scores <- UCell::ScoreSignatures_UCell(expr,features=list(full=genes,
    omit=setdiff(genes,"COL1A2")), name="", maxRank=960,
    BPPARAM=BiocParallel::SerialParam())
scores <- as.data.frame(scores[match(cells$barcode,rownames(scores)),,drop=FALSE])
stopifnot(identical(cells$barcode,rownames(scores)),
          max(abs(scores$full - cells$collagens_ucell)) < 1e-10)
grid <- hexbin(cells$x_um,cells$y_um,xbins=37L,xbnds=c(0,4000),
               ybnds=c(0,4000),shape=1,IDs=TRUE)
avg <- data.table(hex_id=grid@cID,full=scores$full,omit=scores$omit)[,
  .(n=.N,full=mean(full),omit=mean(omit)),by=hex_id]
h <- merge(merge(hex,avg,by="hex_id"),baselines[,.(hex_id,fibroblast_density)],by="hex_id")
stopifnot(nrow(h)==1545,all(h$n==h$n_cells),max(abs(h$full-h$collagens_ucell))<1e-10)
rho <- function(a,b) cor(a,b,method="spearman")
calc <- function(i) {
  shg <- h$shg_percentile[i]
  den <- h$fibroblast_density[i]
  full <- h$cosmx_percentile[i]
  omit <- h$omit[i]
  c(full_shg=rho(full,shg),omit_shg=rho(omit,shg),density_shg=rho(den,shg),
    full_density=rho(full,den),omit_density=rho(omit,den),full_omit=rho(full,omit),
    delta_shg=rho(omit,shg)-rho(full,shg),delta_density=rho(omit,den)-rho(full,den),
    omit_minus_density_shg=rho(omit,shg)-rho(den,shg))
}
obs <- calc(seq_len(nrow(h)))
stopifnot(abs(obs["full_shg"]-0.489245181950197)<1e-10)
blocks <- split(seq_len(nrow(h)),paste(floor(h$x_um/400),floor(h$y_um/400),sep="_"))
set.seed(15240532)
boot <- replicate(1000,calc(unlist(blocks[sample(names(blocks),length(blocks),replace=TRUE)],
                                   use.names=FALSE)))
ci <- t(apply(boot,1,quantile,probs=c(.025,.975),names=FALSE,type=8))
tab <- data.frame(contrast=names(obs),estimate=as.numeric(obs),lower=ci[,1],upper=ci[,2])
print(tab,digits=5)
fwrite(tab,file.path(out,"col1a2_overlap.csv"))
