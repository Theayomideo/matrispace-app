# Run from the app repository root.
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(UCell)})
if (file.exists('helpers.R')) {
  app <- new.env(parent=globalenv())
  source('helpers.R',local=app,encoding='UTF-8')
  lines <- readLines('global.R',encoding='UTF-8')
  eval(parse(text=lines[seq(grep('^ecm_ucell_signatures <-',lines),length(lines))]),app)
  ref <- as.list(app)
  pairs <- setNames(ref$lr_db,c('Gene1','Gene2','Interaction_Type','Source','Database_Score'))
} else {
  pkgload::load_all('.',quiet=TRUE)
  app <- asNamespace('matrispace.app')
  ref <- app$initialize_reference_data()
  pairs <- ref$matrisome_pairs
}
stopifnot(length(ref$mobj)==1027, nrow(pairs)==16876,
          !any(pairs$Gene1==pairs$Gene2),
          !anyDuplicated(paste(pairs$Gene1,pairs$Gene2)),
          length(intersect(unique(c(pairs$Gene1,pairs$Gene2)),ref$mobj))==929,
          !any(c('F8','USH2A') %in% unlist(ref$matrisome_feature_signatures[1:6])),
          all(c('F8','USH2A') %in% ref$matrisome_feature_signatures$perivascular),
          'F8' %in% ref$matrisome_feature_signatures$hemostasis,
          'USH2A' %in% ref$matrisome_feature_signatures$basement_membrane)
cat('PASS canonical reference, pair coverage and retained functional signatures\n')
counts <- matrix(c(1,2,3,4,5,6,7,8),2,
                 dimnames=list(c('LGALSL','COL1A1'),paste0('s',1:4)))
obj <- CreateSeuratObject(counts)
obj <- SetAssayData(obj,layer='data',new.data=log1p(counts))
obj[['PAIR']] <- CreateAssayObject(data=counts*10)
obj$group <- c('1 A','1 A','B-cells','B-cells')
means <- app$mean_matrisome_pair_scores(obj,'PAIR','group')
stopifnot(identical(colnames(means),c('1 A','B-cells')),
          isTRUE(all.equal(unname(means),matrix(c(20,30,60,70),2))))
obj$group <- 'one'
stopifnot(identical(dim(app$mean_matrisome_pair_scores(obj,'PAIR','group')),c(2L,1L)))
cat('PASS arithmetic means, group labels and single-group shape\n')
translated <- app$translate_gene_signatures(obj,list(alias='HSPC159',missing='MISSING'))
stopifnot(identical(translated,list(alias='LGALSL')))
obj$collagens_score <- 999
obj[['Collagens']] <- 999
scored <- app$add_matrisome_feature_scores(obj,list(collagens='COL1A1'))
expected <- UCell::ScoreSignatures_UCell(GetAssayData(obj,layer='data'),
             features=list(collagens='COL1A1'), BPPARAM=BiocParallel::SerialParam())
stopifnot(isTRUE(all.equal(unname(scored$collagens_score),unname(expected[,1]))),
          identical(unname(scored[['Collagens']][,1]),unname(scored$collagens_score)))
cat('PASS signature translation and cached score refresh\n')

missing <- app$add_matrisome_feature_scores(scored,list(collagens="MISSING"))
stopifnot(!any(c("collagens_score","Collagens") %in% names(missing[[]])))
cat("PASS unavailable signatures do not retain stale scores\n")
