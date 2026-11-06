# Run from the repository root.
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
root <- normalizePath('.',winslash='/')
pkgload::load_all(root,quiet=TRUE)
off <- asNamespace('matrispace.app')
for (expr in parse(file.path(root,'dev/check_figure_exports.R'))) {
  if(is.call(expr) && identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name('fixture'))) eval(expr)
}
obj <- fixture()
ref <- off$initialize_reference_data()
set.seed(61)
genes <- unique(c(ref$matrisome$gene,paste0('OTHER',1:50)))
counts <- matrix(rpois(length(genes)*ncol(obj),5),length(genes),
                 dimnames=list(genes,colnames(obj)))
obj[['Spatial']] <- CreateAssayObject(counts=counts)
obj$cell_type <- rep(c('A','B'),length.out=ncol(obj))
obj <- SCTransform(obj,assay='Spatial',verbose=FALSE)
outdir <- tempfile('matrispace_integration_')
dir.create(outdir,showWarnings=FALSE)
check <- function(ok,label) {stopifnot(isTRUE(ok));cat('PASS',label,'\n')}
shiny::testServer(off$app_server, {
  session$setInputs(annot='cell_type',sel1='matrisome gene',sel2='matrisome gene',
                   gene1='COL1A1',gene2='DCN',selectedCellTypes_main=c('A','B'),
                   export_preset='publication',violin_plot=TRUE,independent_y=FALSE,
                   coex_palette='Classic',coex_alpha=0.7,coex_pt_size=3,
                   main_pt_size=2,subcategories_pt_size=2,families_pt_size=2,ecm_pt_size=2,
                   leak_correct=FALSE,leak_alpha=0.1,leak_k=6)
  prepared <- off$prepare_uploaded_object(obj,ref$matrisome_feature_signatures,ref$ecm_ucell_signatures)
  check('ECM Glycoproteins' %in% names(prepared[[]]),'fresh preprocessing completes')
  check('ecm_domain_annotation' %in% names(prepared[[]]),'offline ScType annotation completes')
  rv$active_seurat_object <- prepared
  rv$active_sample_name <- 'integration'
  rv$data_is_loaded <- TRUE
  session$flushReact()
  session$setInputs(run_button=1)
  check(isTRUE(rv$feature_analysis_done),'feature analysis through server')
  session$setInputs(run_matrisome_profile=1)
  check(isTRUE(rv$matrisome_analysis_done),'matrisome profiling through server')
  session$setInputs(run_matrisome_pair_analysis=1)
  check(!is.null(rv$matrisome_pair_results) && all(is.finite(rv$matrisome_pair_results$means)), 'pair analysis through server')
  snap <- capture_export_snapshot()
  check(identical(snap$feature$annotation_col,'cell_type') && !is.null(snap$feature$sel1),'snapshot records feature inputs')
  for(preset in c('publication','current')) {
    snap$preset <- preset
    for(mode in c('all','plots')) {
      file <- file.path(outdir,paste0(preset,'_',mode,'.zip'))
      status <- if(mode=='all') build_export_all(file,snap) else build_export_plots(file,snap)
      check(!status$has_errors,paste(preset,mode,'ZIP has no failed figures'))
      files <- zip::zip_list(file)$filename
      check(all(vapply(c('export_manifest.csv','export_settings.json','README.txt'),function(x) any(endsWith(files,x)),logical(1))),paste(preset,mode,'manifest files'))
      check(all(vapply(c('.svg','.png','.pdf'),function(x) any(endsWith(files,x)),logical(1))),paste(preset,mode,'three figure formats'))
      if(mode=='all') check(any(endsWith(files,'spot_metadata.csv')) && any(endsWith(files,'matrisome_scores.csv')), 'all ZIP includes metadata and scores')
    }
  }
  session$setInputs(prepare_export_plots=1)
  check(file.exists(prepared_export_plots()$path),'prepared download succeeds')
  serve_prepared_export(prepared_export_plots,file.path(outdir,'download.zip'))
  check(file.exists(file.path(outdir,'download.zip')),'prepared download copied successfully')
  corrected <- off$matrispace_spotclean(prepared,genes=ref$matrisome$gene)
  check(corrected$applied && corrected$n_genes>0,'leakage correction runs')
})
cat('All integration checks passed; artifacts:', outdir, '\n')
