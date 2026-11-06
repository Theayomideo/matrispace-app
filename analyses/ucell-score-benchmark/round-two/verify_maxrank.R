library(data.table)
# Independently reconstruct both sensitivity metrics from the saved per-spot scores.
args <- commandArgs(trailingOnly = TRUE)
bench <- if (length(args)) args[1] else "results"
rho <- function(a, b) cor(a, b, method = "spearman")
scores <- fread(file.path(bench,"strength_maxrank_spot_scores.csv.gz"))
print(names(scores))
ref <- fread(file.path(bench,"strength_maxrank_comparison.csv"))
spotcol <- intersect(c("spot","spot_id","barcode"),names(scores))
scorecol <- intersect(c("score","ucell","ucell_score"),names(scores))
stopifnot(length(spotcol)==1,length(scorecol)==1)
for (i in seq_len(nrow(ref))) {
  r <- ref[i]
  s <- scores[signature == r$signature]
  a <- s[maxRank==1500]
  b <- s[maxRank==r$tested_maxRank][match(a[[spotcol]],get(spotcol))]
  stopifnot(nrow(a)==4992,identical(a[[spotcol]],b[[spotcol]]),!anyDuplicated(a[[spotcol]]))
  top <- function(d) d[[spotcol]][head(order(d[[scorecol]],decreasing=TRUE),ceiling(nrow(d)*0.1))]
  ta <- top(a); tb <- top(b)
  stopifnot(abs(rho(a[[scorecol]],b[[scorecol]]) - r$spearman)<1e-10,
            abs(length(intersect(ta,tb))/length(union(ta,tb)) - r$top10_jaccard)<1e-10)
}
print(ref[,.(signature,tested_maxRank,spearman,top10_jaccard)])
