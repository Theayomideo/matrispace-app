### MatriSpace breast cancer case study: deconvoluting TME heterogeneity ###

# Data are from Janesick et al 2023 (https://www.nature.com/articles/s41467-023-43458-x)
# We commonly use these in our tests, so a minimal Surat object for these data is already available and will be loaded at the beginning
# The same data piece is available via MatriSpace online, under the identifier "Ductal carcinoma in situ, Invasive carcinoma (Breast) 1 10x"
# Also, as the Seurat object sports the low-res tissue picture, we will download the original high-res from https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast/Visium Spatial/Spatial Imaging Data

library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(corrplot)
library(scales)
library(UCell)
library(FNN)
library(ggsci)
library(pheatmap)
library(effects)

seu <- readRDS("visium.seurat.breast.rds")
seu <- NormalizeData(seu)
df <- as.data.frame(fread("MatriSpace_Metadata_Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.csv"))
coords <- as.data.frame(seu@images[[1]]@boundaries$centroids@coords)
rownames(coords) <- rownames(seu@meta.data) # you can check that MatriSpace metadata are ordered as the original Seurat metadata with identical(rownames(coords),df$V1)
seu@images[[1]]@scale.factors$spot <- 180

# 1. General tissue overview
SpatialFeaturePlot(seu,"nCount_Spatial") #overall depth
SpatialFeaturePlot(seu,"nFeature_Spatial") #cell complexity
SpatialFeaturePlot(seu,"ERBB2") #general tumor
SpatialFeaturePlot(seu,"SCGB2A2") #DCIS #1
SpatialFeaturePlot(seu,"CPB1") #DCIS #2
SpatialFeaturePlot(seu,"CDH2") #IDC
SpatialFeaturePlot(seu,"SFRP2") #stroma
SpatialFeaturePlot(seu,"KRT17") #myoepithelial
SpatialFeaturePlot(seu,"FABP4") #adipose
SpatialFeaturePlot(seu,"CD4") #immune

# 2. Matrisome scoring
# we start by profiling basement membrane and interstitial niche scores using the same gene lists as in MatriSpace
# expectedly, there is quite an overlap between the two scores, as they tend to crowd the same (microenvironmental) space. However, areas of strong reciprocal expression (above the median for one score and not the other) tend to be much less contiguous than areas with both high or low scores, indicating the specificity of niche signals
g <- c("COLQ","HMCN1","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LAMB1","LAMB2","LAMB3","LAMB4","LAMC1",
       "LAMC2","LAMC3","NID1","NID2","NPNT","NTN1","NTN3","NTN4","NTN5","NTNG1","NTNG2","PAPLN","USH2A",
       "COL15A1","COL18A1","COL28A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6","COL6A1","COL6A2",
       "COL6A3","COL6A5","COL6A6","AGRN","HSPG2") #BM niche genes
g2 <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL5A3","ELN","FBN1","FBN2","FN1","PCOLCE","PCOLCE2",
        "TNC","FMOD","LOX","MMP2","P4HA1","P4HA2","PLOD1","PLOD2","TGM2","TGFB1") #interstitial niche genes
seu <- AddModuleScore_UCell(seu,list(BM_signature=g,I_signature=g2))
SpatialFeaturePlot(seu,"BM_signature_UCell")
SpatialFeaturePlot(seu,"I_signature_UCell")
nscore <- ifelse(seu@meta.data$BM_signature_UCell>median(seu@meta.data$BM_signature_UCell) & seu@meta.data$I_signature_UCell<=median(seu@meta.data$I_signature_UCell),"BM",
                 ifelse(seu@meta.data$BM_signature_UCell<=median(seu@meta.data$BM_signature_UCell) & seu@meta.data$I_signature_UCell>median(seu@meta.data$I_signature_UCell),"I","overlap"))
seu@meta.data$niche <- nscore
ggplot(seu@meta.data,aes(BM_signature_UCell,I_signature_UCell))+
  geom_point(aes(color=niche)) +
  theme_minimal() + xlab("BM signature") + ylab("Interstitial signature") +
  stat_smooth(color="red")
cor.test(seu@meta.data$BM_signature_UCell,seu@meta.data$I_signature_UCell)

knn <- get.knn(coords, k = 6)

cross_contiguity <- function(flagA, flagB) {
  neigh_mat <- knn$nn.index
  
  neighbor_B <- sapply(1:nrow(neigh_mat), function(i) {
    mean(flagB[neigh_mat[i,]])
  })
  
  mean(neighbor_B[flagA])
}

q1 <- seu@meta.data$BM_signature_UCell >= quantile(seu@meta.data$BM_signature_UCell, 0.5) &
  seu@meta.data$I_signature_UCell <= quantile(seu@meta.data$I_signature_UCell, 0.5)
q2 <- seu@meta.data$I_signature_UCell >= quantile(seu@meta.data$I_signature_UCell, 0.5) &
  seu@meta.data$BM_signature_UCell <= quantile(seu@meta.data$BM_signature_UCell, 0.5) 
q3 <- seu@meta.data$I_signature_UCell >= quantile(seu@meta.data$I_signature_UCell, 0.5) &
  seu@meta.data$BM_signature_UCell >= quantile(seu@meta.data$BM_signature_UCell, 0.5)
q4 <- seu@meta.data$I_signature_UCell <= quantile(seu@meta.data$I_signature_UCell, 0.5) &
  seu@meta.data$BM_signature_UCell <= quantile(seu@meta.data$BM_signature_UCell, 0.5)

cross_contiguity(q1,q2)
cross_contiguity(q1,q3)
cross_contiguity(q1,q4)
cross_contiguity(q2,q3)
cross_contiguity(q2,q4)
cross_contiguity(q3,q4)

set.seed(1234)
perm <- replicate(10000, {
  perm_flag <- sample(q2)
  cross_contiguity(q1, perm_flag)
})
mean(perm <= s1) #0

id <- ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"SCGB2A2",]>=2,"DCIS1",
             ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"CPB1",]>=2,"DCIS2",
                    ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"CDH2",]>=2,"IDC",
                           ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"SFRP2",]>=2,"Stroma",
                                  ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"KRT17",]>=2,"Myoepi",
                                         ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"FABP4",]>=2,"Adipose",
                                                ifelse(seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"CD4A",]>=2,"Immune",NA)))))))
             
t1 <- as.data.frame(table(id,q1))
t2 <- as.data.frame(table(id,q2))
t3 <- as.data.frame(table(id,q3))
t4 <- as.data.frame(table(id,q4))
names(t1)[2] <- "category"
names(t2)[2] <- "category"
names(t3)[2] <- "category"
names(t4)[2] <- "category"
t1$signature <- "q1"
t2$signature <- "q2"
t3$signature <- "q3"
t4$signature <- "q4"
fin <- bind_rows(t1,t2,t3,t4)
fin <- fin[fin$category=="TRUE",]
ggplot(fin,aes(id,Freq,fill=signature))+
  geom_bar(stat="identity",position="fill") +
  theme_bw() + xlab("") + ylab("Dominant niche signal (frequency)")

# 3. Matrisome profiles at the border and core of the invasive tumor
# we start by defining a quick function to generate robust scores as in the online version of MatriSpace
rs <- function(obj,
               genes){
  available_genes <- intersect(genes,rownames(obj@assays$Spatial$counts))
  base_counts <- colSums(obj@assays$Spatial$counts[available_genes, , drop = FALSE])
  robust <- (base_counts - median(base_counts)) / IQR(base_counts)
}
# next, we can take advantage of Seurat clusters to define different areas within the IDC lesion. For example, comparing cluster 2 vs. 5
seu@meta.data$clusters <- df$seurat_clusters
IDC_clusters <- df$seurat_clusters
IDC_clusters[seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%"CDH2",]<2] <- NA
seu@meta.data$IDC_clusters <- IDC_clusters
seu@meta.data$IDC_clusters[!seu@meta.data$IDC_clusters%in%c(2,5)] <- NA
SpatialDimPlot(seu,"IDC_clusters") + scale_fill_bmj()

g <- c("F2","F8","F9","F12","F13","KNG","PLG","FGA","FGB","FGG","VTN","VWF","FN1")
g2 <- c("COLQ", "HMCN1", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2",
        "LAMB3", "LAMB4", "LAMC1", "LAMC2", "LAMC3", "NID1", "NID2", "NPNT", "NTN1", "NTN3", "NTN4", "NTN5", 
        "NTNG1", "NTNG2", "PAPLN", "USH2A", "COL15A1", "COL18A1", "COL28A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", 
        "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "AGRN", "HSPG2", "F2", "F8", "F9", "F12", "F13", "KNG", "PLG", "FGA", "FGB", "FGG", 
        "VTN", "VWF", "FN1", "ANGPT2", "VEGFA", "VEGFB"
)
g3 <- matrisome_feature_signatures$`ecm-affiliated_proteins`

hs <- rs(seu,g)
hs2 <- rs(seu,g2)
hs3 <- rs(seu,g3)
seu@meta.data$gs <- hs
seu@meta.data$gs[is.na(seu@meta.data$IDC_clusters)] <- NA
seu@meta.data$gs2 <- hs2
seu@meta.data$gs2[is.na(seu@meta.data$IDC_clusters)] <- NA
seu@meta.data$gs3 <- hs3
seu@meta.data$gs3[is.na(seu@meta.data$IDC_clusters)] <- NA
SpatialFeaturePlot(seu,"gs")
SpatialFeaturePlot(seu,"gs2")
SpatialFeaturePlot(seu,"gs3")
dd <- na.omit(data.frame(v=hs,region=seu@meta.data$IDC_clusters))
dd2 <- na.omit(data.frame(v=hs2,region=seu@meta.data$IDC_clusters))
dd3 <- na.omit(data.frame(v=hs3,region=seu@meta.data$IDC_clusters))
ddf <- bind_rows(dd,dd2,dd3)
ddf$region <- ifelse(ddf$region==5,"border","core")
ddf$region <- factor(ddf$region,levels=c("core","border"))
ddf$signature <- c(rep("hemo",nrow(dd)),rep("pv",nrow(dd2)),rep("ecma",nrow(dd3)))
ggplot(ddf,aes(region,v,fill=region)) +
  geom_bar(stat="identity") +
  scale_fill_bmj() +
  theme_bw() + xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_wrap(~signature,scales = "free",nrow=1)
for(i in unique(ddf$signature)){
  z <- ddf[ddf$signature==i,]
  print(wilcox.test(z$v~z$region))
}

ll <- list()
for(i in g3){
  s <- rs(seu,i)
  d <- data.frame(v1=s,v2=seu@meta.data$IDC_clusters)
  a <- aggregate(d$v1,list(d$v2),mean)
  if(any(is.na(a$x))){
    next
  }else{
    w <- wilcox.test(d$v1~d$v2)
    a <- abs((a[1,2]/a[2,2])) * sign(a[1,2])
    ll[[i]] <- data.frame(gene=i,diff=a,p=w$p.value) 
  }
}
ll <- bind_rows(ll)
ll <- ll[order(-ll$diff),]
ll$p.adj <- p.adjust(ll$p,"BH")
ll <- ll[ll$p.adj<0.05,]

ll2 <- list()
for(i in unique(c(g,g2))){
  s <- rs(seu,i)
  d <- data.frame(v1=s,v2=seu@meta.data$IDC_clusters)
  a <- aggregate(d$v1,list(d$v2),mean)
  if(any(is.na(a$x))){
    next
  }else{
    w <- wilcox.test(d$v1~d$v2)
    a <- abs((a[2,2]/a[1,2])) * sign(a[2,2])
    ll2[[i]] <- data.frame(gene=i,diff=a,p=w$p.value) 
  }
}
ll2 <- bind_rows(ll2)
ll2 <- ll2[order(-ll2$diff),]
ll2$p.adj <- p.adjust(ll2$p,"BH")
ll2 <- ll2[ll2$p.adj<0.05,]

ll$group <- "ECM.affiliated"
ll2$group <- "Vascular"

df_all <- bind_rows(ll,ll2)
df_all$diff[df_all$group=="Vascular"] <- -1*df_all$diff[df_all$group=="Vascular"]

df_all <- df_all %>%
  group_by(gene) %>%
  mutate(mean_diff = mean(diff)) %>%
  ungroup() %>%
  arrange(mean_diff)

df_all$gene <- factor(df_all$gene, levels = unique(df_all$gene))

ggplot(df_all, aes(x = gene, y = diff, group = group, color = group)) +
  geom_segment(aes(xend = gene, y = 0, yend = diff),
               position = position_dodge(width = 0.7),
               size = 3) +
  geom_point(position = position_dodge(width = 0.7), size = 4) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  scale_color_bmj()+
  theme_bw(base_size = 14) +
  guides(x= guide_axis(angle=90)) +
  theme(axis.text.x = element_text(lineheight=0.75))

# given the clear differences in the regulation of these signatures (where practically all genes are one-sided), we suspect differential transcription factor control and test for this hypothesis as done before for MatriCom
tf <- as.data.frame(fread("Curate.txt"))
tf$gene <- gsub(",","",tf$gene)
tf1 <- tf[tf$gene%in%ll$gene,]
tf2 <- tf[tf$gene%in%ll2$gene,]
tf1 <- unique(tf1$TF)
tf2 <- unique(tf2$TF)
int <- intersect(tf1,tf2)
tf1 <- tf1[!(tf1%in%int)]
tf2 <- tf2[!(tf2%in%int)]
tf <- c(tf1,tf2)
k <- seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%tf,]
nm <- seu@meta.data[!is.na(seu$IDC_clusters),]
nm <- rownames(nm)
nm2 <- na.omit(seu@meta.data[seu$IDC_clusters=="2",])
nm2 <- rownames(nm2)
k <- k[,colnames(k)%in%nm]
k <- cor(as.matrix(k))
ann <- data.frame(v=ifelse(rownames(k)%in%nm2,"core","border"))
rownames(ann) <- rownames(k)
bmj_cols <- pal_bmj()(2)
ann$v <- factor(ann$v,levels=c("core","border"))
ann_colors <- list(
    v = setNames(
      bmj_cols[1:length(levels(ann$v))],
      levels(ann$v)
    )
)
pheatmap(k,
         annotation_row = ann,
         annotation_colors = ann_colors)

ll3 <- list()
for(i in tf){
  s <- rs(seu,i)
  d <- data.frame(v1=s,v2=seu@meta.data$IDC_clusters)
  a <- aggregate(d$v1,list(d$v2),mean)
  if(any(is.na(a$x))){
    next
  }else{
    w <- wilcox.test(d$v1~d$v2)
    if(a[1,2]>a[2,2]){
      a <- abs((a[1,2]/a[2,2])) * sign(a[1,2])
    }else{
      a <- abs((a[1,2]/a[2,2])) * sign(a[2,2])
    }
    ll3[[i]] <- data.frame(gene=i,diff=a,p=w$p.value) 
  }
}
ll3 <- bind_rows(ll3)
ll3 <- ll3[order(-ll3$diff),]
ll3$p.adj <- p.adjust(ll3$p,"BH")
ll3 <- ll3[ll3$p.adj<0.05,]
k <- seu@assays$Spatial$data[rownames(seu@assays$Spatial$data)%in%tf,]
k <- k[,colnames(k)%in%nm]
ex1 <- data.frame(v=k[rownames(k)%in%"MYCN"],a=ann$v)
ex2 <- data.frame(v=k[rownames(k)%in%"IKZF1"],a=ann$v)

ggplot(ex1, aes(x = v, fill = a, color = a)) +
  geom_density(alpha = 0.3, size = 1) +
  scale_fill_bmj() +
  scale_color_bmj() +
  theme_bw() +
  labs(
    x = "Value (v)",
    y = "Density",
    fill = "Region",
    color = "Region"
  )

wilcox.test(v ~ a, data = ex1)

ggplot(ex2, aes(x = v, fill = a, color = a)) +
  geom_density(alpha = 0.3, size = 1) +
  scale_fill_bmj() +
  scale_color_bmj() +
  theme_bw() +
  labs(
    x = "Value (v)",
    y = "Density",
    fill = "Region",
    color = "Region"
  )

wilcox.test(v ~ a, data = ex2)

# 4. Selected matrisome features with strong locoregional clustering
# the original stduy from Janesick et al. shows an intra-adipose small triple-positive tumor lesion (DCIS #2-type) different from the rest of the tumor (double-positive)
# we start by marking this area
SpatialFeaturePlot(seu,"ERBB2")
SpatialFeaturePlot(seu,"ESR1")
SpatialFeaturePlot(seu,"PGR")

# then we notice that it's clearly less proliferative than rest of the tumor, in an area of high DDR1 and a clear collagen switch
SpatialFeaturePlot(seu,"MKI67")
SpatialFeaturePlot(seu,"DDR1")
SpatialFeaturePlot(seu,"COL1A1")
SpatialFeaturePlot(seu,"COL3A1")
s1 <- rs(seu,c("COL3A1","DDR1"))
s2 <- rs(seu,c("COL1A1","DDR1"))
s1[s1<0] <- 0
s2[s2<0] <- 0
s3 <- s1-s2
s3[s3<0] <- 0
seu@meta.data$ss <- s3
SpatialFeaturePlot(seu,"ss")

# is this a dormant niche supported by COLIII/DDR1?
df <- data.frame(
  MKI67 = seu@assays$Spatial$data["MKI67",],
  s1 = s1,
  s2 = s2)

fit <- glm(MKI67 ~ s1 + s2, data = df)
plot(allEffects(fit),color=c("blue","red"))

eff_s1 <- as.data.frame(Effect("s1", fit))
eff_s2 <- as.data.frame(Effect("s2", fit))
colnames(eff_s1)[colnames(eff_s1) == "s1"] <- "predictor"
colnames(eff_s2)[colnames(eff_s2) == "s2"] <- "predictor"
coefs <- summary(fit)$coefficients
p_s1 <- signif(coefs["s1", "Pr(>|t|)"], 3)
p_s2 <- signif(coefs["s2", "Pr(>|t|)"], 3)
eff_s1$type <- "COL3A1–DDR1 (s1)"
eff_s2$type <- "COL1A1–DDR1 (s2)"
eff_all <- rbind(eff_s1, eff_s2)
ggplot(eff_all, aes(x = predictor, y = fit, color = type, fill = type)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  theme_bw(base_size = 14) +
  labs(
    x = "ECM–DDR1 score",
    y = "Predicted MKI67",
    color = "",
    fill = ""
  ) +
  scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
  scale_fill_manual(values = c("#2C7BB6", "#D7191C"))

   # + annotate("text", x = Inf, y = Inf,
   #         label = paste0("s1 p = ", p_s1, "\n",
   #                        "s2 p = ", p_s2),
   #         hjust = 2.75, vjust = 1.5,
   #         size = 4)
