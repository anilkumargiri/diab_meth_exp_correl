setwd("/Users/anilkuma/Documents/methylation_paper_T2DM/DNA_methylation_T2DM/Submitted/New_analysis/GSE114065/")
cpg_id=c("cg19693031","cg16171771","cg01178710","cg04673737","cg11024682","cg24458314","cg15815318","cg00574958","cg23372795")
meth=fread("GSE114134_series_matrix.txt", skip=65,select = 1:2)

meth1=meth[meth$ID_REF%in%cpg_id,]

meth1=meth[meth$ID_REF%in%c("cg01178710"),]
RNA=fread("GSE114065_processed_RNAseq.txt")
gene=c("ENSG00000265972","ENSG00000170248","ENSG00000072310","ENSG00000048052","ENSG00000139970","ENSG00000110090","ENSG00000127129","ENSG00000204060","ENSG00000284415","ENSG00000095981","ENSG00000107521","ENSG00000172987","ENSG00000119943","ENSG00000010803","ENSG00000124780","ENSG00000164627","ENSG00000204060")

RNA1=RNA[RNA$V1%in%gene,]

write.csv(RNA1,"RNA1.csv")

cpg=read.csv("meth.csv",header = T,row.names = 1)
rna=read.csv("RNA1.csv",header = T,row.names = 1)
a=intersect(colnames(cpg),colnames(rna))
rna_a=rna[,a]
cpg_a=cpg[,a]

pdf("qqplot_new.pdf")
qq(dat$P,col="blue")
dev.off()

cor.test(as.numeric(rna_a[4,]),as.numeric(cpg_a[5,]))
data123=rbind(rna_a,cpg_a)
write.csv(data123,"comb_data.csv")

pdf("scatterplot_PDCD6IP.pdf")
ggplot(data123, aes(x = ENSG00000170248, y = cg16171771)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()



pdf("scatterplot_TXNIP.pdf")
ggplot(data123, aes(x = ENSG00000265972, y = cg19693031)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()


pdf("scatterplot_HDAC9.pdf")
ggplot(data123, aes(x = ENSG00000048052, y = cg24458314)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()



pdf("scatterplot_SREBF1.pdf")
ggplot(data123, aes(x = ENSG00000072310, y = cg11024682)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()


pdf("scatterplot_CPT1A.pdf")
ggplot(data123, aes(x = ENSG00000110090, y = cg00574958)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()


pdf("scatterplot_SCMH1.pdf")
ggplot(data123, aes(x = ENSG00000010803, y = cg04673737)) +
  geom_point(aes(color = "red")) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1)+theme_classic()
dev.off()



meth_sen=read.table("Methylatioon_sensitivty.txt",header = T,row.names = 1,fill = T)
meth_sen_40=meth_sen[meth_sen$Age>=40,]


meth_sen_40_BMI25=meth_sen_40[meth_sen_40$BMI<25,]
meth_sen_40_BMI25_30=subset(meth_sen_40, BMI>24.98 & BMI<30)
meth_sen_40_BMI30=meth_sen_40[meth_sen_40$BMI>30,]


meth_sen_40_BMI25_40=subset(meth_sen_40_BMI25,Age>39 & Age<50 )
meth_sen_40_BMI25_50=subset(meth_sen_40_BMI25,Age>49 & Age<60 )
meth_sen_40_BMI25_60=subset(meth_sen_40_BMI25,Age>59 & Age<80 )


meth_sen_40_BMI25_30_40=subset(meth_sen_40_BMI25_30,Age>39 & Age<50 )
meth_sen_40_BMI25_30_50=subset(meth_sen_40_BMI25_30,Age>49 & Age<60 )
meth_sen_40_BMI25_30_60=subset(meth_sen_40_BMI25_30,Age>59 & Age<80 )


meth_sen_40_BMI30_40=subset(meth_sen_40_BMI30,Age>39 & Age<50 )
meth_sen_40_BMI30_50=subset(meth_sen_40_BMI30,Age>49 & Age<60 )
meth_sen_40_BMI30_60=subset(meth_sen_40_BMI30,Age>59 & Age<80 )


meth_sen_40_BMI25_40_1=meth_sen_40_BMI25_40[meth_sen_40_BMI25_40$Pheno==1,]
meth_sen_40_BMI25_40_2=meth_sen_40_BMI25_40[meth_sen_40_BMI25_40$Pheno==2,]

a=colMeans(data.matrix(meth_sen_40_BMI25_40_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_40_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI25_40.csv")


meth_sen_40_BMI25_40_1=meth_sen_40_BMI25_40[meth_sen_40_BMI25_40$Pheno==1,]
meth_sen_40_BMI25_40_2=meth_sen_40_BMI25_40[meth_sen_40_BMI25_40$Pheno==2,]

a=colMeans(data.matrix(meth_sen_40_BMI25_40_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_40_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI25_40.csv")


meth_sen_40_BMI25_50_1=meth_sen_40_BMI25_50[meth_sen_40_BMI25_50$Pheno==1,]
meth_sen_40_BMI25_50_2=meth_sen_40_BMI25_50[meth_sen_40_BMI25_50$Pheno==2,]

wilcox.test(meth_sen_40_BMI25_50_1$cg19693031,meth_sen_40_BMI25_50_2$cg19693031)
wilcox.test(meth_sen_40_BMI25_50_1$cg16171771,meth_sen_40_BMI25_50_2$cg16171771)
wilcox.test(meth_sen_40_BMI25_50_1$cg04673737,meth_sen_40_BMI25_50_2$cg04673737)
wilcox.test(meth_sen_40_BMI25_50_1$cg23372795,meth_sen_40_BMI25_50_2$cg23372795)
wilcox.test(meth_sen_40_BMI25_50_1$cg24458314,meth_sen_40_BMI25_50_2$cg24458314)
wilcox.test(meth_sen_40_BMI25_50_1$cg01178710,meth_sen_40_BMI25_50_2$cg01178710)
wilcox.test(meth_sen_40_BMI25_50_1$cg00574958,meth_sen_40_BMI25_50_2$cg00574958)
wilcox.test(meth_sen_40_BMI25_50_1$cg15815318,meth_sen_40_BMI25_50_2$cg15815318)
wilcox.test(meth_sen_40_BMI25_50_1$cg11024682,meth_sen_40_BMI25_50_2$cg11024682)


a=colMeans(data.matrix(meth_sen_40_BMI25_50_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_50_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI25_50.csv")



meth_sen_40_BMI25_60_1=meth_sen_40_BMI25_60[meth_sen_40_BMI25_60$Pheno==1,]
meth_sen_40_BMI25_60_2=meth_sen_40_BMI25_60[meth_sen_40_BMI25_60$Pheno==2,]

a=colMeans(data.matrix(meth_sen_40_BMI25_60_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_60_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)

wilcox.test(meth_sen_40_BMI25_60_1$cg19693031,meth_sen_40_BMI25_60_2$cg19693031)
wilcox.test(meth_sen_40_BMI25_60_1$cg16171771,meth_sen_40_BMI25_60_2$cg16171771)
wilcox.test(meth_sen_40_BMI25_60_1$cg04673737,meth_sen_40_BMI25_60_2$cg04673737)
wilcox.test(meth_sen_40_BMI25_60_1$cg23372795,meth_sen_40_BMI25_60_2$cg23372795)
wilcox.test(meth_sen_40_BMI25_60_1$cg24458314,meth_sen_40_BMI25_60_2$cg24458314)
wilcox.test(meth_sen_40_BMI25_60_1$cg01178710,meth_sen_40_BMI25_60_2$cg01178710)
wilcox.test(meth_sen_40_BMI25_60_1$cg00574958,meth_sen_40_BMI25_60_2$cg00574958)
wilcox.test(meth_sen_40_BMI25_60_1$cg15815318,meth_sen_40_BMI25_60_2$cg15815318)
wilcox.test(meth_sen_40_BMI25_60_1$cg11024682,meth_sen_40_BMI25_60_2$cg11024682)



write.csv(dat1,"colmeans_BMI25_60.csv")



meth_sen_40_BMI25_30_40_1=meth_sen_40_BMI25_30_40[meth_sen_40_BMI25_30_40$Pheno==1,]
meth_sen_40_BMI25_30_40_2=meth_sen_40_BMI25_30_40[meth_sen_40_BMI25_30_40$Pheno==2,]
wilcox.test(meth_sen_40_BMI25_30_40_1$cg19693031,meth_sen_40_BMI25_30_40_2$cg19693031)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg04673737,meth_sen_40_BMI25_30_40_2$cg04673737)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg16171771,meth_sen_40_BMI25_30_40_2$cg16171771)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg23372795,meth_sen_40_BMI25_30_40_2$cg23372795)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg24458314,meth_sen_40_BMI25_30_40_2$cg24458314)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg01178710,meth_sen_40_BMI25_30_40_2$cg01178710)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg00574958,meth_sen_40_BMI25_30_40_2$cg00574958)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg15815318,meth_sen_40_BMI25_30_40_2$cg15815318)
wilcox.test(meth_sen_40_BMI25_30_40_1$cg11024682,meth_sen_40_BMI25_30_40_2$cg11024682)


meth_sen_40_BMI25_30_50_1=meth_sen_40_BMI25_30_50[meth_sen_40_BMI25_30_50$Pheno==1,]
meth_sen_40_BMI25_30_50_2=meth_sen_40_BMI25_30_50[meth_sen_40_BMI25_30_50$Pheno==2,]


wilcox.test(meth_sen_40_BMI25_30_50_1$cg19693031,meth_sen_40_BMI25_30_50_2$cg19693031)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg04673737,meth_sen_40_BMI25_30_50_2$cg04673737)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg16171771,meth_sen_40_BMI25_30_50_2$cg16171771)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg23372795,meth_sen_40_BMI25_30_50_2$cg23372795)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg24458314,meth_sen_40_BMI25_30_50_2$cg24458314)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg01178710,meth_sen_40_BMI25_30_50_2$cg01178710)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg00574958,meth_sen_40_BMI25_30_50_2$cg00574958)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg15815318,meth_sen_40_BMI25_30_50_2$cg15815318)
wilcox.test(meth_sen_40_BMI25_30_50_1$cg11024682,meth_sen_40_BMI25_30_50_2$cg11024682)



a=colMeans(data.matrix(meth_sen_40_BMI25_30_50_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_30_50_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI25_30_50.csv")



meth_sen_40_BMI25_30_60_1=meth_sen_40_BMI25_30_60[meth_sen_40_BMI25_30_60$Pheno==1,]
meth_sen_40_BMI25_30_60_2=meth_sen_40_BMI25_30_60[meth_sen_40_BMI25_30_60$Pheno==2,]

a=colMeans(data.matrix(meth_sen_40_BMI25_30_60_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI25_30_60_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI25_30_60.csv")




meth_sen_40_BMI30_40_1=meth_sen_40_BMI30_40[meth_sen_40_BMI30_40$Pheno==1,]
meth_sen_40_BMI30_40_2=meth_sen_40_BMI30_40[meth_sen_40_BMI30_40$Pheno==2,]

wilcox.test(meth_sen_40_BMI30_40_1$cg19693031,meth_sen_40_BMI30_40_2$cg19693031)
wilcox.test(meth_sen_40_BMI30_40_1$cg04673737,meth_sen_40_BMI30_40_2$cg04673737)
wilcox.test(meth_sen_40_BMI30_40_1$cg16171771,meth_sen_40_BMI30_40_2$cg16171771)
wilcox.test(meth_sen_40_BMI30_40_1$cg23372795,meth_sen_40_BMI30_40_2$cg23372795)
wilcox.test(meth_sen_40_BMI30_40_1$cg24458314,meth_sen_40_BMI30_40_2$cg24458314)
wilcox.test(meth_sen_40_BMI30_40_1$cg01178710,meth_sen_40_BMI30_40_2$cg01178710)
wilcox.test(meth_sen_40_BMI30_40_1$cg00574958,meth_sen_40_BMI30_40_2$cg00574958)
wilcox.test(meth_sen_40_BMI30_40_1$cg15815318,meth_sen_40_BMI30_40_2$cg15815318)
wilcox.test(meth_sen_40_BMI30_40_1$cg11024682,meth_sen_40_BMI30_40_2$cg11024682)




a=colMeans(data.matrix(meth_sen_40_BMI30_40_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI30_40_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI30_40.csv")


meth_sen_40_BMI30_50_1=meth_sen_40_BMI30_50[meth_sen_40_BMI30_50$Pheno==1,]
meth_sen_40_BMI30_50_2=meth_sen_40_BMI30_50[meth_sen_40_BMI30_50$Pheno==2,]

wilcox.test(meth_sen_40_BMI30_50_1$cg19693031,meth_sen_40_BMI30_50_2$cg19693031)
wilcox.test(meth_sen_40_BMI30_50_1$cg04673737,meth_sen_40_BMI30_50_2$cg04673737)
wilcox.test(meth_sen_40_BMI30_50_1$cg16171771,meth_sen_40_BMI30_50_2$cg16171771)
wilcox.test(meth_sen_40_BMI30_50_1$cg23372795,meth_sen_40_BMI30_50_2$cg23372795)
wilcox.test(meth_sen_40_BMI30_50_1$cg24458314,meth_sen_40_BMI30_50_2$cg24458314)
wilcox.test(meth_sen_40_BMI30_50_1$cg01178710,meth_sen_40_BMI30_50_2$cg01178710)
wilcox.test(meth_sen_40_BMI30_50_1$cg00574958,meth_sen_40_BMI30_50_2$cg00574958)
wilcox.test(meth_sen_40_BMI30_50_1$cg15815318,meth_sen_40_BMI30_50_2$cg15815318)
wilcox.test(meth_sen_40_BMI30_50_1$cg11024682,meth_sen_40_BMI30_50_2$cg11024682)




a=colMeans(data.matrix(meth_sen_40_BMI30_50_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI30_50_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI30_50.csv")

meth_sen_40_BMI30_60_1=meth_sen_40_BMI30_60[meth_sen_40_BMI30_60$Pheno==1,]
meth_sen_40_BMI30_60_2=meth_sen_40_BMI30_60[meth_sen_40_BMI30_60$Pheno==2,]

a=colMeans(data.matrix(meth_sen_40_BMI30_60_1),na.rm = T)
b=colMeans(data.matrix(meth_sen_40_BMI30_60_2),na.rm = T)
d=a-b
dat1=data.frame(a,b,d)
write.csv(dat1,"colmeans_BMI30_60.csv")


data=read.csv("all_info_504_withmethpc.csv",header=T,row.names = 1)

betares_batchadjusted_converted=betares_batchadjusted_converted[rownames(data),]
betares_batchadjusted_converted_sele=betares_batchadjusted_converted[,cpg_id]
beta=data.matrix(data[,55:63])
load('BMIQ_normalized_504_qcdone_outlierfixed.Rdata')
beta_illumina_504=beta_illumina512[,rownames(data)]
BMIQ_normalized_504=BMIQ_normalized_504[,rownames(data)]
betares_BMI <- apply(BMIQ_normalized_504,1,function(x){residuals(lm(unlist(x)~data$Age+data$Sex+data$alcohol+data$Smoking+data$CD8T+data$CD4T+data$NK+data$Bcell+data$Mono+data$Gran+data$PC_1+data$PC_2+data$PC_3+data$PC_4+data$PC_5+data$PC1+data$PC2))})##This function returns

betares_BMI <- apply(BMIQ_normalized_504_cpg,1,function(x){summary(lm(unlist(x)~data$Pheno+data$Age+data$Sex+data$alcohol+data$Smoking+data$CD8T+data$CD4T+data$NK+data$Bcell+data$Mono+data$Gran+data$PC_1+data$PC_2+data$PC1+data$PC2))$coefficients})##This function returns


BMIQ_normalized_504_cpg=BMIQ_normalized_504[c("cg19693031","cg11024682","cg00574958","cg08309687"),]
saveRDS(betares,"age_sex_smoking_alcohol_bsconv_cellcomp_chip_lot_vegnonveg_23012023.RDS")
saveRDS(betares_pc1,"age_sex_smoking_alcohol_bsconv_cellcomp_chip_lot_vegnonveg_pc1_23012023.RDS")
saveRDS(betares_BMI,"age_sex_smoking_almmcohol_bsconv_cellcomp_chip_lot_vegnonveg_BMI_23012023.RDS")


betares_cpg=betares[,cpg_id]
betares_BMI=betares_BMI[rownames(data),]

p <- apply(betares_BMI,2,function(x){wilcox.test(unlist(x)~data$Pheno)$p.value})##This function returns

p <- apply(betares_BMI,2,function(x){summary(glm(pheno$Pheno~unlist(x)))$coefficients[2,]})##This function returns

betares_BMI_cpg=betares_BMI[,c("cg19693031","cg16171771","cg01178710","cg04673737","cg11024682","cg24458314","cg15815318","cg00574958","cg23372795")]

test1_cpg=c()
for (i in 1:ncol(betares_cpg) ){
  a=wilcox.test(betares_cpg[1:253,i],betares_cpg[254:504,i] )
  test=t.test(betares_cpg[1:253,i],betares_cpg[254:504,i] )
  test1_cpg=rbind(test1_cpg,c(colnames(betares_cpg)[i], a$p.value,test$p.value,test$estimate))
 print(i)
}

test1_cpg=test1[test1[,1]%in%cpg_id,]
illumina_cpg=beta_illumina512[cpg_id,]


pvalue1 <-as.numeric(p)
chisq1 <- qchisq(1-pvalue1,1)
lambda1 = median(chisq1)/qchisq(0.5,1)
lambda1



pca=prcomp(t(BMIQ_normalized_504))
cpg=fread("cpg_17_rep.txt")
betares_batchadjusted_converted_methpcadjusted_selec_cpg=betares_batchadjusted_converted_methpcadjusted[,cpg_id]

betares_BMI_17cpg=betares_BMI[,rownames(cpg)]
for ( i in 1:17){
  
}
data1=data[rownames(betares_BMI),]
metanew=fread("meta_new1.csv")
signi_cpg=BMIQ_normalized_504[rownames(cpg),]
pheno=fread("all_info_504_withmethpc.csv",header=T)
rownames(pheno)=pheno$CHIP
pheno_control=pheno[pheno$Pheno==2,]
pheno_case=pheno[pheno$Pheno==1,]
signi_cpg_cont=signi_cpg[,rownames(pheno_control)]
signi_cpg_case=signi_cpg[,rownames(pheno_case)]
saveRDS(signi_cpg_case,"signi_cpg_case.RDS")


dataExpr=as.data.frame(t(lnames))









net = blockwiseConsensusModules(
  multiExpr, power = 5, minModuleSize = 10, deepSplit = 2,
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)


powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();


consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]; 


moduleColors
sizeGrWindow(8,6);
pdf(file = "ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()


> table(net$colors)

0  1  2  3 
22 69 19 16

pdf("dendogram_analysis.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                                    "Module colors",
                                dendroLabels = FALSE, hang = 0.03,
                                     addGuide = TRUE, guideHang = 0.05)
dev.off()

nSets = 2
setLabels = c("control", "case")
shortLabels = c("control", "case")
multiExpr = vector(mode = "list", length = nSets)
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(signi_cpg_cont)))
names(multiExpr[[1]]$data) = rownames(signi_cpg_cont)
multiExpr[[2]] = list(data = as.data.frame(t(signi_cpg_case)))
names(multiExpr[[2]]$data) = rownames(signi_cpg_case)
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);


Traits = vector(mode="list", length = nSets);
allTraits=pheno[,c(1,35:55)]
rownames(allTraits)=allTraits$CHIP
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, pheno$CHIP);
  Traits[[set]] = list(data = allTraits[traitRows,-1]);
  a=allTraits[traitRows, ]
  rownames(Traits[[set]]$data) = a$CHIP
}
collectGarbage()



# one step network construction
net = blockwiseModules(datExpr, power = 5,
                       TOMType = "unsigned", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)


plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                                       "Module colors",
                                      dendroLabels = FALSE, hang = 0.03,
                                      addGuide = TRUE, guideHang = 0.05)

colorsFemale = net$colors
multiColor = list(cont = colorsFemale)

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )




ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);






modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();



          medianRank.pres medianRank.qual Zsummary.pres Zsummary.qual
0                 5             5.0         -0.39         -4.30
0.1               4             4.0          6.80          0.78
1                 2             1.5         11.00         10.00
2                 3             3.0          4.70          2.90
3                 1             1.5          4.30          3.70












pheno_control1=pheno_control[,c(36:55)]
moduleTraitCor = cor(MEs, pheno_control1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

pdf("heatmap_cor_trait_eigen.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = names(pheno_control1),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = paste("Module-trait relationships"))
dev.off()



system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )


betares_BMI_newcpg=betares_BMI[rownames(newcpg),]


panc_ip <- FindAllInteractionPrograms(panc_id, iterate.threshold = 3, group.by = "celltype", assay = "alra", sim_threshold = 0.1)