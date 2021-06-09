#Read and subset the clinical data
clin <- read.table("TCGA_PCPG_clinical_molSubtype.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
for(i in 1:nrow(clin)){
alls <- unlist(strsplit(clin[i,1],"-"))
clin[i,1] <- paste(alls[1],alls[2],alls[3],alls[4],sep=".")
}
clinpseudo <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia")
clincortical <- subset(clin, mRNA.Subtype.Clusters=="Cortical admixture")
clinwnt <- subset(clin, mRNA.Subtype.Clusters=="Wnt-altered")
clinkinase <- subset(clin, mRNA.Subtype.Clusters=="Kinase signaling")
clinsdhb <- subset(clin, !is.na(SDHB.Germline.Mutation))
clinsdhd <- subset(clin, !is.na(SDHD.Germline.Mutation))
clinsdhx <- unique(c(clinsdhb[,1],clinsdhd[,1]))
clinaggressive <- subset(clin, Clinically.Aggressive.and.or.Metastatic == "Yes" | Clinically.Aggressive.and.or.Metastatic == "yes")
clinnonaggressive <- subset(clin, Clinically.Aggressive.and.or.Metastatic == "no" | Clinically.Aggressive.and.or.Metastatic == "No")

clin <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia" | mRNA.Subtype.Clusters=="Cortical admixture" | mRNA.Subtype.Clusters=="Wnt-altered" | mRNA.Subtype.Clusters=="Kinase signaling" | (!is.na(SDHB.Germline.Mutation)) | (!is.na(SDHD.Germline.Mutation)), select=c("SampleID","Clinically.Aggressive.and.or.Metastatic","SDHB.Germline.Mutation","ATRX.Somatic.Mutation","Tumor.Location"))

#Process the miRNA-Seq data from GDC PCPG cohort 

pcpg <- read.table("PCPG_GDC_miRCountData.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pcpg <- na.omit(pcpg)
rownames(pcpg) <- pcpg[,1]
pcpg <- pcpg[,-1]
clinpcpg <- colnames(pcpg)
spseudo <- clinpseudo[,1]
swnt  <- clinwnt[,1]
skinase  <- clinkinase[,1]
scortical <- clincortical[,1]
ssdhb <- clinsdhb[,1]
ssdhd <- clinsdhd[,1]
ssdhx <- unique(c(ssdhb,ssdhd))
#ssdhx <- ssdhb
snormal <- c("TCGA.P8.A5KD.11A","TCGA.SQ.A6I4.11A","TCGA.P8.A5KC.11A")
spseudo <- setdiff(spseudo,ssdhx)
#pcpglinc <- subset(pcpg, rownames(pcpg) %in% c(linc[,1],"TERT"))
pcpglinc <- pcpg
normalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% snormal)
sdhxexpr <- subset(pcpglinc, select=names(pcpglinc) %in% ssdhx)
nonsdhxexpr <- subset(pcpglinc, select=!(names(pcpglinc) %in% ssdhx))
pseudoexpr <- subset(pcpglinc, select=names(pcpglinc) %in% spseudo)
nonpseudoexpr <- subset(pcpglinc, select=!(names(pcpglinc) %in% spseudo))
wntexpr <- subset(pcpglinc, select=names(pcpglinc) %in% swnt)
kinaseexpr <- subset(pcpglinc, select=names(pcpglinc) %in% skinase)
corticalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% scortical)
dataexp <- cbind(sdhxexpr,pseudoexpr, wntexpr, kinaseexpr, corticalexpr, normalexpr)
#dataexp <- cbind(pseudoexpr,nonpseudoexpr)
ms <- numeric(ncol(sdhxexpr))
#nms <- numeric(ncol(nonpseudoexpr))
nms <- numeric(ncol(pseudoexpr))
ws <- numeric(ncol(wntexpr))
ks <- numeric(ncol(kinaseexpr))
cs <- numeric(ncol(corticalexpr))
ns <- numeric(ncol(normalexpr))
for(i in 1:ncol(sdhxexpr))
ms[i] <- 1
#for(i in 1:ncol(nonpseudoexpr))
#nms[i] <- 2
for(i in 1:ncol(pseudoexpr))
nms[i] <- 2
for(i in 1:ncol(wntexpr))
ws[i] <- 3
for(i in 1:ncol(kinaseexpr))
ks[i] <- 4
for(i in 1:ncol(corticalexpr))
cs[i] <- 5
for(i in 1:ncol(normalexpr))
ns[i] <- 6
group <- factor(c(ms,nms,ws,ks,cs,ns))
#group <- factor(c(ms,nms))
y <- edgeR::DGEList(counts=dataexp,group=group)
y <- edgeR::calcNormFactors(y, na.rm=TRUE)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
logcpm <- edgeR::cpm(y, prior.count=2, log=TRUE)
dataexp <- t(scale(t(logcpm)))
dataexpmir <- dataexp

dataexpmir <- t(scale(t(logcpm)))
aggrexpr <- subset(dataexpmir, select=colnames(dataexpmir) %in% clinaggressive[,1])
nonaggrexpr <- subset(dataexpmir, select=colnames(dataexpmir) %in% clinnonaggressive[,1])
dataexpaggr <- cbind(aggrexpr,nonaggrexpr)

#Cluster the miRNAs differentially expressed in SDHx and non-SDHx mutated pseudohypoxia subtype
mirsdhxup <- read.table("PCPG_SDHX_mirexpUp.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#mirsdhxup <- subset(mirsdhxup, logFC < -2)
mirpseudoup <- read.table("PCPG_nonSDHX-pseudo_mirexpUp.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#mirpseudoup <- subset(mirpseudoup, logFC < -2)
mirsdhxdown <- read.table("PCPG_SDHX_mirexpDown.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#mirsdhxdown <- subset(mirsdhxdown, logFC > 2)
mirpseudodown <- read.table("PCPG_nonSDHX-pseudo_mirexpDown.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#mirpseudodown <- subset(mirpseudodown, logFC > 2)
dfmirpseudo <- unique(c(mirsdhxup[,1],mirpseudoup[,1],mirsdhxdown[,1],mirpseudodown[,1]))
dataexp <- t(scale(t(dataexp)))
dataexp <- subset(dataexp, rownames(dataexp) %in% dfmirpseudo)
corpath <- cor(t(dataexp), method="spearman")
#corpath[is.na(corpath)] <- 0
d <- as.dist(1-corpath)
corpath2 <- cor(dataexp, method="spearman")
#corpath2[is.na(corpath2)] <- 0
d2 <- as.dist(1-corpath2)
hc <- hclust(d2, method = "complete", members = NULL)
hr <- hclust(d, method = "complete", members = NULL)
clustrow <- cutree(hr, k=4)
clustrow[clustrow==1] <- "#00FFFF"
clustrow[clustrow==2] <- "#E41A1C"
clustrow[clustrow==3] <- "#865b98"
clustrow[clustrow==4] <- "#FFDEAD"
png("miRClusters_tcga.png", width=6400, height=1200, res=300)
plot(hr, cex=0.8)
dev.off()
hmcols<-colorRampPalette(c("midnightblue","blue4","mediumblue","blue","dodgerblue","white","red","red1","red2","red3","darkred"))(256)
#hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
clustcol.height = c(rep("#00FFFF",ncol(sdhxexpr)),rep("#25A309",ncol(pseudoexpr)),rep("#E2EF28",ncol(wntexpr)), rep("#865b98",ncol(kinaseexpr)),rep("#FF7256",ncol(corticalexpr)),rep("#000000",ncol(normalexpr)))
#colab <- c(rep("Pseudohypoxia",ncol(pseudoexpr)),rep("WntAltered",ncol(wntexpr)),rep("KinaseSignaling",ncol(kinaseexpr)),rep("CorticalAdmixture",ncol(corticalexpr)))
png("miRExpression_ExtendedMolecular_Subtypes_clusters.png", width=3400, height=9200, res=300)
par(mai=c(5,8,2,5), pin=c(3,3))
gplots::heatmap.2(as.matrix(dataexp), scale='none', Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=hmcols, dendrogram="both", cexRow=1.4, trace="none", offsetRow=0.001, ColSideColors=clustcol.height, RowSideColors=clustrow, labCol=FALSE, margins=c(5,12),lwid=c(0.5,4),lhei=c(0.3,4))
dev.off()


#Process the RNA-Seq data from GDC PCPG cohort HTseq count pipeline

pcpg <- read.table("PCPG_GDC_CountData_all_nr.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pcpg <- na.omit(pcpg)
rownames(pcpg) <- pcpg[,1]
pcpg <- pcpg[,-1]
clinpcpg <- colnames(pcpg)
spseudo <- clinpseudo[,1]
swnt  <- clinwnt[,1]
skinase  <- clinkinase[,1]
scortical <- clincortical[,1]
ssdhb <- clinsdhb[,1]
ssdhd <- clinsdhd[,1]
ssdhx <- unique(c(ssdhb,ssdhd))
#ssdhx <- ssdhb
snormal <- c("TCGA.P8.A5KD.11A","TCGA.SQ.A6I4.11A","TCGA.P8.A5KC.11A")
spseudo <- setdiff(spseudo,ssdhx)
#pcpglinc <- subset(pcpg, rownames(pcpg) %in% c(linc[,1],"TERT"))
pcpglinc <- pcpg
normalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% snormal)
sdhxexpr <- subset(pcpglinc, select=names(pcpglinc) %in% ssdhx)
pseudoexpr <- subset(pcpglinc, select=names(pcpglinc) %in% spseudo)
wntexpr <- subset(pcpglinc, select=names(pcpglinc) %in% swnt)
kinaseexpr <- subset(pcpglinc, select=names(pcpglinc) %in% skinase)
corticalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% scortical)
dataexp <- cbind(sdhxexpr,pseudoexpr, wntexpr, kinaseexpr, corticalexpr, normalexpr)
ms <- numeric(ncol(sdhxexpr))
nms <- numeric(ncol(pseudoexpr))
ws <- numeric(ncol(wntexpr))
ks <- numeric(ncol(kinaseexpr))
cs <- numeric(ncol(corticalexpr))
#ns <- numeric(ncol(normalexpr))
for(i in 1:ncol(sdhxexpr))
ms[i] <- 1
for(i in 1:ncol(pseudoexpr))
nms[i] <- 2
for(i in 1:ncol(wntexpr))
ws[i] <- 3
for(i in 1:ncol(kinaseexpr))
ks[i] <- 4
for(i in 1:ncol(corticalexpr))
cs[i] <- 5
#for(i in 1:ncol(normalexpr))
#ns[i] <- 6
group <- factor(c(ms,nms,ws,ks,cs))
y <- edgeR::DGEList(counts=dataexp,group=group)
y <- edgeR::calcNormFactors(y, na.rm=TRUE)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
logcpm <- edgeR::cpm(y, prior.count=2, log=TRUE)
dataexp <- t(scale(t(logcpm)))
aggrexpr <- subset(dataexp, select=colnames(dataexp) %in% clinaggressive[,1])
nonaggrexpr <- subset(dataexp, select=colnames(dataexp) %in% clinnonaggressive[,1])
dataexpaggr <- cbind(aggrexpr,nonaggrexpr)
datamirtarget <- logcpm


#Get the target lncRNAs and mRNAs for miRNAs differentially expressed in SDHx and non-SDHx mutated pseudohypoxia
mirtargetlncsdhx <- read.table("SDHx-miR-lncRNA_targets.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
datamirtarget <- subset(dataexp, rownames(dataexp) %in% unique(mirtargetlncsdhx[,2]))
sdhxmirex <- c("hsa-mir-182","hsa-mir-183","hsa-mir-96","hsa-mir-490","hsa-mir-383","hsa-mir-891a")
pseudomirex <- c(sdhxmirex,"hsa-mir-6715a","hsa-mir-551a","hsa-mir-935","hsa-mir-7151","hsa-mir-1236","hsa-mir-3193","hsa-mir-582","hsa-mir-504","hsa-mir-582","hsa-mir-34c","hsa-mir-34b","hsa-mir-592","hsa-mir-208b","hsa-mir-642a","hsa-mir-6716","hsa-mir-3912","hsa-mir-670","hsa-mir-3929","hsa-mir-483","hsa-mir-210","hsa-mir-6844","hsa-mir-4691","hsa-mir-4742")
dataexpmir1 <- subset(dataexpmir, rownames(dataexpmir) %in% pseudomirex, select=colnames(datamirtarget))

#calculate expression correlation of the miRNAs and target lncRNAs
mirsdtargetcor <- data.frame(microRNA="m", Gene_name="g", expression.correlation=0, correlation.pvalue=1)
for(k in 1:nrow(dataexpmir1)){
data1 <- as.numeric(dataexpmir1[k,])
slcor <- numeric(nrow(datamirtarget))
slcorp <- numeric(nrow(datamirtarget))
#miR target correlation analysis
for(i in 1:nrow(datamirtarget)){
data2 <- datamirtarget[i,]
if(length(data2)>2){
data2 <- as.numeric(data2)
if(length(data1)>=4 && length(data2)>=4){
genemat <- as.matrix(c(data1,data2),nrow=length(data1),ncol=2, bycol=T)
cv <- cor(data1,data2,method="pearson")
if(!is.na(cv)){
slcor[i] <- as.numeric(cor.test(data1,data2,method="pearson")$estimate)
slcorp[i] <- cor.test(data1,data2,method="pearson")$p.value
}
}
else{
slcor[i] <- 0
slcorp[i] <- 0
}
}
}
mirsdtargetcor <- rbind(mirsdtargetcor, data.frame(microRNA=rep(rownames(dataexpmir1)[k], nrow(datamirtarget)), Gene_name=rownames(datamirtarget), expression.correlation=slcor, correlation.pvalue=slcorp))
}
mirsdtargetcor <- mirsdtargetcor[-1,]
mirsdtargetcorsig <- subset(mirsdtargetcor, correlation.pvalue<0.01 & expression.correlation< -0.33)


#Select the miRNAs and target lncRNAs related to metastatic status using general linear model
mirsdtargetcorsig <- read.table("SDHx-miR-lncRNA_corsig.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
lncsd <- unique(as.character(subset(mirsdtargetcorsig, microRNA %in% c("hsa-mir-182","hsa-mir-183","hsa-mir-383","hsa-mir-96"))[,2]))
mirsdtargetcorsig[,1] <- as.character(mirsdtargetcorsig[,1])
mirsdtargetcorsig[,2] <- as.character(mirsdtargetcorsig[,2])
dataexpmirsig <- subset(dataexpmir, rownames(dataexpmir) %in% mirsdtargetcorsig[,1], select=colnames(datamirtarget))
dataexpmirsig <- data.frame(sampleID=colnames(dataexpmirsig), t(dataexpmirsig))
dataexptargetsig <- subset(datamirtarget, rownames(datamirtarget) %in% mirsdtargetcorsig[,2])
dataexptargetsig <- data.frame(sampleID=colnames(dataexptargetsig), t(dataexptargetsig))
dataexpmirtargetsig <- dplyr::inner_join(dataexpmirsig,dataexptargetsig)
dataexpmirsig <- subset(dataexpmirtargetsig, select=colnames(dataexpmirtargetsig) %in% colnames(dataexpmirsig))
dataexptargetsig <- subset(dataexpmirtargetsig, select=colnames(dataexpmirtargetsig) %in% colnames(dataexptargetsig))
datamat <- dataexpmirtargetsig
colnames(datamat)[1] <- "SampleID"
clindata <- subset(clin, select=c(SampleID, Clinically.Aggressive.and.or.Metastatic))
clindataexp <- dplyr::inner_join(clindata,datamat)
rownames(clindataexp) <- clindataexp[,1]
clindataexp <- clindataexp[,-1]
clindataexp$Clinically.Aggressive.and.or.Metastatic[clindataexp$Clinically.Aggressive.and.or.Metastatic=="Yes" | clindataexp$Clinically.Aggressive.and.or.Metastatic=="yes"] <- 1
clindataexp$Clinically.Aggressive.and.or.Metastatic[clindataexp$Clinically.Aggressive.and.or.Metastatic=="No" | clindataexp$Clinically.Aggressive.and.or.Metastatic=="no"] <- 0
clindataexp$Clinically.Aggressive.and.or.Metastatic <- as.numeric(clindataexp$Clinically.Aggressive.and.or.Metastatic)
#GLM
f1 <- as.formula(paste("Clinically.Aggressive.and.or.Metastatic ~", paste(colnames(clindataexp)[c(2:ncol(clindataexp))], collapse= "+")))
glmtargetset <- glm(f1, data=clindataexp, family=gaussian())
dfglm <- as.data.frame(summary(glmtargetset)[[13]])
colnames(dfglm)[4] <- "Pr"
selectmirtarset <- subset(dfglm, Pr < 0.1)
for(i in 1:length(lncsd))
lncsd[i] <- gsub("-", ".", lncsd[i])
selectmirtarset <- subset(dfglm, Pr < 0.1 & rownames(dfglm) %in% c("hsa.mir.182","hsa.mir.183","hsa.mir.383",lncsd))


#Perform multivariate analysis with relevant clinical parameters, miRNAs and lncRNAs to find candidates associated with metastasis-free survival
dataexpmir1 <- subset(dataexpmir, rownames(dataexpmir) %in% c("hsa-mir-182","hsa-mir-183", "hsa-mir-96", "hsa-mir-383"))
dataexpmir1 <- data.frame(sampleId=colnames(dataexpmir1), t(dataexpmir1))

selectar <- read.table("mir_lnc_target_GLM.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dataexp1 <- subset(dataexp, rownames(dataexp) %in% c("LINC00877","MAGI2-AS3","USP3-AS1","AC009312.1"))
dataexptarget1 <- data.frame(sampleId=colnames(dataexp1), t(dataexp1))

dataexp2 <- subset(dataexp, rownames(dataexp) %in% c("TERT","HIF1A-AS2"))
dataexptarget2 <- data.frame(sampleId=colnames(dataexp2), t(dataexp2))

dataexpmirtarget1 <- dplyr::inner_join(dataexpmir1,dataexptarget1)
dataexpmirtarget1 <- dplyr::inner_join(dataexpmirtarget1,dataexptarget2)

clin <- read.table("TCGA_PCPG_clinical_molSubtype.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
for(i in 1:nrow(clin)){
alls <- unlist(strsplit(clin[i,1],"-"))
clin[i,1] <- paste(alls[1],alls[2],alls[3],alls[4],sep=".")
}
clin <- within(clin, SDHB.Germline.Mutation[!is.na(SDHB.Germline.Mutation)] <- 1.0)
clin <- within(clin, SDHB.Germline.Mutation[is.na(SDHB.Germline.Mutation)] <- 0.0)
clin <- within(clin, ATRX.Somatic.Mutation[ATRX.Somatic.Mutation != 0] <- 1.0)
clin <- within(clin, ATRX.Somatic.Mutation[ATRX.Somatic.Mutation == 0] <- 0.0)
clin <- within(clin, Tumor.Location[Tumor.Location != "adrenal"] <- 1.0)
clin <- within(clin, Tumor.Location[Tumor.Location == "adrenal"] <- 0.0)
clin <- within(clin, Normetanephrine.secreting[Normetanephrine.secreting == "Yes"] <- 1.0)
clin <- within(clin, Normetanephrine.secreting[Normetanephrine.secreting != 1.0] <- 0.0)
clin <- within(clin, Norepinephrine.Secreting[Norepinephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Norepinephrine.Secreting[Norepinephrine.Secreting != 1.0] <- 0.0)
clin <- within(clin, Epinephrine.Secreting[Epinephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Epinephrine.Secreting[Epinephrine.Secreting != 1.0] <- 0.0)
clin <- within(clin, Metanephrine.Secreting[Metanephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Metanephrine.Secreting[Metanephrine.Secreting != 1.0] <- 0.0)
clin <- within(clin, Methoxytyramine.Secreting[Methoxytyramine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Methoxytyramine.Secreting[Methoxytyramine.Secreting != 1.0] <- 0.0)
clin <- within(clin, Dopamine.Secreting[Dopamine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Dopamine.Secreting[Dopamine.Secreting != 1.0] <- 0.0)

clin$SDHB.Germline.Mutation <- as.numeric(clin$SDHB.Germline.Mutation)
clin$ATRX.Somatic.Mutation <- as.numeric(clin$ATRX.Somatic.Mutation)
clin$Tumor.Location <- as.numeric(clin$Tumor.Location)
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting)
clin$Norepinephrine.Secreting <- as.numeric(clin$Norepinephrine.Secreting)
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting)
clin$Metanephrine.Secreting <- as.numeric(clin$Metanephrine.Secreting)
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting)
clin$Dopamine.Secreting <- as.numeric(clin$Dopamine.Secreting)
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting) + as.numeric(clin$Norepinephrine.Secreting)
clin$Normetanephrine.secreting[clin$Normetanephrine.secreting > 1.0] <- 1.0
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting)
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting) + as.numeric(clin$Metanephrine.Secreting)
clin$Epinephrine.Secreting[clin$Epinephrine.Secreting > 1.0] <- 1.0
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting)
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting) + as.numeric(clin$Dopamine.Secreting)
clin$Methoxytyramine.Secreting[clin$Methoxytyramine.Secreting > 1.0] <- 1.0
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting)
clinpseudo <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia")
clincortical <- subset(clin, mRNA.Subtype.Clusters=="Cortical admixture")
clinwnt <- subset(clin, mRNA.Subtype.Clusters=="Wnt-altered")
clinkinase <- subset(clin, mRNA.Subtype.Clusters=="Kinase signaling")
clinsdhb <- subset(clin, !is.na(SDHB.Germline.Mutation))
clinsdhd <- subset(clin, !is.na(SDHD.Germline.Mutation))
clinret <- subset(clin, !is.na(RET.Germline.Mutation))
clinretsom <- subset(clin, RET.Somatic.Mutation!=0)
clinret <- unique(rbind(clinret,clinretsom))
clinnf1 <- subset(clin, !is.na(NF1.Germline.Mutation))
clinnf1som <- subset(clin, NF1.Somatic.Mutation!=0)
clinnf1 <- unique(rbind(clinnf1,clinnf1som))
clinvhl <- subset(clin, !is.na(VHL.Germline.Mutation))
clinvhlsom <- subset(clin, VHL.Somatic.Mutation!=0)
clinvhl <- unique(rbind(clinvhl,clinvhlsom))
clinepas1 <- subset(clin, EPAS1.Somatic.Mutation!=0)
clinhras <- subset(clin, HRAS.Somatic.Mutation!=0)
colnames(clin)[26] <- "DFS_MONTHS"
clin[,26] <- as.numeric(clin[,26])/30
colnames(clin)[16] <- "DFS_STATUS"
clin <- within(clin, DFS_STATUS[DFS_STATUS == 'yes' | DFS_STATUS == 'Yes'] <- 'Recurred/Progressed')
clin <- within(clin, DFS_STATUS[DFS_STATUS == 'no' | DFS_STATUS == 'No'] <- 'DiseaseFree')
colnames(clin)[6] <- "AGE"
rownames(clin) <- clin[,1]
clin <- clin[,-1]
datac <- data.frame(clin,sampleId=rownames(clin))
clinical <- subset(datac, select=c(DFS_STATUS,DFS_MONTHS,AGE,sampleId,SDHB.Germline.Mutation,ATRX.Somatic.Mutation,Tumor.Location,Normetanephrine.secreting,Epinephrine.Secreting,Methoxytyramine.Secreting))
datametasig <- dataexpmirtarget1
clinical <- dplyr::inner_join(clinical,datametasig)
rownames(clinical) <- clinical$sampleId
clinical <- clinical[,-4]
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- '0')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == 'NA'] <- '0')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
f1 <- as.formula(paste("survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~", paste(colnames(clinical)[4:ncol(clinical)], collapse= "+")))
cc <- survival::coxph.control(eps = 1e-40, toler.chol = .Machine$double.eps^0.55, iter.max = 10000000, toler.inf = 1e-5, outer.max = 10000)
coxdiff <- survival::coxph(f1, data=clinical, control=cc)
