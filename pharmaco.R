library(ggplot2)
library(ROCR)

normalized<-function(y) {
x<-y[!is.na(y)]
x<-(x - min(x)) / (max(x) - min(x))
y[!is.na(y)]<-x
return(y) }

#Read in ROC data files
ROC <- read.csv2("ROC_data.csv",sep="\t",header=TRUE)
ROC34 <- read.csv2("ROC_data_cat34.csv",sep="\t",header=TRUE)

#Polyphen first
poly_df <- data.frame(ROC$PolyPhen.Score,ROC$Labels)
poly_df <- poly_df[!is.na(poly_df$ROC.PolyPhen.Score),]
poly_pred <- prediction(as.numeric(as.character(poly_df$ROC.PolyPhen.Score)),poly_df$ROC.Labels)
poly_perf <- performance(poly_pred,"tpr","fpr")
plot(poly_perf,main="PolyPhen ROC Plot Type A Enriched Variants")

#AUC
poly_auc <- performance(poly_pred, measure = "auc")
poly_auc@y.values[[1]]

#MCC
poly_mcc <- performance(poly_pred,"mat")
poly_mcc@y.values[[1]][which.max(poly_mcc@y.values[[1]])]

poly_df <- data.frame(ROC34$PolyPhen.Score,ROC34$Labels)
poly_df <- poly_df[!is.na(poly_df$ROC34.PolyPhen.Score),]
poly_pred <- prediction(as.numeric(as.character(poly_df$ROC34.PolyPhen.Score)),poly_df$ROC34.Labels)
poly_perf <- performance(poly_pred,"tpr","fpr")
plot(poly_perf,main="PolyPhen ROC Plot Type B Enriched Variants")

#AUC
poly_auc <- performance(poly_pred, measure = "auc")
poly_auc@y.values[[1]]

#MCC
poly_mcc <- performance(poly_pred,"mat")
poly_mcc@y.values[[1]][which.max(poly_mcc@y.values[[1]])]

#Mutation predictor
mutpred_df <- data.frame(ROC$MutPred.Score,ROC$Labels)
mutpred_df <- mutpred_df[!is.na(mutpred_df$ROC.MutPred.Score),]
mutpred_pred <- prediction(as.numeric(as.character(mutpred_df$ROC.MutPred.Score)),mutpred_df$ROC.Labels)
mutpred_perf <- performance(mutpred_pred,"tpr","fpr")
plot(mutpred_perf,main="MutPred ROC Plot Type A Enriched Variants")

#AUC
mutpred_auc <- performance(mutpred_pred, measure = "auc")
mutpred_auc@y.values[[1]]

#MCC
mutpred_mcc <- performance(mutpred_pred,"mat")
mutpred_mcc@y.values[[1]][which.max(mutpred_mcc@y.values[[1]])]

mutpred_df <- data.frame(ROC34$MutPred.Score,ROC34$Labels)
mutpred_df <- mutpred_df[!is.na(mutpred_df$ROC34.MutPred.Score),]
mutpred_pred <- prediction(as.numeric(as.character(mutpred_df$ROC34.MutPred.Score)),mutpred_df$ROC34.Labels)
mutpred_perf <- performance(mutpred_pred,"tpr","fpr")
plot(mutpred_perf,main="MutPred ROC Plot Type B Enriched Variants")

#AUC
mutpred_auc <- performance(mutpred_pred, measure = "auc")
mutpred_auc@y.values[[1]]

#MCC
mutpred_mcc <- performance(mutpred_pred,"mat")
mutpred_mcc@y.values[[1]][which.max(mutpred_mcc@y.values[[1]])]


revel_df <- data.frame(ROC$REVEL.Score,ROC$Labels)
revel_df <- revel_df[!is.na(revel_df$ROC.REVEL.Score),]
revel_pred <- prediction(as.numeric(as.character(revel_df$ROC.REVEL.Score)),revel_df$ROC.Labels)
revel_perf <- performance(revel_pred,"tpr","fpr")
plot(revel_perf,main="REVEL ROC Plot Type A Enriched Variants")

#AUC
revel_auc <- performance(revel_pred, measure = "auc")
revel_auc@y.values[[1]]

#MCC
revel_mcc <- performance(revel_pred,"mat")
revel_mcc@y.values[[1]][which.max(revel_mcc@y.values[[1]])]


#REVEL
revel_df <- data.frame(ROC34$REVEL.Score,ROC34$Labels)
revel_df <- revel_df[!is.na(revel_df$ROC34.REVEL.Score),]
revel_pred <- prediction(as.numeric(as.character(revel_df$ROC34.REVEL.Score)),revel_df$ROC34.Labels)
revel_perf <- performance(revel_pred,"tpr","fpr")
plot(revel_perf,main="REVEL ROC Plot Type B Enriched Variants")

#AUC
revel_auc <- performance(revel_pred, measure = "auc")
revel_auc@y.values[[1]]

#MCC
revel_mcc <- performance(revel_pred,"mat")
revel_mcc@y.values[[1]][which.max(revel_mcc@y.values[[1]])]

#SIFT (need to flip significance of labels; 0 is most damaging)
sift_df <- data.frame(ROC$SIFT.Score,ROC$LabelsFlip)
sift_df <- sift_df[!is.na(sift_df$ROC.SIFT.Score),]
sift_pred <- prediction(as.numeric(as.character(sift_df$ROC.SIFT.Score)),sift_df$ROC.LabelsFlip)
sift_perf <- performance(sift_pred,"tpr","fpr")
plot(sift_perf,main="SIFT ROC Plot Type A Enriched Variants")

#AUC
sift_auc <- performance(sift_pred, measure = "auc")
sift_auc@y.values[[1]]

#MCC
sift_mcc <- performance(sift_pred,"mat")
sift_mcc@y.values[[1]][which.max(sift_mcc@y.values[[1]])]

sift_df <- data.frame(ROC34$SIFT.Score,ROC34$LabelsFlip)
sift_df <- sift_df[!is.na(sift_df$ROC34.SIFT.Score),]
sift_pred <- prediction(as.numeric(as.character(sift_df$ROC34.SIFT.Score)),sift_df$ROC34.LabelsFlip)
sift_perf <- performance(sift_pred,"tpr","fpr")
plot(sift_perf,main="SIFT ROC Plot Type B Enriched Variants")

#AUC
sift_auc <- performance(sift_pred, measure = "auc")
sift_auc@y.values[[1]]

#MCC
sift_mcc <- performance(sift_pred,"mat")
sift_mcc@y.values[[1]][which.max(sift_mcc@y.values[[1]])]

#CADD (Need to normalise to range 0-1 for scores)
cadd <- normalized(as.integer(ROC$CADD.Phred.Score))
cadd_df <- data.frame(cadd,ROC$Labels)
cadd_df <- cadd_df[!is.na(cadd_df$cadd),]
cadd_pred <- prediction(as.numeric(as.character(cadd_df$cadd)),cadd_df$ROC.Labels)
cadd_perf <- performance(cadd_pred,"tpr","fpr")
plot(cadd_perf,main="CADD ROC Plot Type A Enriched Variants")
cadd_auc <- performance(cadd_pred, measure = "auc")
cadd_auc@y.values[[1]]
cadd_mcc <- performance(cadd_pred,"mat")
cadd_mcc@y.values[[1]][which.max(cadd_mcc@y.values[[1]])]

cadd <- normalized(as.integer(ROC34$CADD.Phred.Score))
cadd_df <- data.frame(cadd,ROC34$Labels)
cadd_df <- cadd_df[!is.na(cadd_df$cadd),]
cadd_pred <- prediction(as.numeric(as.character(cadd_df$cadd)),cadd_df$ROC34.Labels)
cadd_perf <- performance(cadd_pred,"tpr","fpr")
plot(cadd_perf,main="CADD ROC Plot Type B Enriched Variants")
cadd_auc <- performance(cadd_pred, measure = "auc")
cadd_auc@y.values[[1]]
cadd_mcc <- performance(cadd_pred,"mat")
cadd_mcc@y.values[[1]][which.max(cadd_mcc@y.values[[1]])]

#Mutation Assessor (also normalise)
mutass <- normalized(as.integer(ROC34$MutationAssessor.Score))
mutass_df <- data.frame(mutass,ROC34$Labels)
mutass_df <- mutass_df[!is.na(mutass_df$mutass),]
mutass_pred <- prediction(as.numeric(as.character(mutass_df$mutass)),mutass_df$ROC34.Labels)
mutass_perf <- performance(mutass_pred,"tpr","fpr")
plot(mutass_perf,main="MutAss ROC Plot Type B Enriched Variants")
mutass_auc <- performance(mutass_pred, measure = "auc")
mutass_auc@y.values[[1]]
mutass_mcc <- performance(mutass_pred,"mat")
mutass_mcc@y.values[[1]][which.max(mutass_mcc@y.values[[1]])]

mutass <- normalized(as.integer(ROC34$MutationAssessor.Score))
mutass_df <- data.frame(mutass,ROC34$Labels)
mutass_df <- mutass_df[!is.na(mutass_df$mutass),]
mutass_pred <- prediction(as.numeric(as.character(mutass_df$mutass)),mutass_df$ROC34.Labels)
mutass_perf <- performance(mutass_pred,"tpr","fpr")
plot(mutass_perf,main="MutAss ROC Plot Type B Enriched Variants")
mutass_auc <- performance(mutass_pred, measure = "auc")
mutass_auc@y.values[[1]]
mutass_mcc <- performance(mutass_pred,"mat")
mutass_mcc@y.values[[1]][which.max(mutass_mcc@y.values[[1]])]

#Box plots
boxdata <- read.csv2("boxplot_input.tsv",header=TRUE,sep="\t")

#CADD boxplot
cadd_box <- data.frame(as.numeric(as.character(boxdata$CADD.Phred.Score)),boxdata$PharmGKB.evidence.level)
#Filter out NAs
cadd_box <- cadd_box[!is.na(cadd_box$as.numeric.as.character.boxdata.CADD.Phred.Score..),]
#box-whiskers plot
ggplot(cadd_box, aes(x=cadd_box$boxdata.PharmGKB.evidence.level,y=cadd_box$as.numeric.as.character.boxdata.CADD.Phred.Score..)) + geom_boxplot(width=0.4) + ggtitle("a) CADD") + xlab("Category") + ylab("Score") + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()

#Same for polyphen2 (added ylim)
poly_box <- data.frame(as.numeric(as.character(boxdata$PolyPhen.Score)),boxdata$PharmGKB.evidence.level)
poly_box <- poly_box[!is.na(poly_box$as.numeric.as.character.boxdata.PolyPhen.Score..),]
ggplot(poly_box, aes(x=poly_box$boxdata.PharmGKB.evidence.level,y=poly_box$as.numeric.as.character.boxdata.PolyPhen.Score..)) + geom_boxplot(width=0.4) + ggtitle("b) Polyphen2") + xlab("Category") + ylab("Score") + ylim(0,1) + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()


sift_box <- data.frame(as.numeric(as.character(boxdata$SIFT.Score)),boxdata$PharmGKB.evidence.level)
sift_box <- sift_box[!is.na(sift_box$as.numeric.as.character.boxdata.SIFT.Score..),]
ggplot(sift_box, aes(x=sift_box$boxdata.PharmGKB.evidence.level,y=sift_box$as.numeric.as.character.boxdata.SIFT.Score..)) + geom_boxplot(width=0.4) + ggtitle("c) SIFT") + xlab("Category") + ylab("Score") + ylim(0,1) + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()

ma_box <- data.frame(as.numeric(as.character(boxdata$MutationAssessor.Score)),boxdata$PharmGKB.evidence.level)
ma_box <- ma_box[!is.na(ma_box$as.numeric.as.character.boxdata.MutationAssessor.Score..),]
ggplot(ma_box, aes(x=ma_box$boxdata.PharmGKB.evidence.level,y=ma_box$as.numeric.as.character.boxdata.MutationAssessor.Score..)) + geom_boxplot(width=0.4) + ggtitle("d) MutationAssessor") + xlab("Category") + ylab("Score") + ylim(-3,5) + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()

mb_box <- data.frame(as.numeric(as.character(boxdata$MutPred.Score)),boxdata$PharmGKB.evidence.level)
mb_box <- mb_box[!is.na(mb_box$as.numeric.as.character.boxdata.MutPred.Score..),]
ggplot(mb_box, aes(x=mb_box$boxdata.PharmGKB.evidence.level,y=mb_box$as.numeric.as.character.boxdata.MutPred.Score..)) + geom_boxplot(width=0.4) + ggtitle("e) MutationPredictor") + xlab("Category") + ylab("Score") + ylim(0,1) + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()

revel_box <- data.frame(as.numeric(as.character(boxdata$REVEL.Score)),boxdata$PharmGKB.evidence.level)
revel_box <- revel_box[!is.na(revel_box$as.numeric.as.character.boxdata.REVEL.Score..),]
ggplot(revel_box, aes(x=revel_box$boxdata.PharmGKB.evidence.level,y=revel_box$as.numeric.as.character.boxdata.REVEL.Score..)) + geom_boxplot(width=0.4) + ggtitle("f) REVEL") + xlab("Category") + ylab("Score") + ylim(0,1) + stat_boxplot(geom ='errorbar', width = 0.3) + theme_bw()

#Adjusted p-value calculations
# t-test to evaluate significance of differences between
# area under the curve and MCC for the 6 predictors
# for type A and type B variants

# data taken from table 2

perf <- data.frame(row.names=c("auc_typeA","mcc_typeA","auc_typeB","mcc_typeB"))
perf["PolyPhen2"] = c(0.852, 0.489, 0.397, 0.00338)
perf["MutPred"] = c( 0.975, 0.788, 0.682, 0.300)
perf["REVEL"] = c( 0.942, 0.794, 0.498 ,0.106)
perf["SIFT"] = c( 0.774, 0.358, 0.410, -0.0400)
perf["CADD"] = c( 0.728, 0.321, 0.461, 0.118)
perf["MutationAssessor"] = c( 0.763, 0.396, 0.427, 0.0764)

# perform paired t-test (1-sided)
auc_tt<-t.test(x=t(perf[1,]),y=t(perf[3,]),pair=TRUE, alternative="greater")
mcc_tt<-t.test(x=t(perf[2,]),y=t(perf[4,]),pair=TRUE, alternative="greater")

# bonferroni adjusted difference between type a and type b over all 6 methods
print("AUC ttest (Bonf. Corrected)")
print(p.adjust(auc_tt$p.value, "bonferroni",6))
# MCC corrected p-value comes in at 0.004 - robustly confirming your observation
print("MCC ttest (Bonf. Corrected)")
print(p.adjust(mcc_tt$p.value, "bonferroni",6))

