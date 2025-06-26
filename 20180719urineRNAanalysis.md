Urine expression data classifier
============================

Path to actual code is /mnt/DATA4/kevhu/scripts/20180715randoForest.R, /home/kevhu/scripts/20180808urineClassifierCleaner.R & /home/kevhu/scripts/20180808urineClassifierGaneshGG2.R . However this is a lot cleaner to look at and better organized

Setup steps and packages needed for the analysis. Expression data is divided by KLK3 * 100,000 for normalization. Similar normalization is done in the original Mips paper. Table a is from Andi and b is from Scott. Same information, but Scott's annotations for cancer non-cancer vs what I'm using which is division by grade groups. I also use Scott's data because he updated it with serum PSA measurements
```{r, echo=TRUE}
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(VSURF)
library(parallel)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(stringr)
library(glmnet)
library(Rcpp)

a <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/AC17.2ExtDes2_mastermatrix_allsamples_byKLK3.xlsx")
b <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/20189_07_26_AC17.2ExtDes2_mastermatrix.xlsx")
```

Steps to first remove grade group 2 samples, afterwards I will set aside a set of samples for cross-validation later. Plus some steps to setup the work the VSURF format which is a list of two things: (1) data frame with samples and variables as rows and columns respsectively (2) vector of classes as factors. 

```{r}
a <- a[-which(a$`Grade group` == 2),]
a.1 <- a[1:109,3:86]
rownames(a.1) <- as.vector(a$Contig_ID[1:109])
ids <- a$`Grade group bin`[1:109]

b <- b[match(rownames(a.1), b$Contig_ID),]
PSA <- b$LAB_RESNUM

set.seed(2222)
g1 <- which(ids == "Benign/GG1")
g2 <- which(ids == "GG2-5")

g1.subset <- sample(g1, length(g1)/3)
g2.subset <- sample(g2, length(g2)/3)
combined.subset <- c(g1.subset,g2.subset)

held.out <- a.1[combined.subset,]
a.2 <- a.1[-combined.subset,]

ids.2 <- factor(ids[-combined.subset])
ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")
```

Running the random forest model and evaluating the results. I only use variables from thresholding step, because the next two steps are recurvsive feature elimination (RFE) .. which I'm not a fan of, and fitting to a random forest. Diaz-Urute 2011(?) popularized random forests into RFE in their microarray paper. I think using random forest to calculate VI from MSE decrease followed by CART on the standard deviations is enough - CART on standard deviation of VI is not common. There aren't too many ways of selecting from the ranked list of VIs though, so this a more objective way of choosing variables. Again, the purpose of the random forest was just to reduce correlated feutures i.e genes because correlated features cause models to overfit. Note a specific type of seed for the parallel processing prior to random forest will make the steps reprodicble - need to reset seed after every run for it to work as intended. Also the VSURF function leverages the parallel package and the parallel=TRUE runs on our server for parallel processing

```{r,  fig.keep= "all"}
testRandom <- list(a.2, ids.2)
names(testRandom) <- c("x","y")
set.seed(2222, kind = "L'Ecuyer-CMRG")
vsurf.parallel <- VSURF(testRandom$x, testRandom$y, mtry = 100, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf.parallel)
colnames(a.1)[vsurf.parallel$varselect.interp]
plot(vsurf.parallel)

```

Below I load the genes found to be informative in the random forest and subset them. One can explore the nature of these genes by using caret's (Max Kuhn) featurePlot. Note, featurePlot is lattice (package) based, so graphical parameters paased need to be lattice ones. I am not too familiar with lattice graphical parameters, so I also include ggplot code for boxplots

```{r, fig.keep= "all"}

genesThres <- colnames(a.2)[vsurf.parallel$varselect.thres]
a.2.genesThres <- a.2[,genesThres]
a.2.genesThres <- a.2.genesThres[, order(colnames(a.2.genesThres))]

featurePlot(x = a.2.genesThres,
            y = ids.3,
            plot = "box",
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")),
            par.strip.text=list(cex=.6),
            pch= "|")


load("/mnt/DATA4/kevhu/urineRNA/20180807randoForestGenes.Robj")

a.2.genes30 <- a.2[,genes30]
a.2.genes30 <- a.2.genes30[, order(colnames(a.2.genes30))]

a.2.log <- log2(a.1[-combined.subset,]+1)
a.2.log.30 <- cbind(ids.3,a.2.log[,genes30])
a.2.log.30.melt <- melt(data = a.2.log.30,id.vars = "ids.3")
colnames(a.2.log.30.melt)[1] <- "Classes"

ggplot(data = a.2.log.30.melt, aes(x = variable, y = value, colour = Classes)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free") + theme_bw() + ylab(label = "log2 normalized expresssion") + 
  ggtitle(label = "Urine panel: 30 selected genes") + xlab(label = "") +
  theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5),panel.grid = element_blank(), legend.title = element_text(hjust = 0.5)) +
  theme(legend.key.size = unit(0.5, "inches"), title = element_text(face = "bold"), strip.text.x = element_text(size = 7.5))



```


Example of regularized regression which is does MLE-based penalty. Done using the glmnet package which is elastic net. Elastic net allows for penalty that is note purely either Lasso or ridge. Note x has to be just a matrix of predictors, whereas y is the response in the form of factors. I also think you can only tune it to maximize AUC with two class classification ... multi-class is not supported, becuase I guess since there is n (# of classes) choose 2 AUC curves needed to tune and it gets complex .... The estimated parameters of alpha and lambda should be almost converge

```{r}
a.3 <- cbind(ids.3,a.2.genes30)

test.cv <- cv.glmnet(x = matrix(unlist(a.3[,2:ncol(a.3)]), ncol = ncol(a.3)-1), y=a.3[,1], family = c("binomial"),nfolds = 5,
                     type.measure = "auc", maxit=1000000, alpha = 1)
```

Using caret for five-fold cross validation. Afterwards, I can use the held out set to do additional cross-validation. Note held out set was not used for random forest variable selection. Data was split 2/3 training and 1/3 testing. Caret does the auto conversion needed for one's response and predictors needed when it passes through to glmnet. Glmnet standardizes the predictor variables, not reponse variables, as a default setting so no preprocessing needs to be done. Standardizing data such that your predictors have mean=0 and var=1 is pretty standard as it rids the intercept or zeroes it for regularized regression. If not, your scales may be off and the shrinking of variables can be skewed. Reference Tibshirani 2002 for original Lasso and Ridge paper. Also I use log2 normalized data as it seems more stable. Weights are also needed since the data set i.e Cancer vs non-cancer is unbalanced. Typically 3 ways to deal with this: upsample to majority class, downsample to minority class, or provide class weights. All three methods should typically yield similar results.


```{r}
a.7 <- a.2.log.30
a.7 <- a.7[order(a.7$ids.3),]
ids.4 <- ids.3
ids.3 <- ids.3[order(ids.3)]


train_control.combined <- trainControl(method = "cv", number = 5, savePredictions = TRUE,
                                       classProbs = TRUE, summaryFunction = twoClassSummary)
model_weights <- ifelse(a.7$ids.3 == "Cancer",
                        (1/table(a.7$ids.3)[1]) * 0.5,
                        (1/table(a.7$ids.3)[2]) * 0.5)

set.seed(2222)
fiveFoldModel.log <- train(ids.3 ~., data = a.7, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.log
```

Below I am making ROC curves from pROC package. Done b/c caret utilizes the roc function from pROC when maximizing for AUC. Then I do validation on the testing cohort. Note the samples needed to be reordered for training b/c of the order of the weights - not needed for testing set. The levels need not be in the same order as the pROC package should detect this, shoot a warning and relevel for you

```{r}
####alpha = 0.55 and lambda = 0.04765286.
### makign auc for testing
testingPred.log <- fiveFoldModel.log$pred
testingPred.log <- testingPred.log[which(testingPred.log$lambda > 0.04765),]
testingPred.log <- testingPred.log[which(testingPred.log$alpha == 0.55),]
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer))
### testing set for log

ids.held.out <- ids[combined.subset]
ids.held.out <- str_replace_all(ids.held.out, "GG2-5", "Cancer")
ids.held.out <- str_replace_all(ids.held.out, "Benign/GG1", "NoCancer")
ids.held.out <- as.factor(ids.held.out)
held.out.log <- log2(held.out[,genes30] + 1)

prob.log <- predict(fiveFoldModel.log, newdata = held.out.log)
confusionMatrix(data = factor(prob.log), ids.held.out, positive = "Cancer")
probROC.log <- predict(fiveFoldModel.log, type = c("prob"),newdata = held.out.log)
rocGraph.log <- roc(ids.held.out ~ probROC.log$`Cancer`)
plot.roc(rocGraph.log, main = "PSA with model ROC")
rocGraph.log$auc
```

Next just redid the same model fitting and performance tuning for PSA, 29 gene + PSA and Mips

```{r}

PSA.training <- PSA[-combined.subset]
PSA.training <- PSA.training[order(ids.4)]
a.5 <- data.frame(ids.3,PSA.training)
a.5$ones <- rep(1,nrow(a.5))


set.seed(2222)
fiveFoldModel.psa <- train(ids.3 ~., data = a.5, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.psa

### makign auc for testing
testingPred.psa <- fiveFoldModel.psa$pred
testingPred.psa <- testingPred.psa[which(testingPred.psa$lambda > 0.0310),]
testingPred.psa <- testingPred.psa[which(testingPred.psa$alpha < 0.11),]
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer))
###

held.out.30.psa <- data.frame(PSA[combined.subset])
held.out.30.psa <- cbind(ids.held.out,held.out.30.psa)
held.out.30.psa$ones <- rep(1, nrow(held.out.30.psa))
colnames(held.out.30.psa) <- colnames(a.5)

prob.psa <- predict(fiveFoldModel.psa, newdata = held.out.30.psa)
confusionMatrix(data = factor(prob.psa), ids.held.out, positive = "Cancer")
probROC.psa <- predict(fiveFoldModel.psa, type = c("prob"),newdata = held.out.30.psa)
rocGraph.psa <- roc(ids.held.out ~ probROC.psa$`Cancer`)
plot(rocGraph.psa)
rocGraph.psa$auc

Combined

a.6 <- cbind(a.7[,1],PSA.training, a.7[,2:ncol(a.7)])
colnames(a.6)[1:2] <- c("ids.3","PSA")

set.seed(2222)
fiveFoldModel.combined <- train(ids.3 ~., data = a.6, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.combined

### making auc for testing
testingPred.combined <- fiveFoldModel.combined$pred
testingPred.combined <- testingPred.combined[which(testingPred.combined$lambda > 0.0476),]
testingPred.combined <- testingPred.combined[which(testingPred.combined$alpha == 0.55),]
plot.roc(pROC::roc(testingPred.combined$obs ~ testingPred.combined$Cancer))
###

held.out.combined <- cbind(ids.held.out,held.out.30.psa[,2],held.out.log)
colnames(held.out.combined) <- colnames(a.6)
prob.combined <- predict(fiveFoldModel.combined, newdata = held.out.combined)
confusionMatrix(data = factor(prob.combined), ids.held.out, positive = "Cancer")
probROC.combined <- predict(fiveFoldModel.combined, type = c("prob"),newdata = held.out.combined)
rocGraph.combined <- roc(ids.held.out ~ probROC.combined$`Cancer`)
plot(rocGraph.combined)
rocGraph.combined$auc
```

Mips is a bit different since it was done on a different technology ... we can try and make our own using the same number of features i.e PCA3 (lncRNA), TMRPSS2:ERG, and PSA. For TMPRSS2:ERG fusions. I just used sum of all the TMPRSS2:ERG fusions as one variable. In the future I will attempt to do weighted mean by frequency.

```{r}
tmprss2.fusion.training <- log2(rowSums(a.2[,(which(grepl("TMPRSS2-ERG", colnames(a.2))))]) + 1)
tmprss2.fusion.training <- tmprss2.fusion.training[order(ids.4)]
PCA3 <- as.vector(log2(a.2[,which(grepl("PCA3", colnames(a.2)))] + 1))
PCA3 <- PCA3$PCA3.E2E3[order(ids.4)]
mips.genes <- data.frame(cbind(PSA.training, PCA3, tmprss2.fusion.training),stringsAsFactors = FALSE)
colnames(mips.genes) <- c("PSA","PCA3","TMPRSS2:ERG")
mips.genes <- cbind(ids.3, mips.genes) 

set.seed(2222)
fiveFoldModel.mips <- train(ids.3 ~., data = mips.genes, trControl=train_control.combined, method="glmnet",
                                metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.mips

### making auc for testing
testingPred.mips <- fiveFoldModel.mips$pred
testingPred.mips <- testingPred.mips[which(testingPred.mips$lambda > 0.00348),]
testingPred.mips <- testingPred.mips[which(testingPred.mips$alpha == 0.1),]
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer))


held.out.tmrpss <- log2(rowSums(held.out[,(which(grepl("TMPRSS2-ERG", colnames(held.out))))])+1)
held.out.pca <- log2(held.out[,which(grepl("PCA3", colnames(held.out)))] + 1)
held.out.mips <- cbind(ids.held.out,
                       held.out.30.psa[,2],
                       held.out.pca,
                       held.out.tmrpss)
colnames(held.out.mips) <- colnames(mips.genes)
prob.mips <- predict(fiveFoldModel.mips, newdata = held.out.mips)
confusionMatrix(data = factor(prob.mips), ids.held.out, positive = "Cancer")
probROC.mips <- predict(fiveFoldModel.mips, type = c("prob"),newdata = held.out.mips)
rocGraph.mips <- roc(ids.held.out ~ probROC.mips$`Cancer`)
plot(rocGraph.mips)
rocGraph.mips$auc
```

AUC plots were made using base R ..... because pROC makes plots using base R, so it's easier to manipulate this way

```{r}
plot(rocGraph.combined, main = "ROCs: Training n=73; Testing n=36", col = alpha("green",alpha = 0.5), xlim = c(1,0), ylim = c(0,1), lty = 2, asp = FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.combined$obs ~ testingPred.combined$Cancer), col = alpha("green", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(rocGraph.psa, col = alpha("red", alpha = 0.5), main = NULL, ylab = "", xlab="", lty = 2, asp = FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer), col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(rocGraph.log, col = alpha("blue", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer), col = alpha("blue", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(rocGraph.mips, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer), col = alpha("orange", alpha = 0.5), asp = FALSE)
legend("right", legend = c("29 genes - AUC:0.899 ; AUC:0.812", "PSA - AUC:0.647 ; AUC:0.625",
                           "29 genes + PSA - AUC:0.899 ; AUC:0.812", "Mips: - AUC:0.749 ; AUC:0.63",
                           "Training","Testing"),
       col = c("blue", "red","green","orange","black","black"), lty=c(rep(1,5),2), cex = 1.0)

```


Final notes, Scott's annotations performed a bit worse than the model here which was GG1/benign vs GG3-5. The model vs GG2-5 is in another script I linked and not included because it is just repeptitive. 



