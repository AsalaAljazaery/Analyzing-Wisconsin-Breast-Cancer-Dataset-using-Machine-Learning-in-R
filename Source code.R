# load required packages
library(dplyr)
library(stringr) 
library(rpart)
library(rpart.plot)
library(RColorBrewer)
library(rattle)
library(ggplot2)
library(PerformanceAnalytics) # for chart corrlation 
library(psych)# for pairs panels 
library(ROCR) # for ROC 
library(tree) # for tree 
library(corrplot) # for corrplot
library(factoextra)#for PCA
library(gridExtra)# for PCA 


wisconsin <- read.csv('wisconsin.csv',header=TRUE,sep=',',stringsAsFactors = TRUE) #To read the data 
str(wisconsin)
dim(wisconsin) #show data dimension
head(wisconsin)

# Turn the data class into Integer 
levels(wisconsin$Class) <- c(0, 1) # assign the number I want to each level 
wisconsin$Class <- as.integer(as.character(wisconsin$Class))
head(wisconsin)

# Check datatypes 
glimpse(wisconsin)


# Frequency Identification
xtabs(~Class+ Cell.size, data = wisconsin)
xtabs(~Class+ Cl.thickness, data = wisconsin)
xtabs(~Class+  Cell.shape, data = wisconsin)
xtabs(~Class+  Marg.adhesion, data = wisconsin)
xtabs(~Class+  Epith.c.size, data = wisconsin)
xtabs(~Class+  Bare.nuclei, data = wisconsin)
xtabs(~Class+  Bl.cromatin, data = wisconsin)
xtabs(~Class+  Normal.nucleoli, data = wisconsin)
xtabs(~Class+  Mitoses, data = wisconsin)


# Correlation Chart for Means
chart.Correlation(wisconsin[, c(1:10)], histogram=TRUE, col="grey10", pch=1, main="Class Means")

#Create Correlation Chart for SE
pairs.panels(wisconsin[,c(1:10)], method="pearson",
             hist.col = "#1fbbfa", density=TRUE, ellipses=TRUE, show.points = TRUE,
             pch=1, lm=TRUE, cex.cor=1, smoother=F, stars = T, main="SE")

# Check for missing values, only Bare nuclei attribute has  16 missing values 
sapply(wisconsin, function(x) sum(is.na(x)))

# Remove any missing value 
wisconsin <- na.omit(wisconsin)
wisconsin

# The class distribution 
ggplot(wisconsin, aes(x = Class)) +
  geom_bar(aes(fill = "blue")) +
  ggtitle("Distribution of the  datset ") +
  theme(legend.position="none")

# Create Boxplot depicting the class distribution for Cell.shape Attribute  
wisconsin$Class <- as.factor(wisconsin$Class)
levels(wisconsin$Class) <- c("benign","malignant")
boxplot(Cell.shape ~Class,data = wisconsin)

# Apply T-test 
t.test(Cell.shape~Class, data=wisconsin)

# Correlation
cor(wisconsin$Cell.size  , wisconsin$Marg.adhesion  ,  method = "pearson")
cor(wisconsin$Cell.size  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Mitoses  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Cl.thickness  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Cell.shape  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Marg.adhesion  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Epith.c.size  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Bare.nuclei  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Bl.cromatin  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Normal.nucleoli  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Mitoses  , wisconsin$Class  ,  method = "pearson")
cor(wisconsin$Cell.size, wisconsin$Cell.shape  ,  method = "spearman") 

# Data Corrgram  
corrgram(wisconsin, order=NULL, lower.panel=panel.shade, upper.panel=NULL, text.panel=panel.txt,
         main="Corrgram of the data")

# Split our data into training and test sets 
set.seed(5)
index <- sample(nrow(wisconsin), 0.7 * nrow(wisconsin)) 
wisconsintrain <- wisconsin[index,] 
wisconsintest <- wisconsin[-index,]

head(wisconsintrain)
head(wisconsintest)

# The training dataset distribution  
ggplot(wisconsintrain, aes(x = Class)) +
  geom_bar(aes(fill = "blue")) +
  ggtitle("Distribution of the training datset ") +
  theme(legend.position="none")

# The testing dataset distribution  
ggplot(wisconsintest, aes(x = Class)) +
  geom_bar(aes(fill = "blue")) +
  ggtitle("Distribution of the training datset ") +
  theme(legend.position="none")


# Calculating the Correlation Coeficients and p-values
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.2, 0.1, txt2)
}

# Graph the linear correlation coef
pairs(wisconsintrain, upper.panel = panel.cor)


# plot the diagnostics diagrams 
fit_all <- lm(Class ~ ., data = wisconsin)
par(mfrow=c(2,2))
plot(fit_all)
# Step 3 - multicollinearity diagnostics
imcdiag(fit_all)


#Generalize Linear Model
model_algorithm_1 = model = glm(Class ~ Cl.thickness + 
                                  Cell.size +
                                  Cell.shape +
                                  Marg.adhesion +
                                  Epith.c.size + 
                                  Bare.nuclei  +
                                  Bl.cromatin  +
                                  Normal.nucleoli +
                                  Mitoses,
                                family=binomial(link='logit'), control = list(maxit = 50),data=wisconsintrain)

print(summary(model_algorithm_1))

# Apply the algorithm to the training sample

prediction_training = predict(model_algorithm_1,wisconsintrain, type = "response")
prediction_training = ifelse(prediction_training > 0.5, 1, 0)
error = mean(prediction_training != wisconsintrain$Class)
print(paste('Model Accuracy',1-error))

# the ROC curve 
p = predict(model_algorithm_1, wisconsintrain, type="response")
pr = prediction(p, wisconsintrain$Class)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)










#Generalize Linear Model2
model_algorithm_2 = model = glm(Class ~ Cl.thickness + 
                                  Cell.size +
                                  Cell.shape +
                                  Marg.adhesion +
                                  Epith.c.size + 
                                  Bare.nuclei  +
                                  Bl.cromatin  +
                                  Normal.nucleoli +
                                  Mitoses,
                                family=binomial(link='logit'), control = list(maxit = 50),data=wisconsintrain)

print(summary(model_algorithm_2))


# Apply the algorithm to the training sample

prediction_training = predict(model_algorithm_2,wisconsintrain, type = "response")
prediction_training = ifelse(prediction_training > 0.5, 1, 0)
error = mean(prediction_training != wisconsintrain$Class)
print(paste('Model Accuracy',1-error))

# the ROC curve 
p = predict(model_algorithm_2, wisconsintrain, type="response")
pr = prediction(p, wisconsintrain$Class)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)





# Generalize Linear Model 3
model_algorithm_3 = model = glm(Class ~ Cl.thickness + 
                                  Marg.adhesion +
                                  Bare.nuclei  +
                                  Bl.cromatin,
                                family=binomial(link='logit'), control = list(maxit = 50),data=wisconsintrain)

print(summary(model_algorithm_3))

# Apply the algorithm to the training sample
prediction_training = predict(model_algorithm_3,wisconsintrain, type = "response")
prediction_training = ifelse(prediction_training > 0.5, 1, 0)
error = mean(prediction_training != wisconsintrain$Class)
print(paste('Model Accuracy',1-error))

# the ROC curve 
p = predict(model_algorithm_3, wisconsintrain, type="response")
pr = prediction(p, wisconsintrain$Class)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)




# Accuracy 
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Model Accuracy", auc))

# Apply the algorithm to the testing sample
prediction_testing = predict(model_algorithm_3,wisconsintest, type = "response")
prediction_testing = ifelse(prediction_testing > 0.5, 1, 0)
error = mean(prediction_testing != wisconsintest$Class)
print(paste('Model Accuracy',1-error)) 

# Get the ROC curve for testing sample 
p = predict(model_algorithm_3, wisconsintest, type="response")
pr = prediction(p, wisconsintest$Class)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)


# using confusion matrix
table(wisconsintest$Class, prediction_testing)

# Accuracy 
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Model Accuracy", auc))


#  Decision Tree 1
treemodel <- rpart(Class~., data=wisconsintrain)
plot(treemodel, margin=0.25)
text(treemodel, use.n=T)

# Tree model Prediction 
fancyRpartPlot(treemodel)
prediction <- predict(treemodel, newdata=wisconsintest, type='vector')
table(prediction, wisconsintest$Class)


# Decision Tree 2
model_tree = tree(Class ~  
                    Cell.size +
                    Cell.shape +
                    Bare.nuclei  +
                    Bl.cromatin,
                  data = wisconsintrain)

summary(model_tree)

plot(model_tree, type = "uniform")
text(model_tree, pretty = 0, cex=0.8)

model_tree_pred_train = predict(model_tree, wisconsintrain) # gives the probability for each class
model_tree_pred_test = predict(model_tree, wisconsintest) # gives the probability for each class


# # Try to prune the tree to avoid over fitting
cv.tree(model_tree)

plot(cv.tree(model_tree)) # Seems like a tree of size 5 might be best
model_tree_prune = prune.tree(model_tree, best = 5)
summary(model_tree_prune)

plot(model_tree_prune, type = "uniform")
text(model_tree, pretty = 0, cex=0.8)





# PCA
all_pca <- prcomp(wisconsin[,1:9], cor=TRUE, scale = TRUE)
summary(all_pca)

# PCA 
fviz_eig(all_pca, addlabels=TRUE, ylim=c(0,60), geom = c("bar", "line"), barfill = "pink",  
         barcolor="blue",linecolor = "red", ncp=10)+
  labs(title = " All Variances - PCA",
       x = "Principal Components", y = "% of variances")

all_var <- get_pca_var(all_pca)
all_var
corrplot(all_var$cos2, is.corr=FALSE)

# Variables contributions
p1 <- fviz_contrib(all_pca, choice="var", axes=1, fill="pink", color="grey", top=10)
p2 <- fviz_contrib(all_pca, choice="var", axes=2, fill="skyblue", color="grey", top=10)
grid.arrange(p1,p2,ncol=2)

# cluster
library("devtools")
set.seed(218)
res.all <- kmeans(all_var$coord, centers = 2, nstart = 25)
grp <- as.factor(res.all$cluster)

fviz_pca_var(all_pca, col.var = grp, 
             palette = "jco",
             legend.title = "Cluster")

# PCA biplot
fviz_pca_biplot(all_pca, col.ind = as.factor(wisconsin$Class), col="black",
                palette = "jco", geom = "point", repel=TRUE,
                legend.title="Diagnosis", addEllipses = TRUE)


