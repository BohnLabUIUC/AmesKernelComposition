#This is the final workflow for the phenotypic analysis of kernel compositional traits
# I conducted an analysis of variance between heterotic groups since this wasnt possible between inbreds due to lack of reps
#I used a simple linear model which assumes fixed factors, since we dont have random factors. Heterotic groups are  fixed groups thus considered fixed

Compiled <- read.csv("GrainQualityData_AmesPanel.csv", header = TRUE)

Compiled$Accession <- as.factor(Compiled$Accession)
Compiled$Origin <- as.factor(Compiled$Origin)
Compiled$Group <- as.factor(Compiled$Group)
Compiled$Starch <- as.numeric(Compiled$Starch)
Compiled$Oil <- as.numeric(Compiled$Oil)
Compiled$Protein <- as.numeric(Compiled$Protein)
Compiled$Fib <- as.numeric(Compiled$Fib)
Compiled$Density <- as.numeric(Compiled$Density)
Compiled$Ash <- as.numeric(Compiled$Ash)



## Testing for outliers
boxplot(Compiled[4:9], plot=TRUE)$out


# create detect outlier function

outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25) # calculate first quantile
  Q3 <- quantile(x, probs=.75) # calculate third quantile
  iqr = Q3-Q1           # calculate inter quartile range
  
  upper_limit = Q3 + (iqr*1.5)  #Calculate upper limit
  lower_limit = Q1 - (iqr*1.5)  #calculate lower limit
  
  x > upper_limit | x < lower_limit   # return true or false
}

remove_outliers <- function(Compiled, cols = names()) {  # for loop to traverse in columns vector
  for (col in cols) {
    Compiled <- Compiled[!outliers(Compiled[[col]]),] # remove observation if it satisfies outlier function
  }
  Compiled
}

#Remove outliers 
Compiled2 <- remove_outliers(Compiled, c(4:9))

outliers <- boxplot(Compiled[4:9], plot=TRUE)$out  # save the outliers in a vector
x<-Compiled
Compiled <- x[-which(x$Starch %in% outliers),] #Removing outlier for just one variable
Compiled <- Compiled[-which(Compiled$Protein %in% outliers),]
Compiled <- Compiled[-which(Compiled$Oil %in% outliers),]
Compiled <- Compiled[-which(Compiled$Fib %in% outliers),]
Compiled <- Compiled[-which(Compiled$Ash %in% outliers),]
Compiled <- Compiled[-which(Compiled$Density %in% outliers),]

#Boxplot
boxplot(Compiled2[4:8], 
        col=c("red", "yellow", "green","#999999", "#E69F00", "#56B4E9"), plot=TRUE)$out 
boxplot(Compiled2[9], col= "blue", plot=TRUE)$out 

## Trait Correlation
#Rename Fib to Fiber 
Corr_data <- Compiled
library(dplyr)
Corr_data <- Corr_data %>%
  rename(Fiber = Fib)

##plot correlations #Figure 1B
library("PerformanceAnalytics")
chart.Correlation(Corr_data[4:9], histogram=TRUE, pch=19)  


##Alternative way
library(car)
library(corrplot)
Correlation <- cor(Corr_data[4:9], use = "pairwise.complete.obs")  #for both IL, IN and WS
#png("Trait_Correlations.png", width = 20, height = 20, units = "cm", res = 300)
corplot <- corrplot(Correlation,method = "color",type="upper", order="hclust", #Type = upper,lower, #method=circle,pie,color,number
                    addCoef.col="black", # Add coefficient of correlation
                    diag=FALSE, # hide correlation coefficient on the principal diagonal
                    tl.col="black", tl.srt=45, #Text label color and rotation
                    p.mat=NULL, sig.level = 0.01, insig = "blank")  # Add coefficient of correlation
print(corplot)

### spearman Rank Correlations between our values and grain quality values obtained by Hirsch et al, 2021. 
## Adjusted Blup values were extracted from Table 2 of the paper, and correlated to our estimated blups for each group 

#Estimate trait Blups for each Group using a random effect model
### Mixed Models 

library(lme4)
library(jtools)
library(lmerTest)
library(car)
install.packages("xlsx")
library(xlsx)
library(readr)

mod1 <- lmer(Protein ~ (1|Group), REML = TRUE, data= Compiled)
mod2 <- lmer(Starch ~ (1|Group), REML = TRUE, data= Compiled)
mod3 <- lmer(Oil ~ (1|Group), REML = TRUE, data= Compiled)
mod4 <- lmer(Ash ~ (1|Group), REML = TRUE, data= Compiled)
mod5 <- lmer(Fib ~ (1|Group), REML = TRUE, data= Compiled)
mod6 <- lmer(Density ~ (1|Group), REML = TRUE, data= Compiled)


## extract blups from each model
varComp<-as.data.frame(VarCorr(mod1,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
blupProt= coef(mod1)$Group

varComp<-as.data.frame(VarCorr(mod3,comp="vcov"))
blupOil = coef(mod3)$Group

varComp<-as.data.frame(VarCorr(mod2,comp="vcov"))
blupStr = coef(mod2)$Group

varComp<-as.data.frame(VarCorr(mod5,comp="vcov"))
blupfiber = coef(mod5)$Group

varComp<-as.data.frame(VarCorr(mod4,comp="vcov"))
blupash = coef(mod4)$Group

varComp<-as.data.frame(VarCorr(mod6,comp="vcov"))
blupdensity = coef(mod6)$Group

Ames_Blups <- cbind(blupProt,blupStr,blupOil,blupfiber,blupdensity,blupash)

write.csv(Ames_Blups, "Ames_Blups.csv")

# Test for normality Assumptioms
# qq plot
par(mfrow=c(3,2))
qqPlot(residuals(mod1), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Starch Content") 
qqPlot(residuals(mod2), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Protein Content", id = FALSE) 
qqPlot(residuals(mod3), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Oil Content", id = FALSE) 
qqPlot(residuals(mod4), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Ash Content", id = FALSE) 
qqPlot(residuals(mod5), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Fiber Content", id = FALSE) 
qqPlot(residuals(mod6), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Kernel Density", id = FALSE) 

## Fit a simple linear model to determine subpopulation group effect
### Only model with Group can be used. We have reps within groups. 
#We cant fit a model with genotype since each genotype appears once (one observation)

modelprot <- aov(Protein ~ Group, Compiled)
modelstr <- aov(Starch ~ Group, Compiled)
modeloil <- aov(Oil ~ Group, Compiled)
modelash <- aov(Ash ~ Group, Compiled)
modeldensity <- aov(Density ~ Group, Compiled)
modelfiber <- aov(Fib ~ Group, Compiled)

anova(modelprot)
anova(modelstr)
anova(modeloil)
anova(modelfiber)
anova(modelash)
anova(modeldensity)

## Pairwise Mean Comparisons for each trait. 
library(agricolae)
out_prot <-HSD.test(modelprot,"Group", group=TRUE)
prot_mean <- as.data.frame(out_prot$groups)
prot_mean$Group <- rownames(prot_mean)

out_str  <-HSD.test(modelstr,"Group", group=TRUE)
str_mean <- as.data.frame(out_str$groups)
str_mean$Group <- rownames(str_mean)

out_oil  <-HSD.test(modeloil,"Group", group=TRUE)
oil_mean <- as.data.frame(out_oil$groups)
oil_mean$Group <- rownames(oil_mean)

out_fib  <-HSD.test(modelfiber,"Group", group=TRUE)
fib_mean <- as.data.frame(out_fib$groups)
fib_mean$Group <- rownames(fib_mean)

out_ash  <-HSD.test(modelash,"Group", group=TRUE)
ash_mean <- as.data.frame(out_ash$groups)
ash_mean$Group <- rownames(ash_mean)

out_den  <-HSD.test(modeldensity,"Group", group=TRUE)
den_mean <- as.data.frame(out_den$groups)
den_mean$Group <- rownames(den_mean)

library(dplyr)
Ames_means2 <- list(prot_mean, str_mean, oil_mean, fib_mean, den_mean, ash_mean) %>%
  reduce(full_join, by = "Group")

Ames_means <- Ames_means2[c(3,1,4,6,8,10,12)]
write.csv(Ames_means, "Ames_means")

#The trait blups and means were combined with the Hirsch et al, 2021 Blups from Table 2 to for "TableValues" dataset used in the next step
#The dataset has the Group, Trait, Ames_Blups for each trait and Wisc_Blups for each trait. 
#Group performance/rankig for a particular trait across the studies was evaluated using a spearman ranking test below. 

Data1  <- read.csv("TableValue.csv", header = TRUE)
library(dplyr) 
##correlation with blups
Corr1 <- Data1%>% 
  group_by (Trait) %>% 
  summarise(cor=cor(Ames_Blups,Wisc_Blups, method = "spearman"))
Corr1

##correlation with means
Corr2 <- Data1%>% 
  group_by (Trait) %>% 
  summarise(cor=cor(Ames_Means,Wisc_Blups, method = "spearman"))
Corr2


### Pearson correlation using mean values of the 275 genotypes that where phenotyped in both studies
Wis <- read.csv("Wis-Rawdata.csv", header = TRUE)  #Subset of Hirsch et al, 2021 raw-dataset with only 275 common, from all tested environment 

Wis$Accession <- as.factor(Wis$Accession)
Wis$Genotype <- as.factor(Wis$Genotype)
Wis$Group <- as.factor(Wis$Group)
Wis$Block <- as.factor(Wis$Block)
Wis$Rep <- as.factor(Wis$Rep)
Wis$Env <- as.factor(Wis$Env)
Wis$Star_W <- as.numeric(Wis$Star_W)
Wis$Oil_W <- as.numeric(Wis$Oil_W)
Wis$Prot_W <- as.numeric(Wis$Prot_W)
Wis$Fib_W <- as.numeric(Wis$Fib_W)
Wis$Ash_W <- as.numeric(Wis$Ash_W)

traits <- c("Star_W", "Oil_W", "Prot_W", "Fib_W", "Ash_W")

## Estimate BLUPs for common genotypes in Renk dataset using model used in Renk et al. 2021
# Create an empty list to store models
models <- list()

# Loop through each trait and fit the model
for (trait in traits) {
  formula <- as.formula(paste(trait, "~ (1|Accession) + (1|Env) + (1|Env/Rep) + (1|Env/Rep/Block) + (1|Accession:Env)"))
  models[[trait]] <- lmer(formula, data = Wis, REML = TRUE)
}

#estimate blups
blup_list <- list()

# Loop through traits and extract Genotype BLUPs
for (trait in names(models)) {
  blup <- coef(models[[trait]])$Accession
  colnames(blup) <- trait  # Rename column to trait name
  blup$Accession <- rownames(blup)  # Save Accession as a column
  blup_list[[trait]] <- blup
}

# Merge all BLUPs by Genotype
library(dplyr)
library(purrr)

# Reduce the list of data.frames by joining them on "Accession"
combined_blups <- reduce(blup_list, full_join, by = "Accession")

# Reorder columns to have Genotype first
combined_blups <- combined_blups %>%
  select(Accession, everything())


library(dplyr)
WisMeans <- Wis %>%
  group_by(Accession) %>%
  summarise(across(Prot_W:Star_W, mean)) #Means of the 275 genotypes 
write.csv (WisMeans, "WisMeans.csv") 

##Combining the 2 datsets from the 2 studies 
#Combined Dataset with grain quality trait values from both studies

WisAmesCombined <- merge(Compiled, combined_blups, by = 'Accession') #merge function joins the 2 datasets by the common accession 
write.csv (WisAmesCombined, "Combined2.csv") 

## Pearson Correlation matrix of the observed and literature values for the 275 common genotypes. 

library(Hmisc) #Create the function to generate the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )}
res3<-rcorr(as.matrix(WisAmesCombined[,4:14]),type = "pearson")
round(res3$r, 2)
round(res3$P, 2)
Corr3 <- flattenCorrMatrix(res3$r, res3$P)

write.csv (Corr2, "PearsonCorr.csv")


#Plotting Correlations
library("ggpubr")
Starch <- ggscatter(WisAmesCombined, x = "Starch", y = "Star_W", 
          add = "reg.line", conf.int = TRUE, color = "Group", 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Starch",
          xlab = "Observed [%]", ylab = "Literature [%]")


Protein <- ggscatter(WisAmesCombined, x = "Protein", y = "Prot_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Protein",
                    xlab = "Observed [%]", ylab = "Literature  [%]")

Oil <- ggscatter(WisAmesCombined, x = "Oil", y = "Oil_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Oil",
                 xlab = "Observed [%]", ylab = "Literature  [%]")

Fiber <- ggscatter(WisAmesCombined, x = "Fib", y = "Fib_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Fiber",
                    xlab = "Observed [%]", ylab = "Literature  [%]")
                  

Ash <- ggscatter(WisAmesCombined, x = "Ash", y = "Ash_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Ash",
                    xlab = "Observed [%]", ylab = "Literature  [%]")

library(gridExtra)
grid.arrange(Starch, Protein, Oil, Fiber,Ash, nrow = 3)


table(WisAmesCombined$Group)

#Colored scatters

Starch <- ggscatter(
  WisAmesCombined, 
  x = "Starch", 
  y = "Star_W", 
  color = "Group",              # color dots by group
  add = "reg.line",              # one regression line
  add.params = list(color = "black"), # set line color manually
  conf.int = TRUE, 
  cor.coef = TRUE, 
  cor.method = "pearson",
  title = "Starch",
  xlab = "Observed [%]", 
  ylab = "Literature [%]")+
  theme(legend.position = "none") 


Protein <- ggscatter(
                WisAmesCombined, 
                x = "Protein", 
                y = "Prot_W", 
                color = "Group",              # color dots by group
                add = "reg.line",              # one regression line
                add.params = list(color = "black"), # set line color manually
                conf.int = TRUE, 
                cor.coef = TRUE, 
                cor.method = "pearson",
                title = "Protein",
                xlab = "Observed [%]", 
                ylab = "Literature [%]")+
           theme(legend.position = "none")  


Oil        <- ggscatter(
                     WisAmesCombined, 
                     x = "Oil", 
                     y = "Oil_W", 
                     color = "Group",              # color dots by group
                     add = "reg.line",              # one regression line
                     add.params = list(color = "black"), # set line color manually
                     conf.int = TRUE, 
                     cor.coef = TRUE, 
                     cor.method = "pearson",
                     title = "Oil",
                     xlab = "Observed [%]", 
                     ylab = "Literature [%]")+
            theme(legend.position = "none")  



Fiber      <- ggscatter(
                   WisAmesCombined, 
                   x = "Fib", 
                   y = "Fib_W", 
                   color = "Group",              # color dots by group
                   add = "reg.line",              # one regression line
                   add.params = list(color = "black"), # set line color manually
                   conf.int = TRUE, 
                   cor.coef = TRUE, 
                   cor.method = "pearson",
                   title = "Fiber",
                   xlab = "Observed [%]", 
                   ylab = "Literature [%]")+
           theme(legend.position = "none")  



Ash       <- ggscatter(
                     WisAmesCombined, 
                     x = "Ash", 
                     y = "Ash_W", 
                     color = "Group",              # color dots by group
                     add = "reg.line",              # one regression line
                     add.params = list(color = "black"), # set line color manually
                     conf.int = TRUE, 
                     cor.coef = TRUE, 
                     cor.method = "pearson",
                     title = "Ash",
                     xlab = "Observed [%]", 
                     ylab = "Literature [%]")+
          theme(legend.position = "none")  


library(gridExtra)
grid.arrange(Starch, Protein, Oil, Fiber,Ash, nrow = 3)




