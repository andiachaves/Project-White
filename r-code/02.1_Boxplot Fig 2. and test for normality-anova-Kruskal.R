attach(sites)
print(sites)

#normality test
#Shapiro-Wilkinson normality test
shapiro.test(WP)
#for producing a normal quantile-quantile plot
qqnorm(WP)
qqline(WP)

#If p> than 0.05. Data is normal, if not, not normal

#Bartlest test of homegenty of variances
bartlett.test(WP~Site)
#If p> than 0.05. variances are homogeneous, if not, not homogeneous


#If normal One way ANOVA comparing between the habitats
summary(aov(area ~ time, data=grow))
#To do the Tukey test a posteriori
TukeyHSD(results, conf.level=0.95)
#This is to make the regular boxplot with the two variables the $ connects the file name with the column within the file
boxplot(grow$area ~ grow$time)
qqnorm(Reefvariables$Coral.Density)


#if no-normal kruskal wallis
#To do kruskal Wallis test, pecify each column withn the file so it can be read as a vector and not as data frame
kruskal.test(sites$WP,sites$Site)

library(dunn.test)
#To do Pots-hoc dunn and Kruskal Wallis in one
dunn.testresults<-dunn.test (sites$WP, g=sites$Site, method="bonferroni", kw=TRUE, label=TRUE, wrap=TRUE, rmc=FALSE, alpha=0.05)
boxplot(sites$WP ~ sites$Site, las=2)
text(3,1,"*",cex=0.5)

--------------------------
  #CORAL DENSITY TEST FOR NORMALITY VARIANCE, AND THEN DO A KRUSKAL WALLIS AND BOX PLOT
#Shapiro-Wilkinson normality test
shapiro.test(Coralfacts$Coral.Density)
#for producing a normal quantile-quantile plot
qqnorm(Coralfacts$Coral.Density)
qqline(Coralfacts$Coral.Density)
#Interpret data. If p<0.05 means data is not normal as the null hypothesis is that sampels come from a normal distribution

#Bartlest test of homegenty of variances
bartlett.test(Coralfacts$Coral.Density,Coralfacts$Habitat)

#if no-normal kruskal wallis
#To do kruskal Wallis test, pecify each column withn the file so it can be read as a vector and not as data frame
kruskal.test(Coralfacts$Coral.Density,Coralfacts$Habitat)
dunn.test:

  #To do Pots-hoc dunn and Kruskal Wallis in one
  dunn.testresults<-dunn.test (Coralfacts$Coral.Density, g=Coralfacts$Habitat, method=p.adjustment.methods, kw=TRUE, label=TRUE, wrap=TRUE, rmc=FALSE, alpha=0.05)
#Boxplot
boxplot(Coralfacts$Coral.Density ~ factor(Coralfacts$Habitat, levels=c("NR", "MR", "OR")), xlab="Habitat", ylab= "Coral density")
text(3,1,"*",cex=2)

#SAME FOR DEAD AREA ON CORAL COLONIES

#Shapiro-Wilkinson normality test
shapiro.test(Coralfacts$Dead)
#for producing a normal quantile-quantile plot
qqnorm(Coralfacts$Dead)
qqline(Coralfacts$Dead)
#Interpret data. If p<0.05 means data is not normal as the null hypothesis is that sampels come from a normal distribution

#Bartlest test of homegenty of variances
bartlett.test(Coralfacts$Dead,Coralfacts$Habitat)

#if no-normal kruskal wallis
#To do kruskal Wallis test, pecify each column withn the file so it can be read as a vector and not as data frame
kruskal.test(Coralfacts$Dead,Coralfacts$Habitat)
dunn.test:

  #To do Pots-hoc dunn and Kruskal Wallis in one
  dunn.testresults<-dunn.test (Coralfacts$Dead, g=Coralfacts$Habitat, method=p.adjustment.methods, kw=TRUE, label=TRUE, wrap=TRUE, rmc=FALSE, alpha=0.05)
#Boxplot
boxplot(Coralfacts$Dead ~ factor(Coralfacts$Habitat, levels=c("NR", "MR", "OR")), xlab="Habitat", ylab= "Coral density")
text(3,1,"*",cex=2)

