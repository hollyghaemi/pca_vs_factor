rm(list=ls()) #clears variables
dev.off() #clears all plots


################
## Question 1 ##
################

#packages & libraries
library(psych)

#data clean up
affect<- cbind(affect$EA1, affect$TA1, affect$PA1, 
               affect$NA1, affect$EA2, affect$TA2, 
               affect$PA2, affect$NA2)
affect<- data.frame(affect)

#outlier ID
outlier(affect) #Mahalonobis distance

#Covariance PCA with outliers
affectPC<- princomp(affect)
biplot(affectPC, col = c(2,3), cex= c(.1, .50),
       xlim = c(-.2, .2), 
       xlab = 'First Principal Component', 
       ylab = 'Second principal component', 
       main = "Biplot: Covariance with outliers")
summary(affectPC)
affectPC$loadings

#Covariance PCA without outliers
library(sos); findFn("outlier") #outlier removal
affect.outliersremoved<- affect[!apply(sapply(affect, function(x) abs(scale(x)) >= 2), 1, any), ]
affectPC.nooutliers<- princomp(affect.outliersremoved)
biplot(affectPC.nooutliers, col = c(2,3), cex= c(.1, .50),
       xlim = c(-.2, .2), 
       xlab = 'First Principal Component', 
       ylab = 'Second principal component', 
       main = "Biplot: Covariance without outliers")
summary(affectPC.nooutliers)
affectPC.nooutliers$loadings

sapply(affect, sd)
sapply(affect.outliersremoved, sd)

#Correlation PCA with outliers
affectPC.corr<- princomp(affect, cor = TRUE)
biplot(affectPC.corr, col = c(2,3), cex= c(.1, .50),
       xlim = c(-.2, .2), 
       xlab = 'First Principal Component', 
       ylab = 'Second principal component', 
       main = "Biplot: Correlation with outliers")
summary(affectPC)
affectPC.corr$loadings

#Correlation PCA without outliers
affectPC.corr.noout<- princomp(affect.outliersremoved, cor = TRUE)
biplot(affectPC.corr.noout, col = c(2,3), cex= c(.1, .50),
       xlim = c(-.2, .2), 
       xlab = 'First Principal Component', 
       ylab = 'Second principal component', 
       main = "Biplot: Correlation without outliers")
summary(affectPC.corr.noout)
affectPC.corr.noout$loadings

################
## Question 2 ##
################

age<- c(81, 95, 95, 165, 286, 299, 380, 418, 420, 547, 590, 
        635, 752, 760, 1171, 1277, 1520, 2138, 3626)
RB <- c(2.0, 6.5, 3.6, 1.9, 2.6, 2.9, 1.9, 7.1, 6.4, 6.4, 
        1.8, 6.7, 1.8, 7.3, 1.8, 1.3, 4.0, 0.5, 4.2)
p16<- c(3.07, 1.90, 3.82, 3.74, 5.17, 5.76, 2.40, 3.38, 
        3.37, 4.05, 5.15, 2.67, 3.28, 0.92, 6.56, 0.05, 2.79, 0.00, 4.24)
DLK<- c(308975, 70988, 153061, 596992, 369601, 1119258, 214071, 69511, 
        81457, 64348, 164881, 126016, 567858, 43438, 716260, 94, 31125, 2331, 560208)
Nanog<- c(94, 382, 237, 88, 282, 177, 45, 265, 659, 336, 2012, 3072, 127, 
          698, 392, 15, 454, 33, 340)
Cmyc<- c(6.49, 1.00, 0.00, 0.00, 12.23, 8.76, 5.76, 1.17, 1.88, 0.78, 35.65, 
         0.00, 4.13, 1.77, 12.92, 0.36, 0.62, 0.03, 5.43)
Ezh2<- c(2.76, 7.09, 5.57, 2.47, 1.63, 3.51, 1.41, 3.07, 3.87, 4.76, 9.45, 
         4.35, 1.00, 3.32, 2.90, 3.83, 2.33, 0.17, 1.36)
igf2<-c(11176, 5340, 6310, 7009, 7104, 9342, 3726, 8039, 12583, 6505, 32722, 11763, 
        10283, 11518, 13264, 30, 1163, 66, 21174)

hemangioma<- cbind(age, RB, p16, DLK, Nanog, Cmyc, Ezh2, igf2)
hemangioma<-data.frame(hemangioma)

library(latticeExtra) #marginal plot
marginal.plot(hemangioma)
plot(hemangioma) #correlation
library(car)
par(mfrow=c(2,2))
scatterplot(hemangioma$age ~ hemangioma$RB, main = "Age vs RB") #scatterplot
scatterplot(hemangioma$age ~ hemangioma$p16, main = "Age vs p16")
scatterplot(hemangioma$age ~ hemangioma$DLK, main = "Age vs DLK")
scatterplot(hemangioma$age ~ hemangioma$Nanog, main = "Age vs Nanog")
scatterplot(hemangioma$age ~ hemangioma$Cmyc, main = "Age vs Cmyc")
scatterplot(hemangioma$age ~ hemangioma$Ezh2, main = "Age vs Ezh2")
scatterplot(hemangioma$age ~ hemangioma$igf2, main = "Age vs igf2")

#Factor Analysis
factanal(hemangioma, factors = 1)$loadings #significant result, 1 factor is not sufficient; The p-value is 0.0073 
factanal(hemangioma, factors = 2)$loadings #insignificant result, 2 factors are sufficient; The p-value is 0.0622 
factanal(hemangioma, factors = 3)$loadings #insignificant result, 3 factors are sufficient; The p-value is 0.33
factanal(hemangioma, factors = 4)$loadings #insignificant result, 4 factors are sufficient; The p-value is 0.519
# Can't go past 5 because 5 factors are too many for 8 variables. 

#Factor Analysis without outliers
library(sos); findFn("outlier") #outlier removal
hemangioma.outliersremoved<- hemangioma[!apply(sapply(hemangioma, function(x) abs(scale(x)) >= 2), 1, any), ]
factanal(hemangioma.outliersremoved, factors = 1)$loadings #significant result, 1 factor is not sufficient; The p-value is 0.00583
factanal(hemangioma.outliersremoved, factors = 2)$loadings #insignificant result, 2 factors are sufficient; The p-value is 0.586
factanal(hemangioma.outliersremoved, factors = 3)$loadings #insignificant result, 3 factors are sufficient; The p-value is 0.738 
factanal(hemangioma.outliersremoved, factors = 4)$loadings #insignificant result, 4 factors are sufficient; The p-value is 0.523

