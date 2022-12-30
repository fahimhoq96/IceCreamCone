Ice = read.table("icecreamcone.txt", header = TRUE)
summary(Ice)

#Scatter Plots
plot(Ice)

#histrograms

par(mfrow=c(2,2))
hist(Ice$viscosity, main = "", xlab = "Viscosity", col = "green")
hist(Ice$moisture, main = "", xlab = "moisture", col = "blue")
hist(Ice$protein, main = "", xlab = "protein", col = "blue")
hist(Ice$ash, main = "", xlab = "ash", col = "blue")

#boxplots
par(mfrow=c(1,3))
boxplot(Ice$viscosity, Ice$moisture, xlab = "viscosity vs moisture", col = c("green", "blue"))
boxplot(Ice$viscosity, Ice$protein, xlab = "viscosity vs protein", col = c("green", "blue"))
boxplot(Ice$viscosity, Ice$ash, xlab = "viscosity vs ash", col = c("green", "blue"))
boxplot(Ice$moisture, Ice$protein, Ice$ash, xlab = "moisture vs protein vs ash", col = c("blue","red","black"))

#Regression
m1 = lm(viscosity~moisture, data=Ice)
summary(m1)
m2 = lm(viscosity~protein, data=Ice)
summary(m2)
m3 = lm(viscosity~ash, data=Ice)
summary(m3)
m4 = lm(viscosity~moisture+protein+ash,data=Ice)
summary(m4)
m5 = lm(viscosity~moisture+protein+ash + I(moisture^2) + I(protein^2) + I(ash^2), data=Ice)
summary(m5)
m6 = lm(viscosity~ protein+ash + I(protein^2) + I(ash^2), data=Ice)
summary(m6)
anova(m6,m5)

#Interaction Terms
n1 =  lm(viscosity~moisture+protein+ash + I(moisture*protein)+I(protein*ash)+I(ash*moisture)+I(moisture*ash*protein) ,data=Ice)
summary(n1)
n2 =  lm(viscosity~moisture+protein+ash +I(ash*moisture) ,data=Ice)
summary(n2)
n3 = lm(viscosity~moisture+protein+ash +I(protein*ash) ,data=Ice)
summary(n3)
n4 = lm(viscosity~moisture+protein+ash +I(protein*moisture) ,data=Ice)
summary(n4)
n5= lm(viscosity~moisture+protein+ash +I(moisture*ash*protein) ,data=Ice)
summary(n5)


#Multicolinearity
library(car)
vif(m4)
vif(m5)
vif(m6)


#plots
par(mfrow=c(2,2))
plot(m4)
plot(m5)
plot(m6)
plot(m7)
plot(m8)
plot(m9)
plot(m10)



#Box_Cox Transformations

par(mfrow=c(1,1))
bc5 = boxCox(m6)

m7 = lm(as.numeric(viscosity)^0.25 ~ protein+ash + I(protein^2) + I(ash^2), data=Ice)
summary(m7)
m8 = lm(as.numeric(viscosity)^0.5 ~ protein+ash + I(protein^2) + I(ash^2), data=Ice)
summary(m8)
m9 = lm(as.numeric(viscosity)^0.75 ~ protein+ash + I(protein^2) + I(ash^2), data=Ice)
summary(m9)

m10 = lm(as.numeric(viscosity)^0.5 ~ protein+ash + I(ash^2), data=Ice)
summary(m10)

anova(m10,m8)

#AIC/BIC
AIC(m4,m5,m6)
BIC(m4,m5,m6)



#Added variable plots
par(mfrow=c(1,3))
avPlots(m1)
avPlots(m2)
avPlots(m3)


#Studentised Residuals
library(MASS)
par(mfrow=c(1,3))
plot(studres(m4),main="Model4 Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2)
plot(studres(m5),main="Model5 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 
plot(studres(m6),main="Model6 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 

plot(studres(m7),main="Model7 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 
plot(studres(m8),main="Model8 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 
plot(studres(m9),main="Model9 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 
plot(studres(m10),main="Model10 of Studentized Residuals",xlab="",ylab="Studentized Residuals",pch=16,cex=1.2,ylim=c( -2,2))  
abline(h=0,lty=2) 


#Others
library(car)
outlierTest(m10)


#K-Fold


k=3
n=dim(Ice)[1]
rand=sample(1:n)
Ice$group<-rep(0,dim(Ice)[1])
Ice$group[rand[1:(n/k)]]=1
Ice$group[rand[(n/k)+1:(2*n/k)]]=2
Ice$group[rand[(2*n/k+1):n]]=3
spr=rep(0,k)
msetr=rep(0,k)
for(i in 1:k)
{
  training=Ice[Ice$group!=i,c(1:3,4)]
  test=Ice[Ice$group==i,c(1:3,4)]
  temp=m6
  msetr[i]=deviance(temp)/df.residual(temp)
  spr[i]=sum((test$viscosity-predict(temp,data.frame(test[,c("protein","ash","moisture")]),type="response"))^2)/(n/k)
}
spr
msetr
mean(spr)
mean(msetr)
