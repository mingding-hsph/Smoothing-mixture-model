##Smoothing mixture model using data with poisson distribution

##The script is largely the same as normal distribution. 
##The main difference is model fit of GAMM4, calculation of pearson residual, and obtainment of predicted values

##This tutorial includes R script assuming 2-5 groups of trajectories 

##Parameters you can change
##1) exposure: we named it "age"
##2) outcome: we named it "bmi"
##3) time-varying covariates: these can be included in the "gamm4" model directly
##4) d_inter: loops of iteration of SMM
##5) g: number of groups

##Model output 
##1) Coefficient estimation
##2) Log likelihood (LLK): refer to the LLK file 
##3) Bayesian information criterion (BIC): : refer to the BIC file 
##4) Group classification: refer to the class file
##5) Plot trajectories 

library(MASS) ##generate data
library(nlme)
library(gmodels)  ##cross table
library(mgcv)  ##GAMM
library(gamm4) ##GAMM4

setwd("/udd/nhmid/mixture.model/revision.github")


######################################
############Two groups################
######################################

bmi_guts <- read.csv(file="simulated.poisson.csv", header=TRUE)

colnames(bmi_guts)=c("id" , "age", "bmi", "group_t")  

g=2          ##g is group member

d_inter=20    ##iterative loops   

LLK_2=data.frame(matrix(,nrow=d_inter, ncol=1))
BIC_2=data.frame(matrix(,nrow=d_inter, ncol=1))
class_2=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi1_2=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi2_2=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))

colnames(LLK_2)=cbind("LLK")
colnames(BIC_2)=cbind("BIC")
colnames(class_2)=cbind("id", "class")
colnames(p_bmi1_2)=cbind("id", "predicted value from group 1")
colnames(p_bmi2_2)=cbind("id", "predicted value from group 2")

##initial assignment of groups: we used mean BMI to categorize individuals into two groups

one=rep(1,nrow(bmi_guts))
bmi_guts=cbind(bmi_guts,one)

total_obs=aggregate(one ~ id, data=bmi_guts, FUN=sum)
colnames(total_obs)=c("id", "t_obs")

missing_obs=aggregate(bmi ~ id, data=bmi_guts, function(x) {sum(is.na(x))}, na.action = NULL)
colnames(missing_obs)=c("id", "m_obs")

data_obs=merge(total_obs, missing_obs, by="id")

bmi_total=aggregate(bmi ~ id, data=bmi_guts, na.rm=T, FUN=sum)
colnames(bmi_total)=c("id", "bmi_total")

bmi_two=merge(data_obs, bmi_total, by="id")
bmi_two$obs=bmi_two$t_obs-bmi_two$m_obs

bmi_two$bmi_mean=bmi_two$bmi_total/bmi_two$obs

bmi_two$group=cut(bmi_two$bmi_mean,breaks=c(quantile(bmi_two$bmi_mean, probs=seq(0,1, by=1/2))), labels=c("1","2"), include.lowest=TRUE) 

bmi_guts=merge(bmi_guts, bmi_two, by="id")

###E-M algorithm

for (i in 1:d_inter)

{

####M step

data1=bmi_guts[which(bmi_guts$group==1),]
data2=bmi_guts[which(bmi_guts$group==2),]

b1=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"),  data=data1)  

b2=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"),  data=data2)  


p_1=predict(b1$gam, se.fit=TRUE, newdata=bmi_guts)
p1=p_1$fit

p_2=predict(b2$gam, se.fit=TRUE, newdata=bmi_guts)
p2=p_2$fit

bmi_guts=cbind(bmi_guts,p1,p2)

####E step 
###calculate Pearson residual for each participant  

bmi_guts$y1_hat=exp(bmi_guts$p1)
bmi_guts$y2_hat=exp(bmi_guts$p2)

bmi_guts$var1_varying=((bmi_guts$y1_hat-bmi_guts$bmi)^2)/bmi_guts$y1_hat
bmi_guts$var2_varying=((bmi_guts$y2_hat-bmi_guts$bmi)^2)/bmi_guts$y2_hat

var1_sum=aggregate(var1_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var2_sum=aggregate(var2_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)

var_sum=merge(var1_sum, var2_sum,by="id")

colnames(var_sum)=c("id","var1_sum","var2_sum")

###assign participants to groups with the lowest residual
 
bmi_guts=merge(var_sum, bmi_guts, by="id",,all=TRUE)

bmi_guts$var_min=pmin(bmi_guts$var1_sum, bmi_guts$var2_sum)

bmi_guts$group[bmi_guts$var_min == bmi_guts$var1_sum] <- 1
bmi_guts$group[bmi_guts$var_min == bmi_guts$var2_sum] <- 2

##calculate LLK and BIC 

LLK_2[i,1]=logLik(b1$mer)+logLik(b2$mer)  

df=summary(b1$gam)$edf+sum(summary(b1$gam)$pTerms.df)+1+summary(b2$gam)$edf+sum(summary(b2$gam)$pTerms.df)+1

BIC_2[i,1]=-2*LLK_2[i,1]+log(nrow(bmi_guts))*(df+g-1)

##save the newly classified group for next loop 

bmi_guts <- subset(bmi_guts, select = c(id, age, bmi, group_t, group))

##if any of the groups has no participants, the E-M algorithm will stop, indicating that less number of groups should be assumed.  
##here we use criteria of less than one participants, as each participant has 20 observations

if (as.matrix(table(bmi_guts$group))[1,1]<=20 | as.matrix(table(bmi_guts$group))[2,1]<=20)
{
break
}

}

if (as.matrix(table(bmi_guts$group))[1,1]>20 | as.matrix(table(bmi_guts$group))[2,1]>20)

{
###model estimation

summary(b1)
summary(b2)

##LLK and BIC

LLK_2_final=max(LLK_2[d_inter,])
BIC_2_final=min(BIC_2[d_inter,])

write.csv(LLK_2,file="poisson_LLK_2.csv", na="", col.names =T, row.names=T)
write.csv(BIC_2,file="poisson_BIC_2.csv", na="", col.names =T, row.names=T)

##group classification

class_2[,1]=bmi_guts$id
class_2[,2]=bmi_guts$group

write.csv(class_2,file="poisson_class_2.csv", na="", col.names =T, row.names=T)

###plot figure

p_bmi1_2[,1]=bmi_guts$id
p_bmi2_2[,1]=bmi_guts$id

p_bmi1_2[,2]=predict(b1$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi2_2[,2]=predict(b2$gam,  se.fit=TRUE, newdata=bmi_guts)$fit

class_pbmi=data.frame(cbind(class_2[,2], p_bmi1_2[,2], p_bmi2_2[,2]))

colnames(class_pbmi)=cbind("class","pbmi1_pre","pbmi2_pre")

class_pbmi$pbmi1=exp(class_pbmi$pbmi1_pre)
class_pbmi$pbmi2=exp(class_pbmi$pbmi2_pre)

bmi_final=cbind(bmi_guts, class_pbmi)

bmi_final$pbmi_new <- ifelse(bmi_final$class==1, bmi_final$pbmi1, bmi_final$pbmi2)

final_plot_1=bmi_final[which(bmi_final$class==1),]
final_plot_2=bmi_final[which(bmi_final$class==2),]

out1_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_1, FUN=mean)
out2_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_2, FUN=mean)

jpeg('SMM_Poisson_2.png', width = 7, height = 7, units = 'in', res = 600)

par(mfrow=c(1,1))

plot(main="SMM_Poisson_2", out1_mean_avg$age, out1_mean_avg$pbmi_new, type="l", col="black", 
xlab="Time", ylab="Mean Y",ylim=range(0, 32),  xlim=range(0, 1))

lines(out2_mean_avg$age,out2_mean_avg$pbmi_new, type="l", col="red")

legend(0, 30, legend=c("Group 1", "Group 2"),
       col=c("black", "red"), lty=1:1, cex=0.9)

dev.off()

}

######################################
############Three groups##############
######################################

bmi_guts <- read.csv(file="simulated.poisson.csv", header=TRUE)

colnames(bmi_guts)=c("id" , "age", "bmi", "group_t")  

g=3          ##g is group member
d_inter=20    ##iterative loops  

LLK_3=data.frame(matrix(,nrow=d_inter, ncol=1))
BIC_3=data.frame(matrix(,nrow=d_inter, ncol=1))
class_3=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi1_3=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi2_3=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi3_3=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))


colnames(LLK_3)=cbind("LLK")
colnames(BIC_3)=cbind("BIC")
colnames(class_3)=cbind("id", "class")
colnames(p_bmi1_3)=cbind("id", "predicted value from group 1")
colnames(p_bmi2_3)=cbind("id", "predicted value from group 2")
colnames(p_bmi3_3)=cbind("id", "predicted value from group 3")

##initial assignment of groups: we used mean BMI to categorize individuals into three groups

one=rep(1,nrow(bmi_guts))
bmi_guts=cbind(bmi_guts,one)

total_obs=aggregate(one ~ id, data=bmi_guts, FUN=sum)
colnames(total_obs)=c("id", "t_obs")

missing_obs=aggregate(bmi ~ id, data=bmi_guts, function(x) {sum(is.na(x))}, na.action = NULL)
colnames(missing_obs)=c("id", "m_obs")

data_obs=merge(total_obs, missing_obs, by="id")

bmi_total=aggregate(bmi ~ id, data=bmi_guts, na.rm=T, FUN=sum)
colnames(bmi_total)=c("id", "bmi_total")

bmi_two=merge(data_obs, bmi_total, by="id")
bmi_two$obs=bmi_two$t_obs-bmi_two$m_obs

bmi_two$bmi_mean=bmi_two$bmi_total/bmi_two$obs

bmi_two$group=cut(bmi_two$bmi_mean,breaks=c(quantile(bmi_two$bmi_mean, probs=seq(0,1, by=1/3))), labels=c("1","2","3"), include.lowest=TRUE) 

bmi_guts=merge(bmi_guts, bmi_two, by="id")

###E-M algorithm

for (i in 1:d_inter)

{

##M step

data1=bmi_guts[which(bmi_guts$group==1),]
data2=bmi_guts[which(bmi_guts$group==2),]
data3=bmi_guts[which(bmi_guts$group==3),]

b1=gamm4(bmi ~ s(age),na.action=na.omit,   random = ~ (1|id), family=poisson(link = "log"), data=data1)  

b2=gamm4(bmi ~ s(age),na.action=na.omit,   random = ~ (1|id), family=poisson(link = "log"), data=data2)  

b3=gamm4(bmi ~ s(age),na.action=na.omit,   random = ~ (1|id), family=poisson(link = "log"), data=data3)  


p_1=predict(b1$gam, se.fit=TRUE, newdata=bmi_guts)
p1=p_1$fit

p_2=predict(b2$gam, se.fit=TRUE, newdata=bmi_guts)
p2=p_2$fit

p_3=predict(b3$gam, se.fit=TRUE, newdata=bmi_guts)
p3=p_3$fit


bmi_guts=cbind(bmi_guts,p1,p2,p3)


####E step 
###calculate Pearson residual for each participant  

bmi_guts$y1_hat=exp(bmi_guts$p1)
bmi_guts$y2_hat=exp(bmi_guts$p2)
bmi_guts$y3_hat=exp(bmi_guts$p3)

bmi_guts$var1_varying=((bmi_guts$y1_hat-bmi_guts$bmi)^2)/bmi_guts$y1_hat
bmi_guts$var2_varying=((bmi_guts$y2_hat-bmi_guts$bmi)^2)/bmi_guts$y2_hat
bmi_guts$var3_varying=((bmi_guts$y3_hat-bmi_guts$bmi)^2)/bmi_guts$y3_hat

var1_sum=aggregate(var1_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var2_sum=aggregate(var2_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var3_sum=aggregate(var3_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)

var12_sum=merge(var1_sum, var2_sum,by="id")
var_sum=merge(var12_sum, var3_sum, by="id")

colnames(var_sum)=c("id","var1_sum","var2_sum","var3_sum")

###assign participants to groups with the lowest residual

bmi_guts=merge(var_sum, bmi_guts, by="id",,all=TRUE)

bmi_guts$var_min=pmin(bmi_guts$var1_sum, bmi_guts$var2_sum, bmi_guts$var3_sum)

bmi_guts$group[bmi_guts$var_min == bmi_guts$var1_sum] <- 1
bmi_guts$group[bmi_guts$var_min == bmi_guts$var2_sum] <- 2
bmi_guts$group[bmi_guts$var_min == bmi_guts$var3_sum] <- 3

##calculate LLK and BIC 

LLK_3[i,1]=logLik(b1$mer)+logLik(b2$mer)+logLik(b3$mer)    

df=summary(b1$gam)$edf+sum(summary(b1$gam)$pTerms.df)+1+summary(b2$gam)$edf+sum(summary(b2$gam)$pTerms.df)+1+summary(b3$gam)$edf+sum(summary(b3$gam)$pTerms.df)+1

BIC_3[i,1]=-2*LLK_3[i,1]+log(nrow(bmi_guts))*(df+g-1)

##save the newly classified group for next loop 

bmi_guts <- subset(bmi_guts, select = c(id, age, bmi, group_t, group))

##if any of the groups has no participants, the E-M algorithm will stop, indicating that less number of groups should be assumed.  
##here we use criteria of less than one participants, as each participant has 20 observations

if (as.matrix(table(bmi_guts$group))[1,1]<=20 | as.matrix(table(bmi_guts$group))[2,1]<=20| as.matrix(table(bmi_guts$group))[3,1]<=20)
{
break
}

}

if (as.matrix(table(bmi_guts$group))[1,1]>20 | as.matrix(table(bmi_guts$group))[2,1]>20| as.matrix(table(bmi_guts$group))[3,1]>20)

{
###model estimation

summary(b1)
summary(b2)
summary(b3)

##LLK and BIC

LLK_3_final=max(LLK_3[d_inter,])
BIC_3_final=min(BIC_3[d_inter,])
BIC_3_final=min(BIC_3[d_inter,])

write.csv(LLK_3,file="poisson_LLK_3.csv", na="", col.names =T, row.names=T)
write.csv(BIC_3,file="poisson_BIC_3.csv", na="", col.names =T, row.names=T)

##group classification

class_3[,1]=bmi_guts$id
class_3[,2]=bmi_guts$group

write.csv(class_3,file="poisson_class_3.csv", na="", col.names =T, row.names=T)

###plot figure

p_bmi1_3[,1]=bmi_guts$id
p_bmi2_3[,1]=bmi_guts$id
p_bmi3_3[,1]=bmi_guts$id

p_bmi1_3[,2]=predict(b1$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi2_3[,2]=predict(b2$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi3_3[,2]=predict(b3$gam,  se.fit=TRUE, newdata=bmi_guts)$fit

class_pbmi=data.frame(cbind(class_3[,2], p_bmi1_3[,2], p_bmi2_3[,2], p_bmi3_3[,2]))

colnames(class_pbmi)=cbind("class","pbmi1_pre","pbmi2_pre","pbmi3_pre")

class_pbmi$pbmi1=exp(class_pbmi$pbmi1_pre)
class_pbmi$pbmi2=exp(class_pbmi$pbmi2_pre)
class_pbmi$pbmi3=exp(class_pbmi$pbmi3_pre)

bmi_final=cbind(bmi_guts, class_pbmi)

bmi_final$pbmi_new <- ifelse(bmi_final$class==1, bmi_final$pbmi1, ifelse(bmi_final$class==2,bmi_final$pbmi2, bmi_final$pbmi3))

final_plot_1=bmi_final[which(bmi_final$class==1),]
final_plot_2=bmi_final[which(bmi_final$class==2),]
final_plot_3=bmi_final[which(bmi_final$class==3),]

out1_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_1, FUN=mean)
out2_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_2, FUN=mean)
out3_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_3, FUN=mean)

jpeg('SMM_Poisson_3.png', width = 7, height = 7, units = 'in', res = 600)

par(mfrow=c(1,1))

plot(main="SMM_Poisson_3", out1_mean_avg$age, out1_mean_avg$pbmi_new, type="l", col="black", 
xlab="Time", ylab="Mean Y",ylim=range(0, 32),  xlim=range(0, 1))


lines(out2_mean_avg$age,out2_mean_avg$pbmi_new, type="l", col="red")
lines(out3_mean_avg$age,out3_mean_avg$pbmi_new, type="l", col="green")

legend(0, 30, legend=c("Group 1", "Group 2", "Group 3"),
       col=c("black", "red", "green"),  lty=1:1:1, cex=0.9)

dev.off()

}

######################################
############Four groups##############
######################################

bmi_guts <- read.csv(file="simulated.poisson.csv", header=TRUE)

colnames(bmi_guts)=c("id" , "age", "bmi", "group_t")  

g=4          ##g is group member
d_inter=20    ##iterative loops 


LLK_4=data.frame(matrix(,nrow=d_inter, ncol=1))
BIC_4=data.frame(matrix(,nrow=d_inter, ncol=1))
class_4=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi1_4=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi2_4=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi3_4=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi4_4=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))

colnames(LLK_4)=cbind("LLK")
colnames(BIC_4)=cbind("BIC")
colnames(class_4)=cbind("id", "class")
colnames(p_bmi1_4)=cbind("id", "predicted value from group 1")
colnames(p_bmi2_4)=cbind("id", "predicted value from group 2")
colnames(p_bmi3_4)=cbind("id", "predicted value from group 3")
colnames(p_bmi4_4)=cbind("id", "predicted value from group 4")

##initial assignment of groups: we used mean BMI to categorize individuals into three groups

one=rep(1,nrow(bmi_guts))
bmi_guts=cbind(bmi_guts,one)

total_obs=aggregate(one ~ id, data=bmi_guts, FUN=sum)
colnames(total_obs)=c("id", "t_obs")

missing_obs=aggregate(bmi ~ id, data=bmi_guts, function(x) {sum(is.na(x))}, na.action = NULL)
colnames(missing_obs)=c("id", "m_obs")

data_obs=merge(total_obs, missing_obs, by="id")

bmi_total=aggregate(bmi ~ id, data=bmi_guts, na.rm=T, FUN=sum)
colnames(bmi_total)=c("id", "bmi_total")

bmi_two=merge(data_obs, bmi_total, by="id")
bmi_two$obs=bmi_two$t_obs-bmi_two$m_obs

bmi_two$bmi_mean=bmi_two$bmi_total/bmi_two$obs

bmi_two$group=cut(bmi_two$bmi_mean,breaks=c(quantile(bmi_two$bmi_mean, probs=seq(0,1, by=1/4))), labels=c("1","2","3","4"), include.lowest=TRUE) 

bmi_guts=merge(bmi_guts, bmi_two, by="id")

###E-M algorithm

for (i in 1:d_inter)

{

##M step

data1=bmi_guts[which(bmi_guts$group==1),]
data2=bmi_guts[which(bmi_guts$group==2),]
data3=bmi_guts[which(bmi_guts$group==3),]
data4=bmi_guts[which(bmi_guts$group==4),]

b1=gamm4(bmi ~ s(age),na.action=na.omit, random = ~ (1|id), family=poisson(link = "log"), data=data1)  

b2=gamm4(bmi ~ s(age),na.action=na.omit, random = ~ (1|id), family=poisson(link = "log"), data=data2)  

b3=gamm4(bmi ~ s(age),na.action=na.omit, random = ~ (1|id), family=poisson(link = "log"), data=data3)  

b4=gamm4(bmi ~ s(age),na.action=na.omit, random = ~ (1|id), family=poisson(link = "log"), data=data4)  


p_1=predict(b1$gam, se.fit=TRUE, newdata=bmi_guts)
p1=p_1$fit

p_2=predict(b2$gam, se.fit=TRUE, newdata=bmi_guts)
p2=p_2$fit

p_3=predict(b3$gam, se.fit=TRUE, newdata=bmi_guts)
p3=p_3$fit

p_4=predict(b4$gam, se.fit=TRUE, newdata=bmi_guts)
p4=p_4$fit


bmi_guts=cbind(bmi_guts,p1,p2,p3,p4)


####E step 
###calculate Pearson residual for each participant  

bmi_guts$y1_hat=exp(bmi_guts$p1)
bmi_guts$y2_hat=exp(bmi_guts$p2)
bmi_guts$y3_hat=exp(bmi_guts$p3)
bmi_guts$y4_hat=exp(bmi_guts$p4)

bmi_guts$var1_varying=((bmi_guts$y1_hat-bmi_guts$bmi)^2)/bmi_guts$y1_hat
bmi_guts$var2_varying=((bmi_guts$y2_hat-bmi_guts$bmi)^2)/bmi_guts$y2_hat
bmi_guts$var3_varying=((bmi_guts$y3_hat-bmi_guts$bmi)^2)/bmi_guts$y3_hat
bmi_guts$var4_varying=((bmi_guts$y4_hat-bmi_guts$bmi)^2)/bmi_guts$y4_hat

var1_sum=aggregate(var1_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var2_sum=aggregate(var2_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var3_sum=aggregate(var3_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var4_sum=aggregate(var4_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)

var12_sum=merge(var1_sum, var2_sum,by="id")
var123_sum=merge(var12_sum, var3_sum, by="id")
var_sum=merge(var123_sum, var4_sum, by="id")

colnames(var_sum)=c("id","var1_sum","var2_sum","var3_sum","var4_sum")

###assign participants to groups with the lowest residual

bmi_guts=merge(var_sum, bmi_guts, by="id",,all=TRUE)

bmi_guts$var_min=pmin(bmi_guts$var1_sum, bmi_guts$var2_sum, bmi_guts$var3_sum,bmi_guts$var4_sum)

bmi_guts$group[bmi_guts$var_min == bmi_guts$var1_sum] <- 1
bmi_guts$group[bmi_guts$var_min == bmi_guts$var2_sum] <- 2
bmi_guts$group[bmi_guts$var_min == bmi_guts$var3_sum] <- 3
bmi_guts$group[bmi_guts$var_min == bmi_guts$var4_sum] <- 4

##calculate LLK and BIC 

LLK_4[i,1]=logLik(b1$mer)+logLik(b2$mer)+logLik(b3$mer)+logLik(b4$mer)    

df=summary(b1$gam)$edf+sum(summary(b1$gam)$pTerms.df)+1+summary(b2$gam)$edf+sum(summary(b2$gam)$pTerms.df)+1+summary(b3$gam)$edf+sum(summary(b3$gam)$pTerms.df)+1+summary(b4$gam)$edf+sum(summary(b4$gam)$pTerms.df)+1 

BIC_4[i,1]=-2*LLK_4[i,1]+log(nrow(bmi_guts))*(df+g-1)

##save the newly classified group for next loop 

bmi_guts <- subset(bmi_guts, select = c(id, age, bmi, group_t, group))

##if any of the groups has no participants, the E-M algorithm will stop, indicating that less number of groups should be assumed.  
##here we use criteria of less than one participants, as each participant has 20 observations

if (as.matrix(table(bmi_guts$group))[1,1]<=20 | as.matrix(table(bmi_guts$group))[2,1]<=20| as.matrix(table(bmi_guts$group))[3,1]<=20|as.matrix(table(bmi_guts$group))[4,1]<=20)
{
break
}

}

if (as.matrix(table(bmi_guts$group))[1,1]>20 | as.matrix(table(bmi_guts$group))[2,1]>20| as.matrix(table(bmi_guts$group))[3,1]>20|as.matrix(table(bmi_guts$group))[4,1]>20)

{  
###model estimation

summary(b1)
summary(b2)
summary(b3)
summary(b4)

##LLK and BIC

LLK_4_final=max(LLK_4[d_inter,])
BIC_4_final=min(BIC_4[d_inter,])

write.csv(LLK_4,file="poisson_LLK_4.csv", na="", col.names =T, row.names=T)
write.csv(BIC_4,file="poisson_BIC_4.csv", na="", col.names =T, row.names=T)

##group classification

class_4[,1]=bmi_guts$id
class_4[,2]=bmi_guts$group

write.csv(class_4,file="poisson_class_4.csv", na="", col.names =T, row.names=T)

###plot figure

p_bmi1_4[,1]=bmi_guts$id
p_bmi2_4[,1]=bmi_guts$id
p_bmi3_4[,1]=bmi_guts$id
p_bmi4_4[,1]=bmi_guts$id

p_bmi1_4[,2]=predict(b1$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi2_4[,2]=predict(b2$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi3_4[,2]=predict(b3$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi4_4[,2]=predict(b4$gam,  se.fit=TRUE, newdata=bmi_guts)$fit

class_pbmi=data.frame(cbind(class_4[,2], p_bmi1_4[,2], p_bmi2_4[,2], p_bmi3_4[,2],p_bmi4_4[,2]))

colnames(class_pbmi)=cbind("class","pbmi1_pre","pbmi2_pre","pbmi3_pre","pbmi4_pre")

class_pbmi$pbmi1=exp(class_pbmi$pbmi1_pre)
class_pbmi$pbmi2=exp(class_pbmi$pbmi2_pre)
class_pbmi$pbmi3=exp(class_pbmi$pbmi3_pre)
class_pbmi$pbmi4=exp(class_pbmi$pbmi4_pre)

bmi_final=cbind(bmi_guts, class_pbmi)

bmi_final$pbmi_new <- ifelse(bmi_final$class==1, bmi_final$pbmi1, ifelse(bmi_final$class==2,bmi_final$pbmi2, ifelse(bmi_final$class==3,bmi_final$pbmi3, bmi_final$pbmi4)))

final_plot_1=bmi_final[which(bmi_final$class==1),]
final_plot_2=bmi_final[which(bmi_final$class==2),]
final_plot_3=bmi_final[which(bmi_final$class==3),]
final_plot_4=bmi_final[which(bmi_final$class==4),]

out1_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_1, FUN=mean)
out2_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_2, FUN=mean)
out3_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_3, FUN=mean)
out4_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_4, FUN=mean)

jpeg('SMM_poisson_4.png', width = 7, height = 7, units = 'in', res = 600)

par(mfrow=c(1,1))

plot(main="SMM_poisson_4", out1_mean_avg$age, out1_mean_avg$pbmi_new, type="l", col="black", 
xlab="Time", ylab="Mean Y", ylim=range(0, 1.3),  xlim=range(0, 1))

lines(out2_mean_avg$age,out2_mean_avg$pbmi_new, type="l", col="red")
lines(out3_mean_avg$age,out3_mean_avg$pbmi_new, type="l", col="green")
lines(out4_mean_avg$age,out4_mean_avg$pbmi_new, type="l", col="pink")

legend(0, 30, legend=c("Group 1", "Group 2", "Group 3", "Group 4"),
       col=c("black", "red", "green", "pink"),  lty=1:1:1:1, cex=0.9)

dev.off()

}

######################################
############Five groups##############
######################################

bmi_guts <- read.csv(file="simulated.poisson.csv", header=TRUE)

colnames(bmi_guts)=c("id" , "age", "bmi", "group_t")  

g=5          ##g is group member
d_inter=20    ##iterative loops 

LLK_5=data.frame(matrix(,nrow=d_inter, ncol=1))
BIC_5=data.frame(matrix(,nrow=d_inter, ncol=1))
class_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi1_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi2_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi3_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi4_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))
p_bmi5_5=data.frame(matrix(,nrow=nrow(bmi_guts), ncol=2))

colnames(LLK_5)=cbind("LLK")
colnames(BIC_5)=cbind("BIC")
colnames(class_5)=cbind("id", "class")
colnames(p_bmi1_5)=cbind("id", "predicted value from group 1")
colnames(p_bmi2_5)=cbind("id", "predicted value from group 2")
colnames(p_bmi3_5)=cbind("id", "predicted value from group 3")
colnames(p_bmi4_5)=cbind("id", "predicted value from group 4")
colnames(p_bmi5_5)=cbind("id", "predicted value from group 5")

##initial assignment of groups: we used mean BMI to categorize individuals into three groups

one=rep(1,nrow(bmi_guts))
bmi_guts=cbind(bmi_guts,one)

total_obs=aggregate(one ~ id, data=bmi_guts, FUN=sum)
colnames(total_obs)=c("id", "t_obs")

missing_obs=aggregate(bmi ~ id, data=bmi_guts, function(x) {sum(is.na(x))}, na.action = NULL)
colnames(missing_obs)=c("id", "m_obs")

data_obs=merge(total_obs, missing_obs, by="id")

bmi_total=aggregate(bmi ~ id, data=bmi_guts, na.rm=T, FUN=sum)
colnames(bmi_total)=c("id", "bmi_total")

bmi_two=merge(data_obs, bmi_total, by="id")
bmi_two$obs=bmi_two$t_obs-bmi_two$m_obs

bmi_two$bmi_mean=bmi_two$bmi_total/bmi_two$obs

bmi_two$group=cut(bmi_two$bmi_mean,breaks=c(quantile(bmi_two$bmi_mean, probs=seq(0,1, by=1/5))), labels=c("1","2","3","4","5"), include.lowest=TRUE) 

bmi_guts=merge(bmi_guts, bmi_two, by="id")

###E-M algorithm

for (i in 1:d_inter)

{

data1=bmi_guts[which(bmi_guts$group==1),]
data2=bmi_guts[which(bmi_guts$group==2),]
data3=bmi_guts[which(bmi_guts$group==3),]
data4=bmi_guts[which(bmi_guts$group==4),]
data5=bmi_guts[which(bmi_guts$group==5),]


##M step

b1=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"), data=data1)  

b2=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"), data=data2)  

b3=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"), data=data3)  

b4=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"), data=data4)  

b5=gamm4(bmi ~ s(age),na.action=na.omit,  random = ~ (1|id), family=poisson(link = "log"), data=data5)  


p_1=predict(b1$gam, se.fit=TRUE, newdata=bmi_guts)
p1=p_1$fit

p_2=predict(b2$gam, se.fit=TRUE, newdata=bmi_guts)
p2=p_2$fit

p_3=predict(b3$gam, se.fit=TRUE, newdata=bmi_guts)
p3=p_3$fit

p_4=predict(b4$gam, se.fit=TRUE, newdata=bmi_guts)
p4=p_4$fit

p_5=predict(b5$gam, se.fit=TRUE, newdata=bmi_guts)
p5=p_5$fit

bmi_guts=cbind(bmi_guts,p1,p2,p3,p4,p5)


####E step 
###calculate Pearson residual for each participant  

bmi_guts$y1_hat=exp(bmi_guts$p1)
bmi_guts$y2_hat=exp(bmi_guts$p2)
bmi_guts$y3_hat=exp(bmi_guts$p3)
bmi_guts$y4_hat=exp(bmi_guts$p4)
bmi_guts$y5_hat=exp(bmi_guts$p5)

bmi_guts$var1_varying=((bmi_guts$y1_hat-bmi_guts$bmi)^2)/bmi_guts$y1_hat
bmi_guts$var2_varying=((bmi_guts$y2_hat-bmi_guts$bmi)^2)/bmi_guts$y2_hat
bmi_guts$var3_varying=((bmi_guts$y3_hat-bmi_guts$bmi)^2)/bmi_guts$y3_hat
bmi_guts$var4_varying=((bmi_guts$y4_hat-bmi_guts$bmi)^2)/bmi_guts$y4_hat
bmi_guts$var5_varying=((bmi_guts$y5_hat-bmi_guts$bmi)^2)/bmi_guts$y5_hat

var1_sum=aggregate(var1_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var2_sum=aggregate(var2_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var3_sum=aggregate(var3_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var4_sum=aggregate(var4_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)
var5_sum=aggregate(var5_varying ~ id, data=bmi_guts, na.rm=T, FUN=sum)

var12_sum=merge(var1_sum, var2_sum,by="id")
var123_sum=merge(var12_sum, var3_sum, by="id")
var1234_sum=merge(var123_sum, var4_sum, by="id")

var_sum=merge(var1234_sum, var5_sum, by="id")

colnames(var_sum)=c("id","var1_sum","var2_sum","var3_sum","var4_sum","var5_sum")

###assign participants to groups with the lowest residual
 
bmi_guts=merge(var_sum, bmi_guts, by="id", all=TRUE)

bmi_guts$var_min=pmin(bmi_guts$var1_sum, bmi_guts$var2_sum, bmi_guts$var3_sum,bmi_guts$var4_sum,bmi_guts$var5_sum)

bmi_guts$group[bmi_guts$var_min == bmi_guts$var1_sum] <- 1
bmi_guts$group[bmi_guts$var_min == bmi_guts$var2_sum] <- 2
bmi_guts$group[bmi_guts$var_min == bmi_guts$var3_sum] <- 3
bmi_guts$group[bmi_guts$var_min == bmi_guts$var4_sum] <- 4
bmi_guts$group[bmi_guts$var_min == bmi_guts$var5_sum] <- 5

##calculate LLK and BIC 

LLK_5[i,1]=logLik(b1$mer)+logLik(b2$mer)+logLik(b3$mer)+logLik(b4$mer)+logLik(b5$mer)      

df=summary(b1$gam)$edf+sum(summary(b1$gam)$pTerms.df)+1+summary(b2$gam)$edf+sum(summary(b2$gam)$pTerms.df)+1+summary(b3$gam)$edf+sum(summary(b3$gam)$pTerms.df)+1+summary(b4$gam)$edf+sum(summary(b4$gam)$pTerms.df)+1 +summary(b5$gam)$edf+sum(summary(b5$gam)$pTerms.df)+1 

BIC_5[i,1]=-2*LLK_5[i,1]+log(nrow(bmi_guts))*(df+g-1)

##save the newly classified group for next loop 

bmi_guts <- subset(bmi_guts, select = c(id, age, bmi, group_t, group))

##if any of the groups has no participants, the E-M algorithm will stop, indicating that less number of groups should be assumed.  
##here we use criteria of less than one participants, as each participant has 20 observations

if (as.matrix(table(bmi_guts$group))[1,1]<=20 | as.matrix(table(bmi_guts$group))[2,1]<=20| as.matrix(table(bmi_guts$group))[3,1]<=20|as.matrix(table(bmi_guts$group))[4,1]<=20|as.matrix(table(bmi_guts$group))[5,1]<=20)
{
break
}

}

if (as.matrix(table(bmi_guts$group))[1,1]>20 | as.matrix(table(bmi_guts$group))[2,1]>20| as.matrix(table(bmi_guts$group))[3,1]>20|as.matrix(table(bmi_guts$group))[4,1]>20|as.matrix(table(bmi_guts$group))[5,1]>20)

{

###model estimation

summary(b1)
summary(b2)
summary(b3)
summary(b4)
summary(b5)

##LLK and BIC

LLK_5_final=max(LLK_5[d_inter,])
BIC_5_final=min(BIC_5[d_inter,])

write.csv(LLK_5,file="poisson_LLK_5.csv", na="", col.names =T, row.names=T)
write.csv(BIC_5,file="poisson_BIC_5.csv", na="", col.names =T, row.names=T)

##group classification

class_5[,1]=bmi_guts$id
class_5[,2]=bmi_guts$group

write.csv(class_5,file="poisson_class_5.csv", na="", col.names =T, row.names=T)

###plot figure

p_bmi1_5[,1]=bmi_guts$id
p_bmi2_5[,1]=bmi_guts$id
p_bmi3_5[,1]=bmi_guts$id
p_bmi4_5[,1]=bmi_guts$id
p_bmi5_5[,1]=bmi_guts$id

p_bmi1_5[,2]=predict(b1$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi2_5[,2]=predict(b2$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi3_5[,2]=predict(b3$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi4_5[,2]=predict(b4$gam,  se.fit=TRUE, newdata=bmi_guts)$fit
p_bmi5_5[,2]=predict(b5$gam,  se.fit=TRUE, newdata=bmi_guts)$fit

class_pbmi=data.frame(cbind(class_5[,2], p_bmi1_5[,2], p_bmi2_5[,2], p_bmi3_5[,2],p_bmi4_5[,2], p_bmi5_5[,2]))

colnames(class_pbmi)=cbind("class","pbmi1_pre","pbmi2_pre","pbmi3_pre","pbmi4_pre", "pbmi5_pre")

class_pbmi$pbmi1=exp(class_pbmi$pbmi1_pre)
class_pbmi$pbmi2=exp(class_pbmi$pbmi2_pre)
class_pbmi$pbmi3=exp(class_pbmi$pbmi3_pre)
class_pbmi$pbmi4=exp(class_pbmi$pbmi4_pre)
class_pbmi$pbmi5=exp(class_pbmi$pbmi5_pre)


bmi_final=cbind(bmi_guts, class_pbmi)

bmi_final$pbmi_new <- ifelse(bmi_final$class==1, bmi_final$pbmi1, ifelse(bmi_final$class==2,bmi_final$pbmi2, ifelse(bmi_final$class==3,bmi_final$pbmi3, ifelse(bmi_final$class==4,bmi_final$pbmi4, bmi_final$pbmi5))))

final_plot_1=bmi_final[which(bmi_final$class==1),]
final_plot_2=bmi_final[which(bmi_final$class==2),]
final_plot_3=bmi_final[which(bmi_final$class==3),]
final_plot_4=bmi_final[which(bmi_final$class==4),]
final_plot_5=bmi_final[which(bmi_final$class==5),]

out1_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_1, FUN=mean)
out2_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_2, FUN=mean)
out3_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_3, FUN=mean)
out4_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_4, FUN=mean)
out5_mean_avg=aggregate(pbmi_new ~ age, data=final_plot_5, FUN=mean)

jpeg('SMM_poisson_5.png', width = 7, height = 7, units = 'in', res = 600)

par(mfrow=c(1,1))

plot(main="SMM_poisson_5", out1_mean_avg$age,out1_mean_avg$pbmi_new, type="l", col="black", lty=1,
xlab="Time", ylab="Mean Y",ylim=range(0, 32),  xlim=range(0, 1))

lines(out2_mean_avg$age,out2_mean_avg$pbmi_new, type="l", col="red")
lines(out3_mean_avg$age,out3_mean_avg$pbmi_new, type="l", col="green")
lines(out4_mean_avg$age,out4_mean_avg$pbmi_new, type="l", col="pink")
lines(out5_mean_avg$age,out5_mean_avg$pbmi_new, type="l", col="blue")

legend(0, 30, legend=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"),
       col=c("black", "red", "green", "pink", "blue"),  lty=1:1:1:1:1, cex=0.9)

dev.off()

}
