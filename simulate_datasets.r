
library(dplyr) ##for binding rows 

setwd("/udd/nhmid/mixture.model/revision.github")

######################################################
################normal distribution###################
######################################################


for (sn in 1:1000)  ##simulate 1000 datasets

{

set.seed(sn)

num=100

r=20

b_r=2    ##random effect for each participant which indicates high, medium, or low separation of trajectories

t <- as.matrix(rep(1:r))/r

smooth_1 <- function(x) 0.2*x^11*(10*(1 - x))^6+10*(10*x)^3*(1 - x)^10
smooth_2 <- function(x) -1*x^11*(10*(1 - x))^6-1*(10*x)^3*(1 - x)^10
smooth_3 <- function(x) 1*x^11*(10*(1 - x))^6+18*(10*x)^3*(1 - x)^10
smooth_4 <- function(x) 2.5*x^11*(10*(1 - x))^6-10*(10*x)^3*(1 - x)^10

y1=matrix(,nrow=num*r, ncol=3)
y2=matrix(,nrow=num*r, ncol=3)
y3=matrix(,nrow=num*r, ncol=3)
y4=matrix(,nrow=num*r, ncol=3)

eij <- as.matrix(rnorm(num*r*4,0,1))   ##random effect for each observation for normal distribution 

b1=as.matrix(rnorm(num,0,b_r))
b2=as.matrix(rnorm(num,0,b_r))
b3=as.matrix(rnorm(num,0,b_r))
b4=as.matrix(rnorm(num,0,b_r))

for (i in 1:num)
{
for (j in 1:r)
{
y1[r*(i-1)+j,1]=i
y2[r*(i-1)+j,1]=i+num
y3[r*(i-1)+j,1]=i+num*2
y4[r*(i-1)+j,1]=i+num*3

y1[r*(i-1)+j,3]=j/r
y2[r*(i-1)+j,3]=j/r
y3[r*(i-1)+j,3]=j/r
y4[r*(i-1)+j,3]=j/r

y1[r*(i-1)+j,2]=smooth_1(t[j,]) +b1[i,1]+eij[r*(i-1)+j,]
y2[r*(i-1)+j,2]=smooth_2(t[j,]) +b2[i,1]+eij[2*(r*(i-1)+j),]
y3[r*(i-1)+j,2]=smooth_3(t[j,]) +b3[i,1]+eij[3*(r*(i-1)+j),]
y4[r*(i-1)+j,2]=smooth_4(t[j,]) +b4[i,1]+eij[4*(r*(i-1)+j),]

}

}


dat_1=data.frame(y1)
dat_2=data.frame(y2)
dat_3=data.frame(y3)
dat_4=data.frame(y4)

colnames(dat_1)=cbind("id", "yij","obstime")
colnames(dat_2)=cbind("id", "yij","obstime")
colnames(dat_3)=cbind("id", "yij","obstime")
colnames(dat_4)=cbind("id", "yij","obstime")

dat_1$group=rep(1,nrow(dat_1))
dat_2$group=rep(2,nrow(dat_2))
dat_3$group=rep(3,nrow(dat_3))
dat_4$group=rep(4,nrow(dat_4))

###output datasets

data_normal_all =data.frame(bind_rows(dat_1,dat_2,dat_3,dat_4))

data_normal <- subset(data_normal_all, select = c(id, obstime, yij, group))

##write.csv(data_normal, paste0("normal_medium", sn, ".csv"),  na="", col.names =T, row.names=F)

}


######################################################
################Bernoulli distribution################
######################################################

for (sn in 1:1000)  ##simulate 1000 datasets

{

set.seed(sn)


num=100

r=20

t <- as.matrix(rep(1:r))/r

b_r=0.5  ##random effect for each participant which indicates high, medium, or low separation of trajectories

smooth_1 <- function(x) 0.3*x^11*(10*(1 - x))^6+0.1*(10*x)^3*(1 - x)^10-2
smooth_2 <- function(x) 0.1*x^11*(10*(1 - x))^6+10*(10*x)^3*(1 - x)^10-3
smooth_3 <- function(x) 0*x^11*(10*(1 - x))^6+0*(10*x)^(-1)*(1 - x)^0
smooth_4 <- function(x) -10*x^2*(10*(1 - x))^0+0*(10*x)^3*(1 - x)^10+2

b1=as.matrix(rnorm(num,0,b_r))
b2=as.matrix(rnorm(num,0,b_r))
b3=as.matrix(rnorm(num,0,b_r))
b4=as.matrix(rnorm(num,0,b_r))

y1=matrix(,nrow=num*r, ncol=3)
y2=matrix(,nrow=num*r, ncol=3)
y3=matrix(,nrow=num*r, ncol=3)
y4=matrix(,nrow=num*r, ncol=3)


for (i in 1:num)
{
for (j in 1:r)
{

y1[r*(i-1)+j,1]=i
y2[r*(i-1)+j,1]=i+num
y3[r*(i-1)+j,1]=i+num*2
y4[r*(i-1)+j,1]=i+num*3

y1[r*(i-1)+j,3]=j/r
y2[r*(i-1)+j,3]=j/r
y3[r*(i-1)+j,3]=j/r
y4[r*(i-1)+j,3]=j/r

y1[r*(i-1)+j,2]=smooth_1(t[j,]) +b1[i,1]
y2[r*(i-1)+j,2]=smooth_2(t[j,]) +b2[i,1]
y3[r*(i-1)+j,2]=smooth_3(t[j,]) +b3[i,1]
y4[r*(i-1)+j,2]=smooth_4(t[j,]) +b4[i,1]

}

}


dat_1=data.frame(y1)
dat_2=data.frame(y2)
dat_3=data.frame(y3)
dat_4=data.frame(y4)


colnames(dat_1)=cbind("id", "yij_logit","obstime")
colnames(dat_2)=cbind("id", "yij_logit","obstime")
colnames(dat_3)=cbind("id", "yij_logit","obstime")
colnames(dat_4)=cbind("id", "yij_logit","obstime")


dat_1$pij<-exp(dat_1$yij_logit)/(1+exp(dat_1$yij_logit))
dat_2$pij<-exp(dat_2$yij_logit)/(1+exp(dat_2$yij_logit))
dat_3$pij<-exp(dat_3$yij_logit)/(1+exp(dat_3$yij_logit))
dat_4$pij<-exp(dat_4$yij_logit)/(1+exp(dat_4$yij_logit))


dat_1$group=rep(1,nrow(dat_1))
dat_2$group=rep(2,nrow(dat_2))
dat_3$group=rep(3,nrow(dat_3))
dat_4$group=rep(4,nrow(dat_4))

p1=c(dat_1$pij)
p2=c(dat_2$pij)
p3=c(dat_3$pij)
p4=c(dat_4$pij)

dat_1$yij=rbinom(nrow(dat_1), 1, p1)
dat_2$yij=rbinom(nrow(dat_2), 1, p2)
dat_3$yij=rbinom(nrow(dat_3), 1, p3)
dat_4$yij=rbinom(nrow(dat_4), 1, p4)


###output datasets

data_binomial_all =data.frame(bind_rows(dat_1,dat_2,dat_3,dat_4))

data_binomial <- subset(data_binomial_all, select = c(id, obstime, yij, group))

##write.csv(data_binomial, paste0("binomial_medium", sn, ".csv"),  na="", col.names =T, row.names=F)

}


######################################################
################poisson distribution##################
######################################################


for (sn in 1:1000)  ##simulate 1000 datasets

{

set.seed(sn)


num=100

r=20

b_r=0.1  ##random effect for each participant which indicates high, medium, or low separation of trajectories

t <- as.matrix(rep(1:r))/r

smooth_1 <- function(x) 0.35*x^11*(10*(1 - x))^6+7*(10*x)^3*(1 - x)^10-5
smooth_2 <- function(x) 0.1*x^11*(10*(1 - x))^6+1*(10*x)^3*(1 - x)^10+1
smooth_3 <- function(x) 0*x^11*(10*(1 - x))^6+0*(10*x)^(-1)*(1 - x)^0+2.6
smooth_4 <- function(x) 1*x^1*(10*(1 - x))^0+0*(10*x)^3*(1 - x)^10+2.5

b1=as.matrix(rnorm(num,0,b_r))
b2=as.matrix(rnorm(num,0,b_r))
b3=as.matrix(rnorm(num,0,b_r))
b4=as.matrix(rnorm(num,0,b_r))

y1=matrix(,nrow=num*r, ncol=3)
y2=matrix(,nrow=num*r, ncol=3)
y3=matrix(,nrow=num*r, ncol=3)
y4=matrix(,nrow=num*r, ncol=3)


for (i in 1:num)
{
for (j in 1:r)
{

y1[r*(i-1)+j,1]=i
y2[r*(i-1)+j,1]=i+num
y3[r*(i-1)+j,1]=i+num*2
y4[r*(i-1)+j,1]=i+num*3

y1[r*(i-1)+j,3]=j/r
y2[r*(i-1)+j,3]=j/r
y3[r*(i-1)+j,3]=j/r
y4[r*(i-1)+j,3]=j/r

y1[r*(i-1)+j,2]=smooth_1(t[j,]) +b1[i,1]
y2[r*(i-1)+j,2]=smooth_2(t[j,]) +b2[i,1]
y3[r*(i-1)+j,2]=smooth_3(t[j,]) +b3[i,1]
y4[r*(i-1)+j,2]=smooth_4(t[j,]) +b4[i,1]

}

}



dat_1=data.frame(y1)
dat_2=data.frame(y2)
dat_3=data.frame(y3)
dat_4=data.frame(y4)


colnames(dat_1)=cbind("id", "yij_log","obstime")
colnames(dat_2)=cbind("id", "yij_log","obstime")
colnames(dat_3)=cbind("id", "yij_log","obstime")
colnames(dat_4)=cbind("id", "yij_log","obstime")


dat_1$pij<-exp(dat_1$yij_log)
dat_2$pij<-exp(dat_2$yij_log)
dat_3$pij<-exp(dat_3$yij_log)
dat_4$pij<-exp(dat_4$yij_log)



dat_1$group=rep(1,nrow(dat_1))
dat_2$group=rep(2,nrow(dat_2))
dat_3$group=rep(3,nrow(dat_3))
dat_4$group=rep(4,nrow(dat_4))

p1=c(dat_1$pij)
p2=c(dat_2$pij)
p3=c(dat_3$pij)
p4=c(dat_4$pij)

dat_1$yij=rpois(nrow(dat_1), p1)
dat_2$yij=rpois(nrow(dat_2), p2)
dat_3$yij=rpois(nrow(dat_3), p3)
dat_4$yij=rpois(nrow(dat_4), p4)


###output datasets

data_poisson_all =data.frame(bind_rows(dat_1,dat_2,dat_3,dat_4))

data_poisson <- subset(data_poisson_all, select = c(id, obstime, yij, group))

##write.csv(data_poisson, paste0("poisson_medium", sn, ".csv"),  na="", col.names =T, row.names=F)

}

