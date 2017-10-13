##Data Exploration Script
#09072017

#Load Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tables)
library(tidyr)
library(lme4)
library(DT)
library(lmerTest)
library(lcmm)
library(lavaan)
library(bnlearn)
library(flexmix)
library(glmmADMB)

#Load Data 
final_df2 <- read.table("~/Documents/Analysis/R Scripts/Symp_Data/final_df2.csv",header=TRUE, sep="|")
#final_df<-read.table("final_df2.csv",header=TRUE,sep="|")
table(final_df2$gender)
#final_df<-final_df%>%filter(gender=="MALE"|gender=="FEMALE")

final_long<-final_df2%>%
  select(id, gender, pre_headache_30:pos_vision_365)%>%
  gather(survey,score,-gender,-id)%>%
  separate(survey,c("session","symptom","interval"),"_")%>%
  mutate(interval=ifelse(session=="pre",as.numeric(interval)*-1,as.numeric(interval)))

final_long$interval<-as.character(final_long$interval)
final_long$gender<-factor(final_long$gender)
final_long$score<-as.numeric(final_long$score)
final_long$id<-factor(final_long$id)

means<-function(x){
  mean(x,na.rm=TRUE)
}
sds<-function(x){
  sd(x,na.rm=TRUE)
}
mean_long<-final_long%>%
  group_by(session,symptom,interval)%>%
  summarise(mean=means(score))
mean_long$interval<-as.numeric(mean_long$interval)
#ggplot(subset(mean_long,session="pos"),aes(x=interval,y=mean,colour=factor(symptom)))+geom_point()+stat_smooth(aes(group=symptom),formula=y~s(x,k=3),method="gam",se=FALSE)
ggplot(mean_long,aes(x=interval,y=mean,colour=symptom))+geom_point()+stat_smooth(aes(group=symptom),formula=y~poly(x,2),method="lm",se=FALSE)+
  ggtitle("Average Projections of Symptoms at Given Intervals")+theme(plot.title = element_text(hjust=0.5,vjust=1))+scale_x_continuous(breaks = c(-365,-180,-90,-30,30,90,180,365))

ggplot(subset(mean_long,session=="pos"),aes(y=mean,x=interval,colour=symptom))+geom_point()+stat_smooth(aes(group=symptom),formula=y~poly(x,2),method="lm",se=FALSE)+
  ggtitle("Average Projections of Symptoms at Given Intervals Post TBI")+theme(plot.title = element_text(hjust=0.5,vjust=1))+scale_x_continuous(breaks = c(30,90,180,365))

#Examining Distributions
ggplot(final_long,aes(x=score))+geom_histogram(bins=50) #overall distributions of scores
ggplot(final_long,aes(x=score))+geom_histogram()+facet_grid(~symptom) #distribution grouped by variable (change facet to designated variable)
ggplot(subset(final_long,session=="pos"),aes(x=score))+geom_histogram() #only post tbi screen
final_long$session<-factor(final_long$session,levels=c("pre","pos"))
ggplot(final_long,aes(x=score,fill=session))+geom_histogram(alpha=0.33,position="identity") #only post tbi screen

mean(final_long$score,na.rm=TRUE)
var(final_long$score,na.rm=TRUE)

#Descriptives
stats<-final_long%>%
  group_by(session, symptom, interval)%>%
  summarise(mean=round(means(score),3),sd=round(sds(score),3))%>%
  mutate(stat=paste(mean," (",sd,")",sep=""))

stats_w<-dcast(stats,session+interval~symptom,value.var = "stat")
kable(stats_w[c(7,5,8,6,2,4,1,3),],row.names = FALSE)

comp<-final_df2[complete.cases(final_df2[,c(2,7:78)]),]
comp<-comp[,c(2,7:78)]
cov(comp)
corrplot(cor(comp))
View(cor(comp))

#Model Experimentation
##Plotting Models without covariates and random repeated models 
final_long$interval<-as.numeric(final_long$interval)

p1<-glmm.admb(score~interval,data=final_long)







#NBIN





dep<-final_long%>%
  filter(symptom=="depression"& session=="pos")%>%
  filter(!is.na(score))
table(dep$interval)
trim<-final_long%>%
  filter(!is.na(score))

#ggplot(trim)+geom_point(aes(x=interval,y=score,group=id))
#ggplot(dep)+geom_point(aes(x=interval,y=score,group=id))

mlin<-lcmm(score~interval+gender+interval*gender,random=~interval,subject="id",data=dep,ng=3,mixture=~interval)
postprob(mlin)
save(mlin,file="mlin")
summary(mlin)

#https://www.r-bloggers.com/latent-class-mixed-models-with-graphics/

#Lavaan Data Modeling For Exprimentation

ggplot(dep,aes(x=interval,y=score,colour=factor(id)))+geom_point()+geom_line()

post<-final_df2[,grep("\\pos",names(final_df2))]
post<-post[,grep("\\_depression",names(post))]

post<-post[complete.cases(post[,c(1:36)]),] #subset with complete datasets
post<-as.data.frame(lapply(post,as.numeric)) #set as dataframe with numeric force 

pre<-final_df2[,grep("\\pre",names(final_df2))]

bl<-read.csv(file="/Users/peterhoover/Documents/Black_list.csv",sep=",")
bl<-data.frame(from=c("pos_depression_365","pos_depression_180","pos_depression_90","pos_depression_365","pos_vision_365"),to=c("pos_depression_180","pos_depression_90","pos_depression_30","pos_depression_90","pos_vision_180"))

plot(gs(post,blacklist=bl))


data('NPreg')
set.seed(123)

model_p<-FLXMRglm(score~interval|id,family="poisson")
m1<-flexmix(score~interval+gender+interval*gender,data=dep,k=2,model=FLXRziglm(family="poisson"))
m2<-flexmix(.~x,k=2,model=model_p,data=dep)
summary(m3)
m3<-flexmix(score~interval|id,data=dep,k=2)


formula<-"score~interval | id"
md1<-FLXMRziglm(family="poisson")
class<-FLXPmultinom(~interval)
mp<-flexmix(as.formula(formula),data=dep,k=2,model=md1)##
test<-flexmix(yn~x|id2,data=NPreg,k=2)
plot(refit(m3))

model_n<-FLXMRglm(yn~.+I(x^2),family="poisson")
flexmix(.~x,data=NPreg,k=2,model=model_n)

#Pull from RBloggers for ZIPR

library(flexmix)

set.seed(2013)
class <- FLXPmultinom(~ AGE + ACADMOS + MINORDRG + LOGSPEND)
formula <- "MAJORDRG ~ AGE + ACADMOS + MINORDRG + LOGSPEND"
control <- list(verbose = 10, iter.max = 500, minprior = 0.1, tol = 0.01)

cat("\n### TWO-CLASS FINITE MIXTURE POSSION MODEL ###\n")
mdl1 <- FLXMRglm(family = "poisson")
fit1 <- flexmix(as.formula(formula), data = data[data$CARDHLDR == 1, ], k = 2, model = mdl1, concomitant = class, control = control)
refit1 <- refit(fit1, method = 'optim')
cat("\n=== MODEL THE RESPONSE ===\n")
summary(refit1, which = 'model')
cat("\n=== MODEL THE MIXTURE DISTRIBUTION ===\n")
summary(refit1, which = 'concomitant') 

cat("\n### ZERO-INFLATED POSSION MODEL ###\n")
mdl2 <- FLXMRziglm(family = "poisson")
fit <- flexmix(as.formula(formula), data = data[data$CARDHLDR == 1, ], k = 2 , model = mdl2, concomitant = class, control = control)
refit2 <- refit(fit, method = 'optim')
cat("\n=== MODEL THE RESPONSE ===\n")
summary(refit2, which = 'model')
cat("\n=== MODEL THE MIXTURE DISTRIBUTION ===\n")
summary(refit2, which = 'concomitant')




library(glmmADMB)

#Ensuring working ZIP Poisson (Apply Negative Binomial as well)
data("dmft", package = "flexmix")
Model <- FLXMRziglm(family = "poisson")
Fitted <- flexmix(End ~ log(Begin + 0.5) + Gender + Ethnic + Treatment, 
                  model = Model, k = 2 , data = dmft, 
                  control = list(minprior = 0.01))
summary(refit(Fitted))

