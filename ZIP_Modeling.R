##09.25.2017
#P.Hoover 

#Libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(DT)
library(flexmix)
library(plotly)
library(knitr)
library(corrplot)
#library(glmmADMB)

#Load Data 
#final_df2 <- read.table("~/Documents/Analysis/R Scripts/Symp_Data/final_df2.csv",header=TRUE, sep="|")
#save(final_df2,file="final_df.RData")
#final_df<-read.table("final_df2.csv",header=TRUE,sep="|")
load(file="final_df.RData")
table(final_df2$gender)
#final_df<-final_df%>%filter(gender=="MALE"|gender=="FEMALE")

#Create long version of data for computations
final_long<-final_df2%>%
  select(id, gender, age_group,pre_headache_30:pos_vision_365)%>% #select columns of interest
  gather(survey,score,-gender,-id,-age_group)%>% #set to long version
  separate(survey,c("session","symptom","interval"),"_")%>% #separate variables based on session and symptom
  mutate(interval=ifelse(session=="pre",as.numeric(interval)*-1,as.numeric(interval)))

#Set variable types
final_long$interval<-as.numeric(final_long$interval)
final_long$gender<-factor(final_long$gender)
final_long$score<-as.numeric(final_long$score)
final_long$id<-factor(final_long$id)
str(final_long) #verify variable types

means<-function(x){ #function for dplyr
  mean(x,na.rm=TRUE)}
sds<-function(x){ #function for dplyr
  sd(x,na.rm=TRUE)
}
mean_long<-final_long%>%
  group_by(session,symptom,interval)%>%
  summarise(mean=means(score))
mean_long$interval<-as.numeric(mean_long$interval)

#ggplot(subset(mean_long,session="pos"),aes(x=interval,y=mean,colour=factor(symptom)))+geom_point()+stat_smooth(aes(group=symptom),formula=y~s(x,k=3),method="gam",se=FALSE)
ggplot(mean_long,aes(x=interval,y=mean,colour=symptom))+geom_point()+stat_smooth(aes(group=symptom),formula=y~poly(x,2),method="lm",se=FALSE)+
  ggtitle("Average Projections of Symptoms at Given Intervals")+theme(plot.title = element_text(hjust=0.5,vjust=1))+scale_x_continuous(breaks = c(-365,-180,-90,-30,30,90,180,365))

g<-ggplot(subset(mean_long,session=="pos"),aes(y=mean,x=interval,colour=symptom))+geom_point()+stat_smooth(aes(group=symptom),formula=y~poly(x,2),method="lm",se=FALSE)+
  ggtitle("Average Projections of Symptoms at Given Intervals Post TBI")+theme(plot.title = element_text(hjust=0.5,vjust=1))+scale_x_continuous(breaks = c(30,90,180,365))

ggplotly(g)

#Solely Explore Post-TBI
final_pos<-final_long%>%
  filter(session=="pos")%>%
  filter(!is.na(score))

graph<-final_pos%>%
  #filter(symptom=="depression")%>% 
  filter(id %in%sample(unique(id),1000))

ggplot(graph,aes(y=score,x=interval))+geom_point()+geom_line(aes(group=id))+
  ggtitle("Average Projections of Symptoms at Given Intervals Post TBI")+theme(plot.title = element_text(hjust=0.5,vjust=1))+scale_x_continuous(breaks = c(30,90,180,365))+facet_grid(~symptom)

#Examining Distributions of "Score"
ggplot(final_long,aes(x=score))+geom_histogram(bins=50)+scale_y_continuous(expand=c(0,0)) #overall distributions of scores
ggplot(final_long,aes(x=score))+geom_histogram()+facet_grid(~symptom)+scale_y_continuous(expand=c(0,0)) #distribution grouped by variable (change facet to designated variable)
ggplot(subset(final_long,session=="pos"),aes(x=score))+geom_histogram()+scale_y_continuous(expand=c(0,0))  #only post tbi screen
final_long$session<-factor(final_long$session,levels=c("pre","pos"))
ggplot(final_long,aes(x=score,fill=session))+geom_histogram(alpha=0.33,position="identity")+scale_y_continuous(expand=c(0,0)) #only post tbi screen


mean(final_long$score,na.rm=TRUE) 
var(final_long$score,na.rm=TRUE)

#Verifying Plot functions
ggplot(graph,aes(y=score,x=interval,colour=id,shape=symptom))+geom_line()+theme(legend.position="none")+geom_point()+scale_shape_manual(values=1:9)

#Descriptives
stats<-final_long%>%
  group_by(session, symptom, interval)%>%
  summarise(mean=round(means(score),3),sd=round(sds(score),3))%>%
  mutate(stat=paste(mean," (",sd,")",sep=""))

#Computing table for viewing means and sd 
stats_w<-dcast(stats,session+interval~symptom,value.var = "stat")
kable(stats_w[c(7,5,8,6,2,4,1,3),],row.names = FALSE)

#Simple coviarnace and correlation matrix to explore variables (we'd expect some correlation with a repeated study)
comp<-final_df2[complete.cases(final_df2[,c(2,7:78)]),] #Need to remove incomplete cases so we can explore the covariances
comp<-comp[,c(2,7:78)] #Subset only numerical components 
cov(comp) #Covariance Matrix
corrplot(cor(comp),tl.cex = .25) #covairance plot 
View(cor(comp)) #View correlation matrix 

#Model Experimentation
FLXMRnegbin <- function(formula = . ~ ., theta = NULL, offset = NULL, control = list(reltol = .Machine$double.eps^(1/1.5), maxit = 500)){ 
  .theta <- theta
  nbrefit <- function(x, y, w) {
    fit <- c(glm.fit(x, y, weights = w, offset = offset, family = MASS::negative.binomial(theta)),
             list(call = sys.call(), offset = offset,
                  control = eval(formals(glm.fit)$control),
                  method = "weighted.glm.fit"))
    fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank - is.null(.theta)
    fit$df.residual <- sum(w) - fit$rank - is.null(.theta)
    fit$x <- x
    fit
  }
  z <- methods::new("FLXMRglm", weighted = TRUE, formula = formula,
                    name = "FLXMRglm: negative.binomial", offset = offset,
                    family = "negative.binomial", refit = nbrefit)
  z@preproc.y <- function(x){
    if (ncol(x) > 1L)
      stop(paste("for the", family, "family y must be univariate"))
    x
  }
  z@defineComponent <- if(is.null(.theta)) {
    expression({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      logLik <- function(x, y, ...)
        suppressWarnings(dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE))
      methods::new("FLXcomponent",
                   parameters = list(coef = coef, theta = theta),
                   logLik = logLik,
                   predict = predict,
                   df = df)
    })
  } else {
    as.expression(substitute({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      logLik <- function(x, y, ...)
        suppressWarnings(dnbinom(y, mu = predict(x, ...), size = theta, log = TRUE))
      methods::new("FLXcomponent",
                   parameters = list(coef = coef),
                   logLik = logLik,
                   predict = predict,
                   df = df)
    }, as.environment(list(theta = .theta))))
  }
  z@fit <- function(x, y, w, component){
    if(is.null(component$theta)) {
      df <- ncol(x)
      theta <- if(is.null(.theta)) 1 else .theta
      cf <- glm.fit(x, y, weights = w, family = MASS::negative.binomial(theta),
                    offset = offset, start = component$coef)$coefficients
    } else {
      ## degrees of freedom
      df <- ncol(x) + 1
      ## offset
      if(is.null(offset)) offset <- 0
      ## objective function
      nll <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
        mu <- exp(drop(x %*% beta + offset))
        suppressWarnings(-sum(w * dnbinom(y, mu = mu, size = theta, log = TRUE)))
      }
      ## corresponding gradients
      gr <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
        mu <- exp(drop(x %*% beta + offset))
        gr <- drop(y - mu * (y + theta)/(mu + theta))
        colSums(-w * cbind(gr * x, theta * (digamma(y + theta) - digamma(theta) +
                                              log(theta) + 1 - log(mu + theta) - (y + theta)/(mu + theta))))
      }
      ## starting values from previous iteration
      start <- c(component$coef, component$theta)
      ## if not available: use geometric GLM
      if(length(start) < df) start <- c(
        glm.fit(x, y, weights = w, family = MASS::negative.binomial(1), offset = offset)$coefficients,0
      )
      ## BFGS optimization
      opt <- optim(par = start, fn = nll, gr = gr, method = "BFGS", control = control)
      ## estimated parameters
      cf <- opt$par[-df]
      theta <- exp(opt$par[df])
    }
    with(list(coef = cf, theta = theta, df = ncol(x) + is.null(.theta)),
         eval(z@defineComponent))
  }
  z
} #NB Formula from git (https://github.com/rforge/countreg/blob/master/pkg/R/FLXMRnegbin.R)

#Set Model specification for later use 
zip<-FLXMRziglm(family = "poisson")
nb=FLXMRnegbin()

#Set subset for calculations
set.seed(111)
final_pos<-final_pos%>%
  #filter(symptom=="depression")%>% 
  filter(id %in%sample(unique(id),1000))

#Model Evaluations
LR_test()

control <- list(verbose = NULL, iter.max = 500, minprior = 0.1, tol = 0.01) #Set control variables

#1a  nconditional means (UM) model to verify distribution
#• Intercept-only model (no time)
#• Baseline for unconditional growth model

#install.packages("countreg", repos="http://R-Forge.R-project.org")
library(pscl)
library(countreg)

m1<-flexmix(score~1,data=final_pos,model=zip,k=1)
m2<-flexmix(score~1, data=final_pos,model=nb,k=1)
m3<-hurdle(score~1, data=final_pos,dist="poisson")
m4<-hurdle(score~1,data=final_pos,dist="negbin")

#m5<-flexmix(score~1,data=final_pos,model=hurdle(score~1),k=1)

m6<-zeroinfl(score~1,data=final_pos,dist="negbin")
m7<-hurdle(score~1,data=final_pos,dist="poisson",zero="poisson")

AIC(m1,m2,m3,m4,m6,m7) #Neg bin appears to support the data
BIC(m1,m2,m3,m4,m6,m7)


#1b Unconditional means to determine groups
unc_zip<-flexmix(score~1,data=final_pos,model=zip,k=1)
unc_nb<-flexmix(score~1,data=final_pos,model=nb,k=1)
unc_pois<-flexmix(score~1,data=final_pos,model=FLXMRglm(family="poisson"),k=1)
unc_hr<-hurdle(score~1,data=final_pos,dist="poisson")
unc_hr_nb<-hurdle(score~1,data=final_pos,dist="negbin")
unc_hr_z<-hurdle(score~1,data=final_pos,dist="poisson",zero.dist = "poisson")


BIC(unc_zip,unc_nb,unc_pois,unc_hr,unc_hr_nb)
AIC(unc_zip,unc_nb,unc_hr,unc_hr_nb)

#Negative binomial providing best model 

unc_zip<-stepFlexmix(score~1,data=final_pos,model=zip,k=1:6)
unc_nb<-stepFlexmix(score~1,data=final_pos,model=nb,k=1:6)

BIC(unc_zip,unc_nb)
AIC(unc_zip,unc_nb)


#2Unconditional Growth Models 
#• Growth curve model with effect of time* (no additional covariates)
#• Slopes can be set as randomly varying* or fixed
#• Baseline for conditional growth model

zip_grow<-flexmix(score~interval,data=final_pos,model=zip,k=1)
nb_grow<-flexmix(score~interval,data=final_pos,model=nb,k=1)
pos_grow<-flexmix(score~interval,data=final_pos,model=FLXMRglm(family="poisson"),k=1)
BIC(zip_grow,nb_grow,pos_grow)

#Non-linear model
zip_grow2<-flexmix(score~interval+I(interval^2),data=final_pos,model=zip,k=1)
nb_grow2<-flexmix(score~interval+I(interval^2),data=final_pos,model=nb,k=1)
pos_grow2<-flexmix(score~interval+I(interval^2),data=final_pos,model=FLXMRglm(family="poisson"),k=1)
BIC(zip_grow2,nb_grow2,pos_grow2)

zip_grow3<-flexmix(score~interval+I(interval^3),data=final_pos,model=zip,k=1)
nb_grow3<-flexmix(score~interval+I(interval^3),data=final_pos,model=nb,k=1)
pos_grow3<-flexmix(score~interval+I(interval^3),data=final_pos,model=FLXMRglm(family="poisson"),k=1)
BIC(zip_grow3,nb_grow3,pos_grow3)


#Latent Growth
unc_zip_t<-stepFlexmix(score~interval,data=final_pos,model=zip,k=1:6)
unc_nb_t<-stepFlexmix(score~interval,data=final_pos,model=nb,k=1:6)
unc_pos_t<-stepFlexmix(score~interval,data=final_pos,model=FLXMRglm(family="poisson"),k=1:6)

unc_zip_t<-getModel(unc_zip_t)
unc_nb_t<-getModel(unc_nb_t)
unc_pos_t<-getModel(unc_pos_t)

BIC(unc_zip_t,unc_nb_t,unb_pos_t)


unc_zip_t2<-stepFlexmix(score~interval+I(interval^2),data=final_pos,model=zip,k=1:6)
unc_nb_t2<-stepFlexmix(score~interval+I(interval^2),data=final_pos,model=nb,k=1:6, control=control)
unc_pos_t2<-stepFlexmix(score~interval+I(interval^2),data=final_pos,model=FLXMRglm(family="poisson"),k=1:6)

unc_zip_t2<-getModel(unc_zip_t2)
unc_nb_t2<-getModel(unc_nb_t2)
unc_pos_t2<-getModel(unc_pos_t2)

BIC(unc_zip_t2,unc_nb_t2,unb_pos_t2)



#3Conditional growth  model(s)
ug_zip<-stepFlexmix(score~interval|id,data=final_pos,model=zip,k=1:6)
ug_nb<-stepFlexmix(score~interval|id,data=final_pos,model=nb,k=1:6, control=control)

BIC(getModel(ug_zip))

BIC(nb1)

summary(getModel(ug_zip))
summary(getModel(ug_nb))

BIC(getModel(ug_zip))
BIC(getModel(ug_nb))

#4Conditional growth (CG) model(s)
#• Add time-invariant and/or time-varying
cg_zip<-stepFlexmix(score~interval+symptom|id,data=final_pos,model=zip,k=1:6)

#test_data<-cbind(final_pos,getModel(cg_zip)@cluster)
#write.csv(test_data,file="test_data.csv")


cg_nb<-stepFlexmix(score~interval+symptom|id,data=final_pos,model=nb,k=1:6)
#cg_nb<-stepFlexmix(score~interval+symptom|id,data=final_pos,model=nb,k=1:6,control=control)
cg_nb
summary(getModel(cg_nb))

cg_p<-stepFlexmix(score~interval+symptom|id,data=final_pos,model=FLXMRglm(family="poisson"),k=1:6)
summary(getModel(cg_p))

#Model Extensions
cg_zip_g<-stepFlexmix(score~interval+symptom+gender|id,data=final_pos,model=zip,k=1:6)
BIC(getModel(cg_zip_g))
summary(getModel(cg_zip_g))

cg_nb_g<-stepFlexmix(score~interval+symptom+gender|id,data=final_pos,model=nb,k=1:6)

test<-stepFlexmix(score~interval+symptom+gender|id,data=final_pos,model = FLXMRglm(family = "poisson"),k=1:6)
BIC(getModel(test))
summary(getModel(test))

test2<-stepFlexmix(score~interval+symptom|id,data=final_pos,model = FLXMRglmfix(family = "poisson",fixed=~gender),k=1:6)
BIC(getModel(test2))

test3<-stepFlexmix(score~interval+symptom+gender+age_group|id,data=final_pos,model = FLXMRglm(family = "poisson"),k=1:6)
BIC(getModel(test3))


model<-getModel(test3)


ggplot(graph_temp,(aes(x=factor(interval),y=score,colour=factor(model@cluster))))+geom_point()+geom_line()
head(final_pos)





