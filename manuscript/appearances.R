library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape)

### MY OWN PEARSONS TEST
pearsons <- function(obs) {
  exp = 1/mean(obs)
  chistat = sum((obs - exp)^2 / exp)
  dof = length(obs) - 1
  p.val = pchisq(q=chistat,df=dof)
  return(p.val)
}

rho_func <- function(row) {
  if (row[5] == '1'){
    return(rho_sim(row))
  }
  if (row[5] == '10') {
    return(rho_cm(row))
  }
  return(1)
}

rho_sim <- function(row) {
  row <- as.numeric(row)
  H <- row[2]
  tau <- row[8]
  U <- row[3]
  s <- row[7]
  return(2 * s * H + 2 * (U ^ 2) * (1 - s) * (tau - 1) / s)
}

rho_cm <- function(row) {
  row <- as.numeric(row)
  H <- row[2]
  s <- row[7]
  return(2 * s * H )
}

plot.on.tau <- function(df, title='') {
    q <- ggplot(df, aes(x=tau,y=value,group=variable,color=variable,linetype=variable))+ 
                  geom_line()         
}

## READ DATA
df <- read.csv(file='appearances_pop1e8.csv.gz', header=T)
#df <- subset(df, s==0.05)

## CALC ADAPTATION
dff <- ddply(df, .(G,H,U,beta,pi,pop_size,s,tau, fname), summarize,
            T = unique(T))
adapt.df <- ddply(dff, .(G,H,U,beta,pi,pop_size,s,tau), summarize,
             time.mean = mean(T),
             time.sd = sd(T),
             time.se = sd(T)/sqrt(length(T)))
qplot(x=tau, y=time.mean, data=adapt.df, ymin=time.mean-time.se, ymax=time.mean+time.se, geom=c("point","errorbar"), facets=pi~s)      

## CALC APPEARANCE
appearance.df <- ddply(df, .(G,H,U,beta,pi,pop_size,s,tau), summarize,
             time.mean = mean(dif),
             time.sd = sd(dif),
             time.se = sd(dif)/sqrt(length(dif)),
             prob.mean = 1/mean(dif))
appearance.df <- ddply(appearance.df, .(G,H,U,beta,pi,pop_size,s,tau), transform,
              time.exact = 1/ (pop_size * ((U*beta)^2 * exp(-U/s - U)+2*tau*(U*beta)^2 / s *exp(-U*beta / s -tau*U))),
              time.approx = 1/(pop_size* ((1-tau*U)*2*tau*(U*beta)^2 /s)),
              time.mean.plus.se = time.mean + time.se,
              time.mean.minus.se = time.mean - time.se)
appearance.df <- melt(data=appearance.df,measure.vars=c("time.mean", "time.exact", "time.approx"))

q <- ggplot(appearance.df, aes(x=tau,y=value,group=variable,color=variable,linetype=variable))+ 
  geom_line() +
  geom_point() + 
  geom_errorbar()
q
## CALC FIXATION
dff <- ddply(df, .(G,H,U,beta,pi,pop_size,s,tau,fname), summarize,
             N = length(app),
             p.val = pearsons(N))
fixation.df <- ddply(dff, .(G,H,U,beta,pi,pop_size,s,tau), summarize,
             time.mean = mean(N),
             time.sd = sd(N),
             N = length(N),
             prob.mean = 1/mean(N))
fixation.df$time.se = fixation.df$time.sd/sqrt(fixation.df$N)

fixation.df <- cbind(fixation.df, rho=apply(fixation.df,1,rho_func))

fixation.df = subset(fixation.df, s==0.05)

q <- ggplot(data=subset(fixation.df,pi==1), mapping=aes(x=tau)) +
  geom_line(mapping=aes(y=1/rho), size=1) + 
  geom_point(mapping=aes(y=time.mean), size=2) + 
  geom_line(mapping=aes(y=time.mean), size=0.5) + 
  geom_errorbar(mapping=aes(ymax = time.mean + time.se, ymin=time.mean - time.se)) 

q + facet_grid(facets=s~pi, scales='free')



# http://stats.stackexchange.com/questions/1174/how-can-i-test-if-given-samples-are-taken-from-a-poisson-distribution

p.est <- 1/mean(test.df$N) ## estimate of probability parameter 
(tab.os<-table(test.df$N)) ## table with empirical frequencies
freq.os<-vector()
for(i in 1: length(tab.os)) freq.os[i]<-tab.os[[i]]  ## vector of emprical frequencies
freq.ex<-(dgeom(sort(unique(test.df$N)),prob=p.est)*length(test.df$N)) ## vector of fitted (expected) frequencies
acc <- mean(abs(freq.os-trunc(freq.ex))) ## absolute goodness of fit index acc
acc/mean(freq.os)*100 ## relative (percent) goodness of fit index

h <- hist(test.df$N ,breaks=length(tab.os))
xhist <- c(min(h$breaks),h$breaks)
yhist <- c(0,h$density,0)
xfit <- min(test.df$N):max(test.df$N)
yfit <- dgeom(xfit,prob=p.est)
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), main="Geom density and histogram")
lines(xfit,yfit, col="red")

#Perform the chi-square goodness of fit test 
#In case of count data we can use goodfit() included in vcd package
library(vcd) ## loading vcd package
gf <- goodfit(test.df$N,type= "geometric",method= "MinChisq")
summary(gf)
plot(gf,main="Count data vs Geom distribution")


