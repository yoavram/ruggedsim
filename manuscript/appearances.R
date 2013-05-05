library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape)

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
df <- read.csv(file='../stochastic/appearances_AB0.csv.gz', header=T)
#df <- subset(df, s==0.01)

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

## CALC FIXATION
dff <- ddply(df, .(G,H,U,beta,pi,pop_size,s,tau,fname), summarize,
             N = length(app))
fixation.df <- ddply(dff, .(G,H,U,beta,pi,pop_size,s,tau), summarize,
             time.mean = mean(N),
             time.sd = sd(N),
             time.se = sd(N)/sqrt(length(N)),
             prob.mean = 1/mean(N))
fixation.df <- cbind(fixation.df, rho=apply(fixation.df,1,rho_func))

q <- ggplot(data=subset(fixation.df,pi==1), mapping=aes(x=tau)) +
  geom_line(mapping=aes(y=1/rho), size=1) + 
  geom_point(mapping=aes(y=time.mean), size=2) + 
  geom_line(mapping=aes(y=time.mean), size=0.5) + 
  geom_errorbar(mapping=aes(ymax = time.mean + time.se, ymin=time.mean - time.se)) 

q + scale_y_log10() + 
  geom_ribbon(mapping=aes(y=time.mean, ymin = time.mean - 2*time.se, ymax = time.mean + 2*time.se),alpha=I(0.2))

