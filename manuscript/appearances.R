library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape)

ggplot2::theme_set(theme_bw())

plot.on.tau <- function(df, title='') {
    q <- ggplot(df, aes(x=tau,y=value,group=variable,color=variable,linetype=variable))+ 
                  geom_line() 
                  
}

## READ DATA
df <- read.csv(file='../stochastic/appearances.csv.gz', header=T)
df <- subset(df, pi==1 & s==0.05)

## CALC ADAPTATION
dff <- ddply(df, .(G,H,U,beta,pi,pop_size,s,tau, fname), summarize,
            T = unique(T))
adapt.df <- ddply(dff, .(G,H,U,beta,pi,pop_size,s,tau), summarize,
             time.mean = mean(T),
             time.sd = sd(T),
             time.se = sd(T)/sqrt(length(T)))
             


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

## PLOTS
adaptation.p <- plot.on.tau(adapt.df, 'Adaptation')
appearance.p <- plot.on.tau(appearance.df, 'Appearance')
fixation.p <- plot.on.tau(fixation.df, 'Fixation')




