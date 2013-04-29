library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape)

ggplot2::theme_set(theme_bw())

plot.on.tau <- function(df, title='') {
  q <- ggplot(data=df, mapping=aes(x=tau, y=time.mean, color=s, group=s)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymax = time.mean + time.se, ymin=time.mean - time.se)) + 
    scale_y_log10() +
    labs(x=expression(tau), y="Average waiting time", title=title) 
  return(q)
}

## READ DATA
df <- read.csv(file='appearances.csv.gz', header=T)
df <- subset(df, pi==1)

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




