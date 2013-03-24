library(ggplot2)
library(plyr)
library(stringr)

plot.mean.fitness <- function(filename) {
  data <- read.csv(str_c("output/",filename))
  df <- ddply(data, .(fitness,tick), summarize, count=sum(population))
  ddf <- ddply(df, .(tick), summarize, mean.fitness = weighted.mean(fitness, count))
  return(qplot(x=tick,y=mean.fitness,data=ddf,geom=c("point")))
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
p <- plot.mean.fitness(filename)
png.name <- str_c("output/",str_replace(filename, ".csv.gz", ".png"))
ggsave(filename=png.name, plot=p)
print(png.name)