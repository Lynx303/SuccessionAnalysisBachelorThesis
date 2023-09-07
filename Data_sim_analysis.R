###################################################
# R species interaction analysis
# by Lukas Ullmann
# lukas.ullmann@campus.lmu.de
#
# R version: 4.3.0
# Project: Bachelor's thesis Lukas Ullmann
# Thesis supervisor: Alexander Keller (LMU München)
###################################################


#load required packages
library(ggplot2)
library(stringr)


#clear environment
rm(list = ls())

setwd()
source("./Thesis_functions.R")
#-------------------------------------------------------------------------------

#### create fake data

### Abundances
## plants
plant_list = list()


for (plant in 1:10){
  plant_list[[plant]] <- table(round(rnorm(100, mean = 10+plant, sd = 3), digits=0))
}

plant_list[[31]] <-table(round(runif(10000,min=1,max=30), digits=0))

plants <- data.frame(t(dplyr::bind_rows(plant_list)))
plants[is.na(plants)] <- 0

plants <- plants[,1:10]

colnames(plants) <- paste("Plant", str_pad(1:length(colnames(plants)), 2, pad = "0"), sep="" )


plants <- plants[order(as.numeric(rownames(plants))), ]


## bees
bee_list = list()


for (bee in 1:10){
  bee_list[[bee]] <- table(round(rnorm(100, mean = 10+bee, sd = 3), digits=0))
}

bee_list[[31]] <-table(round(runif(10000,min=1,max=30), digits=0))

bees <- data.frame(t(dplyr::bind_rows(bee_list)))
bees[is.na(bees)]<-0

bees <- bees[,1:10]

colnames(bees) <- paste("Bee",str_pad(1:length(colnames(bees)), 2, pad = "0"), sep="" )

bees <- bees[order(as.numeric(rownames(bees))), ]


## turn into relative data
plants.rel <- t(apply(plants[1:30,], 1, function(value) 
  value/sum(value)
))

plants.rel[is.na(plants.rel)]  <- 0

bees.rel <- t(apply(bees[1:30,], 1, function(value) 
  value/sum(value)
))

bees.rel[is.na(bees.rel)] <- 0


### Networks
net.list.random <- networks.random(1:30, plants.rel, bees.rel)
net.list.nonrandom <- networks.nonrandom(1:30, plants.rel, bees.rel)
net.list.null <- networks.null.random(1:30, 100, plants.rel, bees.rel)
net.list.fullrandom <- networks.fullrandom(1:30, plants.rel, bees.rel)

#-------------------------------------------------------------------------------


###testing if within 95% of random distribution
results <- analysis.summary.direct(plants.rel, bees.rel, net.list.random, 100, 0.05)
results

###plotting
##distance between null distribution mean and partner interaction mean
data <- results$Interaction.results
sig.level <- (1 - as.numeric(results$Parameters[["Significance.level",1]])) * 100 

ggplot(data, aes(x = Distance.between.means, y = Partner1, color = Within.interval)) +
  geom_point(shape = 18, size = 5) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "blue")) +
  geom_segment(aes(xend = 0, yend = Partner1, color = Within.interval), linetype = "dashed", linewidth = 2) +
  geom_vline(xintercept = 0, color = "red", linewidth = 1.5) +
  labs(x = "Distance between means", color = paste0("Within ", sig.level, "% \nInterval")) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 40),
    legend.title = element_text(size=25),
    legend.text = element_text(size = 15),
    axis.text.y = element_text(size = 35),
    axis.text.x = element_text(size = 30)
  )



