###################################################
# Functions for R species interaction analysis
# by Lukas Ullmann
# lukas.ullmann@campus.lmu.de
#
# R version: 4.3.0
# Project: Bachelor's thesis Lukas Ullmann
# Thesis supervisor: Alexander Keller (LMU München)
###################################################


#load required packages
library(igraph)


#Create interaction networks
##nonrandom: only abundance based
networks.nonrandom <- function(timerange, partner1, partner2){
  net.list <- list()
  for (i in timerange){
    p1.matrix = matrix(nrow=length(colnames(partner1)),ncol=length(colnames(partner2)),rep.int(unlist(partner1[i,]),length(colnames(partner2))))
    p2.matrix = t(matrix(nrow=length(colnames(partner2)),ncol=length(colnames(partner1)),rep.int(unlist(partner2[i,]),length(colnames(partner1)))))
    
    rownames(p1.matrix) <- colnames(partner1)
    colnames(p2.matrix) <- colnames(partner2)
     
    nonrandom.interaction <- p1.matrix*p2.matrix
    
    colnames(nonrandom.interaction) <- colnames(partner2)
    
    nonrandom.interaction <- nonrandom.interaction/sum(nonrandom.interaction)
    
    nonrandom.interaction[is.nan(nonrandom.interaction)] <- 0
    
    net.list[[i]] <- nonrandom.interaction
  }
  return(net.list)
}


##random: abundance + random effect
networks.random <- function(timerange, partner1, partner2){
  net.list <- list()
  for (i in timerange){
    net <- as_incidence_matrix(sample_bipartite(length(colnames(partner1)), length(colnames(partner2)), p=0.2))
    
    net <- round(net * (rnorm(n=length(net), 30, 20)+30), digits=0)
    
    rownames(net) <- names(partner1)
    colnames(net) <- names(partner2)
    
    p1.matrix = matrix(nrow=length(colnames(partner1)),ncol=length(colnames(partner2)),rep.int(unlist(partner1[i,]),length(colnames(partner2))))
    p2.matrix = t(matrix(nrow=length(colnames(partner2)),ncol=length(colnames(partner1)),rep.int(unlist(partner2[i,]),length(colnames(partner1)))))
    
    rownames(p1.matrix) <- colnames(partner1)
    colnames(p2.matrix) <- colnames(partner2)
    
    nonrandom.interaction <- p1.matrix*p2.matrix
    random.interaction <- nonrandom.interaction * net
    
    colnames(random.interaction) <- colnames(partner2)
    
    random.interaction <- random.interaction/sum(random.interaction)
    
    random.interaction[is.nan(random.interaction)] <- 0
    
    net.list[[i]] <- random.interaction
  }
  return(net.list)
}


##fullrandom: only random effect
networks.fullrandom <- function(timerange, partner1, partner2){
  net.list <- list()
  for (i in timerange){
    net <- as_incidence_matrix(sample_bipartite(length(colnames(partner1)), length(colnames(partner2)), p=0.2))
    
    net <- round(net * (rnorm(n=length(net), 30, 20)+30), digits=0)
    
    rownames(net) <- names(partner1)
    colnames(net) <- names(partner2)
    
    p1.matrix = matrix(nrow=length(colnames(partner1)),ncol=length(colnames(partner2)),rep.int(unlist(partner1[i,]),length(colnames(partner2))))
    p2.matrix = t(matrix(nrow=length(colnames(partner2)),ncol=length(colnames(partner1)),rep.int(unlist(partner2[i,]),length(colnames(partner1)))))
    
    nonrandom.interaction <- p1.matrix*p2.matrix
    fullrandom.interaction <- sample(nonrandom.interaction) * net
    
    colnames(fullrandom.interaction) <- colnames(partner2)
    rownames(fullrandom.interaction) <- colnames(partner1)
    
    fullrandom.interaction <- fullrandom.interaction/sum(fullrandom.interaction)
    
    fullrandom.interaction[is.nan(fullrandom.interaction)] <- 0
    
    net.list[[i]] <- fullrandom.interaction
  }
  return(net.list)
}


##null: networks for null model
###random: abundance + random effect
networks.null.random <- function(timerange, n_random, partner1, partner2){
  net.list.null <- list()
  fullrandom.interaction <- list()
  
  for (j in 1:n_random){
    fullrandom.interaction <- networks.random(timerange, partner1, partner2)
    
    net.list.null[[j]] <- fullrandom.interaction 
  }
  return(net.list.null)
}

###fullrandom: only random effect
networks.null.fullrandom <- function(timerange, n_random, partner1, partner2){
  net.list.null <- list()
  fullrandom.interaction <- list()
  
  for (j in 1:n_random){
    fullrandom.interaction <- networks.fullrandom(timerange, partner1, partner2)
    
    net.list.null[[j]] <- fullrandom.interaction 
  }
  return(net.list.null)
}


#Calculate modes with time stamps
modes.time <- function(partner1){
  partner1.modes <- apply(partner1, 2, function(column) {
    mode_time <- as.vector(as.numeric(rownames(partner1)[which.max(column)]))
    mode <- as.vector(column[which.max(column)])
    c(Mode = mode, Time = mode_time)
  })
  return(partner1.modes)
}


#Calculate (by interaction strength) weighted interaction partner means
partner.means.weighted <- function(timerange, partner1, partner2, networks){
  partner1.partner.means.weighted <- data.frame(matrix(ncol = length(colnames(partner1)), nrow = length(rownames(partner1))))
  colnames(partner1.partner.means.weighted) <- colnames(partner1)
  
  partner2_means <- modes.time(partner2)
  
  for (i in timerange){
    mean_means_partner1 <- apply(networks[[i]], 1, function(row) {
      non_zero_partner2 <- names(row[row != 0])
      non_zero_means <- partner2_means["Time", non_zero_partner2]
      non_zero_interactions <- row[row != 0]
      sum(as.numeric(non_zero_means) * non_zero_interactions) / sum(row[row != 0])
    })
    partner1.partner.means.weighted[i, ] <- mean_means_partner1
  }
  
  return(partner1.partner.means.weighted)
}


#calculate difference between weighted means of partner1 and all partners
diff.weighted <- function(timerange, partner1, partner2, networks){
  diff.partner1.weighted <- data.frame(matrix(ncol = length(colnames(partner1)), nrow = length(rownames(partner1))))
  colnames(diff.partner1.weighted) <- colnames(partner1)
  
  for (i in timerange){
    diff <- partner.means.weighted(timerange, partner1, partner2, networks)[i,] - modes.time(partner1)["Time",]
    diff.partner1.weighted[i, ] <- diff
  }
  return(diff.partner1.weighted)
}


#calculate a (by relative abundance) weighted mean value for differences
diff.weighted.total <- function(timerange, partner1, partner2, networks){
  diff <- diff.weighted(timerange, partner1, partner2, networks)
  diff.total <- colSums(diff * partner1, na.rm = TRUE) / colSums(partner1, na.rm = TRUE)
  
  return(diff.total)
}


#create null distribution
diff.weighted.null <- function(timerange, n_random, partner1, partner2){
  #replace "networks.null.random" with "networks.null.fullrandom" to remove influence of abundance on distribution
  net.list.null <- networks.null.random(timerange, n_random, partner1, partner2)
  
  totals <- data.frame(matrix(ncol = length(colnames(partner1)), nrow = length(net.list.null)))
  colnames(totals) <- colnames(plants.rel)
  
  for (k in 1:nrow(totals)){
    total <- diff.weighted.total(timerange, partner1, partner2, net.list.null[[k]])
    totals[k, ] <- total
  }
  return(totals)
}


#summarize results
analysis.summary.direct <- function(partner1, partner2, networks, n_random = 100, a = 0.05){
  if((length(rownames(partner1)) == length(rownames(partner2)) & length(rownames(partner2)) == length(networks)) == FALSE){
    stop("WARNING: interaction partners and networks must be of equal length")
  }
  if(any((colnames(partner1) == rownames(networks[[1]])) == FALSE)){
    stop("WARNING: colnames of partner1 must be equal to rownames of networks")
  }
  if(any((colnames(partner2) == colnames(networks[[1]])) == FALSE)){
    stop("WARNING: colnames of partner2 must be equal to colnames of networks")
  }
  
  timerange <- 1:length(networks)
  total_means <- diff.weighted.total(timerange, partner1, partner2, networks)
  totals_random <- diff.weighted.null(timerange, n_random, partner1, partner2)
  interval_borders <- c(a/2, 1-a/2)
  
  interaction.results <- data.frame(matrix(ncol = length(colnames(total_means)), nrow = 0))
  
  for (i in 1:length(names(total_means))) {
    
    partner1_name <- names(total_means[i])
    partner1_value <- unlist(total_means[i])
    distribution <- totals_random[[partner1_name]]
    dist <-  partner1_value - mean(distribution)
    interval <- c(qnorm(a/2, mean(distribution), sd(distribution)), qnorm(1-a/2, mean(distribution), sd(distribution)))
    
    
    interaction.results[i, "Partner1"] <- partner1_name
    interaction.results[i, "Mean.partner.modes.null"] <- as.numeric(mean(totals_random[[i]]))
    interaction.results[i, "Interval.null.lower"] <- as.numeric(interval[1])
    interaction.results[i, "Interval.null.upper"] <- as.numeric(interval[2])
    interaction.results[i, "Mean.partner.modes.sample"] <- as.numeric(partner1_value)
    interaction.results[i, "Distance.between.means"] <- as.numeric(dist)
    interaction.results[i, "Probability"] <- dnorm(partner1_value, mean(distribution), sd(distribution))
    interaction.results[i, "Within.interval"] <- as.character(partner1_value > interval[1] & partner1_value < interval[2])
  }
  
 
  parameters <- data.frame(matrix(ncol = 1, nrow = 0))
  parameters["Time.steps", 1] <- format(length(timerange), scientific = FALSE)
  parameters["N.Partner1", 1] <- format(length(colnames(partner1)), scientific = FALSE)
  parameters["N.Partner2", 1] <- format(length(colnames(partner2)), scientific = FALSE)
  parameters["Significance.level", 1] <- format(a, scientific = FALSE)
  parameters["N.random.networks", 1] <- format(n_random, scientific = FALSE)
  colnames(parameters) <- "Value"
  
  
  output <- list(Means.Partner1 = total_means, Means.null.distribution = totals_random, Parameters = parameters, Interaction.results = interaction.results)
  
  return(output)
}


