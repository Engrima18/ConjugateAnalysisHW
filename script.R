load("Hmwk.RData")

# A.1 ---------------------------------------------------------------------

ex1 <- function(J){
  print(paste(c("The joint distriution sums to", sum(J))))
}


# A.2 ---------------------------------------------------------------------

derive.distr <- function(J){
  p.Y <- apply(J, 1, sum)
  p.Z <- apply(J, 2, sum)
  p.Y.given.Z <- t(t(J)/p.Z)
  p.Z.given.Y <- J/p.Y
  return(list(p.Y.given.Z = p.Y.given.Z, p.Z.given.Y = p.Z.given.Y,
              p.Y = p.Y, p.Z = p.Z))
}

show_distr <- function(J){
  xx <- derive.distr(J)
  p.Y <-  xx$p.Y
  p.Z <-  xx$p.Z
  p.Y.given.Z <- xx$p.Y.given.Z
  p.Z.given.Y <- xx$p.Z.given.Y
  # print marginals
  print(data.frame(p.Y, p.Z))
  
  # print conditionals
  print("----------------------------")
  print("p.Z.given.Y")
  print(p.Z.given.Y)
  
  print("----------------------------")
  print("p.Y.given.Z")
  print(p.Y.given.Z)
}

# A.3 ---------------------------------------------------------------------

verify_distr <- function(J){
  xx <- derive.distr(J)
  p.Y.given.Z <- xx$p.Y.given.Z
  p.Z.given.Y <- xx$p.Z.given.Y
  
  print("verify that each row of the p.Z.given.Y matrix sums to 1:")
  print(apply(p.Z.given.Y, 1, sum))
  
  print("verify that each column of the p.Y.given.Z matrix sums to 1:")
  print(apply(p.Y.given.Z, 2, sum))
}


# A.4 ---------------------------------------------------------------------

# first simulation method
sim1 <- function(J, n){
  probs <- c(J) # flattening the probability matrix
  # well define the support from which we sample
  support <- as.data.frame(t(expand.grid(y=1:3, z=1:3))) 
  samples.list <- as.data.frame(t(sample(support, n, replace=T, prob=probs))) # sample
  rownames(samples.list) <- NULL # reset the row indexes
  return(samples.list)
}

# second simulation method
sim2 <- function(J, n){
  distros <- derive.distr(J) # derive marginals and conditionals
  p.Y <- distros$p.Y # distr. of Y
  p.Z.given.Y <- distros$p.Z.given.Y # distr. of Z|Y
  y.sample <- sample(1:3, n, replace=T, prob=p.Y) # sample y from the marginal of Y
  z.sample <- rep(NA,n) # init the z sample
  for(idx in 1:n){
    y <- y.sample[idx]
    p.Z.given.y <- p.Z.given.Y[,y] # choose the correct conditional distr. 
    # sample from the conditional distr. given the sampled Y=y
    z.sample[idx] <- sample(1:3, 1, prob=p.Z.given.y) 
  }
  sample.list <- data.frame(y=y.sample, z=z.sample)
  return(sample.list)
}


# B.5 ---------------------------------------------------------------------


a = 3.007 # params of the prior
b = 1002.372
# observations
obs = c(1, 13, 27, 43, 73, 75 , 154, 196, 220, 297, 344, 610, 734, 783, 796, 845, 859, 992, 1066, 1471)
a.new = a+length(obs) # updated parameters
b.new = b+sum(obs)

E.psi.given.y <- function(a, b, obs){
  b.star <- b+sum(obs)
  a.star <- a+length(obs)
  return( b.star / (a.star-1) )
}

Var.psi.given.y <- function(a, b, obs){
  b.star <- b+sum(obs)
  a.star <- a+length(obs)
  return( (b.star^2) / ((a.star-1)^2 * (a.star-2)) )
}

bayesian.mean.lifetime <- E.psi.given.y(a, b, obs) # E[psi]
sample.mean <- sum(obs)/length(obs) # sample mean
variance.mean.life <- Var.psi.given.y(a, b, obs) # Var(psi)