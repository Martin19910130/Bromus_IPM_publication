##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                AIC
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()
options(stringsAsFactors = F)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###           read data
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Demography.csv")
seed_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Seed_data.csv")

## Delete the new plants which grow bigger then 30
demo_data <- demo_data[-(which(demo_data$sizet1 >= 30 & demo_data$new_plant == 1)),]

## AIC calculations
AIC_pub <- function(df, dep)
{
  
  if(dep == "logsizet1")
  {
    z1 <- lm(df[, dep] ~ df[, "logsizet0"]) 
    z2 <- lm(df[, dep] ~ df[, "logsizet0"] + df[, "climate"]) 
    z3 <- lm(df[, dep] ~ df[, "logsizet0"] + df[,"management"]) 
    z4 <- lm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] * df[,"management"]) 
    z5 <- lm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] + df[,"management"])
  }
  
  if(dep == "survival" | dep == "flower")
  {
    z1 <- glm(df[, dep] ~ df[, "logsizet0"], family = "binomial")
    z2 <- glm(df[, dep] ~ df[, "logsizet0"] + df[, "climate"], family = "binomial")
    z3 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"management"], family = "binomial") 
    z4 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] * df[,"management"], family = "binomial") 
    z5 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] + df[,"management"], family = "binomial")
  }
  
  if(dep == "number_of_seeds")
  {
    z1 <- glm(df[, dep] ~ df[, "logsizet0"], family = "poisson")
    z2 <- glm(df[, dep] ~ df[, "logsizet0"] + df[, "climate"], family = "poisson")
    z3 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"management"], family = "poisson") 
    z4 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] * df[,"management"], family = "poisson") 
    z5 <- glm(df[, dep] ~ df[, "logsizet0"] + df[,"climate"] + df[,"management"], family = "poisson")
  }
  
  x <- c(AIC(z1), AIC(z2), AIC(z3), AIC(z4),AIC(z5)) 
  delta <- x - min(x)                   
  L <- exp(-0.5 * delta)            
  w <- L/sum(L)                     
  
  w <- t(as.data.frame(w))
  colnames(w) <- c("Null", "climate", "management", "climate*management", "climate+management")
  return(w)  
}

dep <- c("survival", "logsizet1", "flower", "number_of_seeds")
mod_compare <- c()
for(i in 1:length(dep))
{
  mod_compare[i] <- list(AIC_pub(demo_data, dep[i]))
  names(mod_compare)[i] <- paste(dep[i], "model", sep = "_")
}

mod_compare$survival_model
mod_compare$logsizet1_model
mod_compare$flower_model
mod_compare$number_of_seeds_model
