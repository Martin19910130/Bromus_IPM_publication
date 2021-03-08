##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Lambda calculation, Bootstrap
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()
options(stringsAsFactors = F)

require(dplyr)
require(ggplot2)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###         Load functions
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  1. function, size range
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Transforms all values below/above limits in min/max size
## this is to deal with unintentional eviction
x_range <- function(x, pars){
  pmin(pmax(x, pars$L), pars$U)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  2. inverse logit function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create the function to apply the inverse logit
inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  3. Growth function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Growth (transition) from size x to size y
gxy <- function(y,x,pars){
  xb <- x_range(x, pars)
  ## returns a *probability density distribution* for each x value
  return(dnorm(y, 
               mean = pars$grow_b0 + pars$grow_b1*xb, 
               sd   = pars$grow_sd) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  4. Survival at size x function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sx<-function(x,pars){
  xb <- x_range(x, pars)
  ## survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + pars$surv_b1 * xb) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  5. transition function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## transition: Survival * growth
pxy<-function(y,x,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(y,xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 6. Production of flowers
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fx1<-function(x,pars){
  ## probability of flowering of each x size class 
  return( inv_logit(pars$flowprob_b0 + 
                      pars$flowprob_b1 * x) )
}

fx2<-function(x,pars){
  ## number of seeds per flowering plant based on size class
  return( exp( pars$seednum_b0 + 
                 pars$seednum_b1 * x) )
}

## put them together- probability of flowering * seeds per flowering plant
fx3<-function(x,pars){
  xb <- x_range(x, pars)
  return( fx1(xb,pars) * fx2(xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 7. Size distribution of new plants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
newpl <- function(y, pars){
  dnorm( x    = y,
         mean = pars$new_size_mean,
         sd   = pars$new_size_sd )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Unintentional eviction using Merow's MEE article suggestion -----------
##             Kernel function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ker_mee <- function(pars){
  
  # set up necessary variables
  h   <- (pars$U-pars$L)/pars$n             #Bin size
  b   <- pars$L+c(0:pars$n)*h               #Lower boundaries of bins 
  y   <- 0.5*(b[1:pars$n]+b[2:(pars$n+1)])  #Bins' midpoints
  n   <- pars$n
  
  # Fertility matrix
  Fmat            <- matrix(0,(n+1),(n+1))
  
  # April seedlings go in top row
  Fmat[1,2:(n+1)] <- fx3(y,pars) * pars$aprilSL_per_fl
  
  
  # November seedlings go into 200 by 200 plant matrix
  Fmat[2:(n+1),2:(n+1)] <- fx3(y,pars) * pars$novSL_per_fl * pars$SL_surv * newpl(y,pars) * h 
  
  # Growth/survival transition matrix
  Smat            <- matrix(0,(n+1),(n+1))
  Gmat            <- matrix(0,(n+1),(n+1))
  Tmat            <- matrix(0,(n+1),(n+1))
  
  # survival matrix
  Smat[2:(n+1),1] <- pars$SL_surv
  for(i in 2:(n+1)){ 
    Smat[2:(n+1),i] <- sx(y[i-1], pars)
  }
  
  # Growth transitions among cts sizes
  Gmat[2:(n+1),
       2:(n+1)]   <- outer(y,y,gxy,pars) * h
  
  # seedlings stage distribution
  Gmat[2:(n+1),1] <- newpl(y,pars) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat            <- Gmat
  
  # fix eviction for seedligns
  Gmat[2,1]  <- Gmat[2,1]+1-sum(Gmat[2:(n+1),1])
  Tmat[,1]   <- Gmat[,1]*Smat[,1]
  
  # fix eviction of large adults
  for(i in 2:(n+1)){
    Gmat[n,i]  <- Gmat[n,i]+1-sum(Gmat[2:(n+1),i])
    Tmat[,i]   <- Gmat[,i]*Smat[,i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrix
  k_yx            <- Fmat + Tmat     
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y ))
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## bootstrap function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot_lam <- function(ii)
{
  # set the random number generator
  set.seed(ii)
  
  # sample the data
  data <- sample_n(demo_dat, nrow(demo_dat), replace = T)
  data_s <- sample_n(seed_dat, nrow(seed_dat), replace = T)
  
  # calculate the mean seeds 
  mean_seed <- mean(data_s$number_of_seeds, na.rm = T)
  
  # add the mean seeds to the demography data
  data$mean_seed <- ifelse(data$sample_year == 2018, mean_seed, NA)
  data$seed_per_ind <- data$mean_seed * data$number_of_flowers
  
  #--------------------------------------
  # survival
  #--------------------------------------
  sr_mod <- glm(survival ~ logsizet0, data = data, family = binomial())
  
  #--------------------------------------
  # growth
  #--------------------------------------
  gr_mod <- lm(logsizet1 ~ logsizet0, data = data)
  
  #--------------------------------------
  # flower probability
  #--------------------------------------
  fl_mod <- glm(flower ~ logsizet0, data = data, family = binomial())
  
  #--------------------------------------
  # number of seeds per flowering plant
  #--------------------------------------
  data2 <- subset(data, flower == "1")
  sefl_mod <- glm(round(seed_per_ind) ~ logsizet0, data = data2, family = poisson)
  
  #--------------------------------------
  # seedling per seed
  #--------------------------------------
  data3 <- subset(data, subplot <= 3)
  
  # get the seedling count per subplot
  Sl_per_subplot_apr18 <- unique(data3[, c("plot", "subplot", "apr18SL")])
  Sl_per_subplot_nov18 <- unique(data3[, c("plot", "subplot", "nov18SL")])
  Sl_per_subplot_apr19 <- unique(data3[, c("plot", "subplot", "apr19SL")])
  
  # NA's are actually zero seedlings
  Sl_per_subplot_apr18$apr18SL[is.na(Sl_per_subplot_apr18$apr18SL)] <- 0 
  Sl_per_subplot_nov18$nov18SL[is.na(Sl_per_subplot_nov18$nov18SL)] <- 0 
  Sl_per_subplot_apr19$apr19SL[is.na(Sl_per_subplot_apr19$apr19SL)] <- 0
  
  # calculate the sum of seeds per subplot
  seedling_per_seed <- aggregate(data3$seed_per_ind, by = list(data3$plot,
                                                               data3$subplot), FUN = sum, na.rm = T)
  
  Sl_per_subplot_apr18 <- Sl_per_subplot_apr18[order(Sl_per_subplot_apr18$plot,
                                                     Sl_per_subplot_apr18$subplot),]
  Sl_per_subplot_nov18 <- Sl_per_subplot_nov18[order(Sl_per_subplot_nov18$plot,
                                                     Sl_per_subplot_nov18$subplot),]
  Sl_per_subplot_apr19 <- Sl_per_subplot_apr19[order(Sl_per_subplot_apr19$plot, 
                                                     Sl_per_subplot_apr19$subplot),]
  seedling_per_seed <- seedling_per_seed[order(seedling_per_seed$Group.1, 
                                               seedling_per_seed$Group.2),]
  
  Sl_per_subplot_nov18$seeds <- ifelse(Sl_per_subplot_nov18$plot == seedling_per_seed$Group.1 & 
                                         Sl_per_subplot_nov18$subplot == seedling_per_seed$Group.2, 
                                       seedling_per_seed$x, NA)
  Sl_per_subplot_nov18 <- subset(Sl_per_subplot_nov18, seeds > 0)
  
  Sl_per_subplot_apr19$seeds <- ifelse(Sl_per_subplot_apr19$plot == seedling_per_seed$Group.1 & 
                                         Sl_per_subplot_apr19$subplot == seedling_per_seed$Group.2, 
                                       seedling_per_seed$x, NA)
  Sl_per_subplot_apr19 <- subset(Sl_per_subplot_apr19, seeds > 0)
  
  # calculate the turnover of seeds to seedlings
  Sl_per_subplot_nov18$seed_SL_nov18 <- Sl_per_subplot_nov18$nov18SL / Sl_per_subplot_nov18$seeds
  Sl_per_subplot_apr19$seed_SL_apr19 <- Sl_per_subplot_apr19$apr19SL / Sl_per_subplot_apr19$seeds
  seedling_per_seed[seedling_per_seed == Inf] <- NA
  seedling_per_seed[is.na(seedling_per_seed)] <- NA
  
  #--------------------------------
  # Seedling survival
  #--------------------------------
  sl_surv <- sum(data3$new_plant, na.rm = T) / (sum(Sl_per_subplot_apr18$apr18SL) + sum(Sl_per_subplot_nov18$nov18SL))
  
  #--------------------------------
  # new Individual size
  #--------------------------------
  data_new_ind <- subset(data, new_plant == 1)
  new_size_mean <- mean(data_new_ind$logsizet1, na.rm = T)
  new_size_sd <- sd(data_new_ind$logsizet1, na.rm = T)
  
  # store the parameters of the models
  pars <- list(  surv_b0       = summary(sr_mod)$coefficients[1,1],
                 surv_b1       = summary(sr_mod)$coefficients[2,1],
                 grow_b0       = summary(gr_mod)$coefficients[1,1],
                 grow_b1       = summary(gr_mod)$coefficients[2,1],
                 grow_sd       = summary(gr_mod)$sigma,
                 seednum_b0    = summary(sefl_mod)$coefficients[1,1],
                 seednum_b1    = summary(sefl_mod)$coefficients[2,1],
                 flowprob_b0   = summary(fl_mod)$coefficients[1,1],
                 flowprob_b1   = summary(fl_mod)$coefficients[2,1],
                 new_size_mean  = new_size_mean,
                 new_size_sd    = new_size_sd,
                 novSL_per_fl   = mean(Sl_per_subplot_nov18$seed_SL_nov18, na.rm = T),
                 aprilSL_per_fl = mean(Sl_per_subplot_apr19$seed_SL_apr19, na.rm = T),
                 SL_surv        = sl_surv,
                 L = min( c(data$logsizet0,
                            data$logsizet1),
                          na.rm=T),
                 U = max( c(data$logsizet0,
                            data$logsizet1),
                          na.rm=T),
                 L_obs = min( c(data$logsizet0,
                                data$logsizet1),
                              na.rm=T),
                 U_obs = max( c(data$logsizet0,
                                data$logsizet1),
                              na.rm=T),
                 n = 200
  )
  
  # calculate Lambda
  Re(eigen(ker_mee(pars)$k_yx)$value[1])
  
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####           Read data
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_d <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Demography.csv")
dat_s <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Seed_data.csv")

## Delete the new plants which grow bigger then 30
dat_d <- dat_d[-(which(dat_d$sizet1 >= 30 & dat_d$new_plant == 1)),]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            Ambient mowing
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lamb_list <- c()
for(i in 1:length(unique(dat_d$treatment)))
{
  demo_dat <- subset(dat_d, treatment == unique(dat_d$treatment)[i])
  seed_dat <- subset(dat_s, treatment == unique(dat_s$treatment)[i])
  
  lamb_list[i] <- list(sapply(1:1000, boot_lam))
  names(lamb_list)[i] <- unique(dat_d$treatment)[i]
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Plot mean Lambdas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create data frame with mean and Confidence intervalls
Results <- data.frame(lambda_mean = c(mean(lamb_list$ambient_mowing), 
                                      mean(lamb_list$ambient_grazing), 
                                      mean(lamb_list$future_mowing), 
                                      mean(lamb_list$future_grazing)),
                      lambda_U = c(quantile(lamb_list$ambient_mowing, 
                                            probs = c(0.025, 0.975))[2],
                                   quantile(lamb_list$ambient_grazing, 
                                            probs = c(0.025, 0.975))[2],
                                   quantile(lamb_list$future_mowing, 
                                            probs = c(0.025, 0.975))[2],
                                   quantile(lamb_list$future_grazing, 
                                            probs = c(0.025, 0.975))[2]),
                      lambda_L = c(quantile(lamb_list$ambient_mowing, 
                                            probs = c(0.025, 0.975))[1],
                                   quantile(lamb_list$ambient_grazing, 
                                            probs = c(0.025, 0.975))[1],
                                   quantile(lamb_list$future_mowing,
                                            probs = c(0.025, 0.975))[1],
                                   quantile(lamb_list$future_grazing,
                                            probs = c(0.025, 0.975))[1]),
                      treatment = c("ambient_mowing", "ambient_grazing", "future_mowing",
                                    "future_grazing"),
                      landuse = c("mowing", "grazing", "mowing", "grazing"),
                      climate = c("ambient", "ambient", "future", "future"))

## Print and look at results
Results

## Choose color palette
myPalette <- c("#0072B2", "#D55E00")

ggplot(Results, aes(x = climate, y = lambda_mean, color = climate, shape = landuse)) + theme_classic() +
  geom_line(aes(group = landuse , linetype = landuse), col = "black", position = position_dodge(.08)) + 
  scale_shape_manual(values=c(16, 17)) +
  scale_colour_manual(values=myPalette) +
  scale_linetype_manual(values = c("grazing" = "dashed" , "mowing" = "solid")) + 
  geom_point( size = 3, position = position_dodge(width = 0.08)) +
  geom_errorbar(aes(ymax = lambda_U, ymin = lambda_L), width = .17, size = 1.2,
                position = position_dodge(0.08)) +
  theme(text = element_text(size = 20), 
        legend.position = c(0.87, 0.89), 
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())  +
  ylab("Asymptotic population growth rate (Î»)") + xlab("Climate treatment") + 
  guides(color = F)

ggsave("C:\\Users/ma22buky/Documents/Julia_Paper/lambda.pdf", 
       plot = last_plot(), 
       device = cairo_pdf, 
       dpi = 1200, 
       width = 18,
       height = 14, 
       units = "cm")


savers=matrix(0, 1000, 1)

for (i in 1:1000) {
  
  savers[i]=lamb_amb_mow[sample(1:1000, 1)]-lamb_amb_gra[sample(1:1000, 1)]
  
}

sortedsaved=sort(savers)

library(readxl)

"C:\\Users\\Julia\\Desktop\\Uni\\Masterarbeit\\Analyse"
randtestdata <- read_excel("C:\\Users\\Julia\\Desktop\\Uni\\Publication\\Daten\\randtestdata.xlsx")



ambgr <- randtestdata$amb_graz

ambmow <- randtestdata$amb_mow



savers=matrix(0, 1000, 1)

for (i in 1:1000) {
  
  savers[i]=ambmow[sample(1:1000, 1)]-ambgr[sample(1:1000, 1)]
  
}



sortedsaved=sort(savers)
