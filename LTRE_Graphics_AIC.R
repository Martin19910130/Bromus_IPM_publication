####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####          Figures, AIC, LTRE
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
gc()
options(stringsAsFactors = F)

library(ggplot2)
library(dplyr)
library(wiqid)
library(gridExtra)
library(grid)
library(patchwork)

## change the resolution you want your figures stored in
res <- 1500

rbPalette <- c("#0072B2", "#D55E00")

df_binned_prop <- function(df, n_bins, siz_var, rsp_var)
{
  
  size_var <- deparse( substitute(siz_var) )
  resp_var <- deparse( substitute(rsp_var) )
  
  min_size <- min(df[,size_var], na.rm = T) - 0.000001
  max_size <- max(df[,size_var], na.rm = T) + 0.000001
  
  # binned survival probabilities
  h    <- (max_size - min_size) / n_bins
  lwr  <- min_size + (h*c(0:(n_bins-1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2*h)
  
  binned_prop <- function(lwr_x, upr_x, response){
    
    id  <- which(df[,size_var] > lwr_x & df[,size_var] < upr_x)
    tmp <- df[id,]
    
    if( response == 'prob' ){   return( sum(tmp[,resp_var],na.rm=T) / nrow(tmp) ) }
    if( response == 'n_size' ){ return( nrow(tmp) ) }
    
  }
  
  y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
  x_binned <- mid
  y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
  
  data.frame( x = x_binned,
              y = y_binned )
}

rotate <- function(x) t(apply(x, 2, rev))

mylegend <- g_legend <- function(a.gplot)
{
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###           read data
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Demography.csv")
seed_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Seed_data.csv")

## Delete the new plants which grew bigger then 30
demo_data <- demo_data[-(which(demo_data$sizet1 >= 30 & demo_data$new_plant == 1)),]

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###     Add seed data to demography data
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seed_trea <- aggregate(seed_data$number_of_seeds, by = list(seed_data$treatment), FUN = mean)
row.names(seed_trea) <- seed_trea[,1]

demo_data[, "mean_seeds"] <- ifelse(demo_data$sample_year == 2018 & demo_data$treatment == "ambient_grazing", 
                                    seed_trea["ambient_grazing", "x"], NA)
demo_data[, "mean_seeds"] <- ifelse(demo_data$sample_year == 2018 & demo_data$treatment == "ambient_mowing",
                                    seed_trea["ambient_mowing", "x"], demo_data$mean_seeds)
demo_data[, "mean_seeds"] <- ifelse(demo_data$sample_year == 2018 & demo_data$treatment == "future_grazing",
                                    seed_trea["future_grazing", "x"], demo_data$mean_seeds)
demo_data[, "mean_seeds"] <- ifelse(demo_data$sample_year == 2018 & demo_data$treatment == "future_mowing", 
                                    seed_trea["future_mowing", "x"], demo_data$mean_seeds)

demo_data$number_of_seeds <- demo_data$number_of_flowers * round(demo_data$mean_seeds)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    Subset data into 4 different treatments for the bin function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
treats <- unique(demo_data$treatment)
demo_treats <- c()

for(i in 1:length(treats))
{
  demo_treats[i] <- list(assign(treats[i], subset(demo_data, treatment == treats[i])))
  
  rm(list = treats[i])
  names(demo_treats)[i] <- treats[i]
}

## get bins for survival
df_bins_surv <- c()
for(i in 1:length(demo_treats))
{
  df_bins_surv[i] <- list(df_binned_prop(demo_treats[[i]], 10, logsizet0, survival))
  df_bins_surv[[i]]$treatment <- names(demo_treats)[i]
  df_bins_surv[[i]]$climate <- demo_treats[[i]][1, "climate"]
  df_bins_surv[[i]]$management <- demo_treats[[i]][1, "management"]
}

df_bins_surv <- do.call(rbind.data.frame, df_bins_surv)

## get bins for flower probability 
df_bins_flow <- c()
for(i in 1:length(demo_treats))
{
  df_bins_flow[i] <- list(df_binned_prop(demo_treats[[i]], 10, logsizet0, flower))
  df_bins_flow[[i]]$climate <- demo_treats[[i]][1, "climate"]
  df_bins_flow[[i]]$management <- demo_treats[[i]][1, "management"]
}
df_bins_flow <- do.call(rbind.data.frame, df_bins_flow)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###                 Graphics
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create a treatment column
df_bins_surv$treatment <- paste(df_bins_surv$climate, df_bins_surv$management)
df_bins_flow$treatment <- paste(df_bins_flow$climate, df_bins_flow$management)
demo_data$treatment <- paste(demo_data$climate, demo_data$management)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        Survival plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
surv_pl <- ggplot(df_bins_surv, aes(x = x, y = y, shape = treatment, color = treatment, 
                                    linetype = treatment)) + 
  geom_point(size = 1.7) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              mapping = aes(x = logsizet0, y = survival), se = F, data = demo_data, size = 0.5) +
  scale_linetype_manual(values = c("ambient grazing" = "dashed", 
                                   "future grazing" = "dashed",
                                   "ambient mowing" = "solid", 
                                   "future mowing" = "solid")) + 
  scale_color_manual(values = c("ambient grazing" = "#0072B2", 
                                "ambient mowing" = "#0072B2",
                                "future grazing" = "#D55E00",
                                "future mowing" = "#D55E00")) + 
  scale_shape_manual(values = c("ambient grazing" = 1,
                                "future grazing" = 1, 
                                "ambient mowing" = 17,
                                "future mowing" = 17)) + 
  ylab("Plant survival") + xlab("") + theme_classic() +
  scale_y_continuous(breaks = c(0,1), limits = c(0,1)) + 
  scale_x_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) +  
  ggtitle("(A) Plant survival") + 
  theme(legend.position = "", plot.title = element_text(size = 10), 
        aspect.ratio = 0.5) + 
  guides(color = guide_legend(""), linetype = guide_legend(""), shape = guide_legend(""))

surv_pl

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Growth plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grow_pl <- ggplot(demo_data, aes(x = logsizet0, y = logsizet1, shape = treatment, color = treatment, 
                                 linetype = treatment)) + geom_point(size = 1.7) +   
  scale_linetype_manual(values = c("ambient grazing" = "dashed", 
                                   "future grazing" = "dashed",
                                   "ambient mowing" = "solid", 
                                   "future mowing" = "solid")) + 
  scale_color_manual(values = c("ambient grazing" = "#0072B2", 
                                "ambient mowing" = "#0072B2",
                                "future grazing" = "#D55E00",
                                "future mowing" = "#D55E00")) + 
  scale_shape_manual(values = c("ambient grazing" = 1,
                                "future grazing" = 1, 
                                "ambient mowing" = 17,
                                "future mowing" = 17)) + theme_classic() + 
  geom_smooth(method = "lm", se = F, size = 0.5) + 
  ylab(expression(paste("log (size ", italic("t"), " + 1)"))) + xlab("") +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        aspect.ratio = 0.5) + 
  ggtitle("(B) Plant growth") +
  scale_x_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) +
  scale_y_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) + 
  guides(color = guide_legend(""), linetype = guide_legend(""), shape = guide_legend(""))

grow_pl
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      Flower probability plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flpr_pl <- ggplot(df_bins_flow, aes(x = x, y = y, shape = treatment, color = treatment, 
                                    linetype = treatment)) +
  geom_point(size = 1.7) +  
  scale_linetype_manual(values = c("ambient grazing" = "dashed", 
                                   "future grazing" = "dashed",
                                   "ambient mowing" = "solid", 
                                   "future mowing" = "solid")) + 
  scale_color_manual(values = c("ambient grazing" = "#0072B2", 
                                "ambient mowing" = "#0072B2",
                                "future grazing" = "#D55E00",
                                "future mowing" = "#D55E00")) + 
  scale_shape_manual(values = c("ambient grazing" = 1,
                                "future grazing" = 1, 
                                "ambient mowing" = 17,
                                "future mowing" = 17)) + theme_classic() + 
  geom_smooth(data = demo_data, mapping = aes(x = logsizet0, y = flower),
              method = "glm" , method.args = list(family = "binomial"), se = F, size = .5) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) + 
  ylab("Flower probability") + xlab(expression(paste("log (size ", italic("t"), ")"))) +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        aspect.ratio = 0.5) +
  labs(color = "Treatment") + ggtitle("(C) Reproduction probability") + 
  scale_x_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) + 
  guides(color = guide_legend(""), linetype = guide_legend(""), shape = guide_legend(""))

flpr_pl

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      number of seeds per flower plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_data2 <- subset(demo_data, flower == 1)

nrse_pl <- ggplot(demo_data2, aes(x = logsizet0, y = number_of_seeds, color = treatment, 
                                  shape = treatment, linetype = treatment)) + 
  geom_point() + theme_classic() + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = F, size = .5) + 
  ylab("Seeds per plant")  + 
  scale_linetype_manual(values = c("ambient grazing" = "dashed", 
                                   "future grazing" = "dashed",
                                   "ambient mowing" = "solid", 
                                   "future mowing" = "solid")) + 
  scale_color_manual(values = c("ambient grazing" = "#0072B2", 
                                "ambient mowing" = "#0072B2",
                                "future grazing" = "#D55E00",
                                "future mowing" = "#D55E00")) + 
  scale_shape_manual(values = c("ambient grazing" = 1,
                                "future grazing" = 1, 
                                "ambient mowing" = 17,
                                "future mowing" = 17)) +
  xlab(expression(paste("log (size ", italic("t"), ")"))) + 
  theme(legend.position = "none", plot.title = element_text(size = 10),
        aspect.ratio = 0.5) + 
  ggtitle("(D) Seeds per reproductive plant") + labs(color = "Treatment") + 
  scale_y_continuous(breaks=c(0, 1000, 2000)) +
  scale_x_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) + 
  guides(color = guide_legend(""), linetype = guide_legend(""), shape = guide_legend(""))

nrse_pl
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      Seedling per seed plot (only on first 3 subplots)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_data3 <- subset(demo_data, subplot <= 3)
dat_apr18_sl <- unique(demo_data3[, c("plot", "subplot", "climate", "management", "apr18SL")])
dat_nov18_sl <- unique(demo_data3[, c("plot", "subplot", "climate", "management", "nov18SL")])
dat_apr19_sl <- unique(demo_data3[, c("plot", "subplot", "climate", "management", "apr19SL")])

## NA seedlings == no seedlings
dat_apr18_sl$apr18SL[is.na(dat_apr18_sl$apr18SL)] <- 0 
dat_nov18_sl$nov18SL[is.na(dat_nov18_sl$nov18SL)] <- 0 
dat_apr19_sl$apr19SL[is.na(dat_nov18_sl$apr19SL)] <- 0

sl_per_plot_apr18 <- aggregate(dat_apr18_sl$apr18SL, 
                               by = list(dat_apr18_sl$plot, dat_apr18_sl$subplot, dat_apr18_sl$climate, 
                                         dat_apr18_sl$management), 
                               FUN = sum, na.rm = T)
colnames(sl_per_plot_apr18) <-  c("plot", "subplot", "climate", "management", "seedling_apr18")

sl_per_plot_nov18 <- aggregate(dat_nov18_sl$nov18SL, 
                               by = list(dat_nov18_sl$plot, dat_nov18_sl$subplot, dat_nov18_sl$climate, 
                                         dat_nov18_sl$management), 
                               FUN = sum, na.rm = T)
colnames(sl_per_plot_nov18) <-  c("plot", "subplot","climate", "management", "seedling")

sl_per_plot_apr19 <- aggregate(dat_apr19_sl$apr19SL, 
                               by = list(dat_apr19_sl$plot, dat_apr19_sl$subplot, 
                                         dat_apr19_sl$climate, dat_apr19_sl$management), 
                               FUN = sum, na.rm = T)
colnames(sl_per_plot_apr19) <-  c("plot", "subplot","climate", "management", "seedling")

se_per_plot <- aggregate(demo_data3$number_of_seeds, 
                         by = list(demo_data3$plot, demo_data3$subplot, 
                                   demo_data3$climate, demo_data3$management), 
                         FUN = sum, na.rm = T)
colnames(se_per_plot) <-  c("plot", "subplot","climate", "management", "seed")

nov_18 <- right_join(se_per_plot, sl_per_plot_nov18, 
                     by = c("plot", "subplot", "climate", "management"))
apr_19 <- right_join(se_per_plot, sl_per_plot_apr19, 
                     by = c("plot", "subplot", "climate", "management"))

nov_18$sl_per_seed <- nov_18$seedling/nov_18$seed 
nov_18$sl_per_seed[nov_18$sl_per_seed == Inf] <- NA
nov_18[is.na(nov_18)] <- NA

apr_19$sl_per_seed <- apr_19$seedling/apr_19$seed
apr_19$sl_per_seed[apr_19$sl_per_seed == Inf] <- NA
apr_19[is.na(apr_19)] <- NA

## nov 18 seeds to seedling
## use 15 for the se because there are 15 samples per treatment (5*3, cause of subplot)
mean_sl_seed <- aggregate(nov_18$sl_per_seed, by = list(nov_18$climate, nov_18$management), 
                          FUN = mean, na.rm = T)
mean_sl_seed$sd <- aggregate(nov_18$sl_per_seed, by = list(nov_18$climate, nov_18$management), 
                             FUN = sd, na.rm = T)[,"x"]
colnames(mean_sl_seed) <- c("climate", "management", "mean_sl", "sd")
mean_sl_seed$se <- mean_sl_seed$sd / sqrt(length(mean_sl_seed))

## Plot nov 18 seed to seedling
sl_seed_nov18 <- ggplot(mean_sl_seed, aes(x = management, y = mean_sl, fill = climate)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.94)) + 
  geom_errorbar(aes(ymin = mean_sl_seed$mean_sl - mean_sl_seed$se, 
                    ymax = mean_sl_seed$mean_sl + mean_sl_seed$se), 
                width = 0.01, position = position_dodge(0.94)) + 
  geom_point(nov_18, mapping = aes(x = management, y = sl_per_seed, 
                                   shape = management), 
             position = position_dodge(0.94)) + 
  scale_fill_manual(values = rbPalette) +
  theme_classic() + ggtitle("(E) Fall recruitment") +
  theme(legend.position = "none", axis.title.x=element_blank(), 
        plot.title = element_text(size = 10), aspect.ratio = 0.5) +
  ylab("Seedlings / Quadrat") + xlab("Treatment combination") +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05), limits = c(-0.01, 0.25))

## apr 19 seeds to seedling
mean_sl_seed_apr19 <- aggregate(apr_19$sl_per_seed, by = list(apr_19$climate, apr_19$management),
                                FUN = mean, na.rm = T)
mean_sl_seed_apr19$sd <- aggregate(apr_19$sl_per_seed, by = list(apr_19$climate, nov_18$management),
                                   FUN = sd, na.rm = T)[, "x"]
colnames(mean_sl_seed_apr19) <- c("climate", "management", "mean_sl", "sd")
mean_sl_seed_apr19$se <- mean_sl_seed_apr19$sd / sqrt(15)

## Plot seeds to seedling apr 19 
sl_seed_apr19 <- ggplot(mean_sl_seed_apr19, aes(x = management, y = mean_sl, fill = climate)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.94)) + 
  geom_errorbar(aes(ymin = mean_sl_seed_apr19$mean_sl - mean_sl_seed_apr19$se, 
                    ymax = mean_sl_seed_apr19$mean_sl + mean_sl_seed_apr19$se), 
                width = 0.01, position = position_dodge(0.94)) + 
  geom_point(apr_19, mapping = aes(x = management, y = sl_per_seed, 
                                   shape = management), 
             position = position_dodge(0.94)) + 
  scale_fill_manual(values = rbPalette) +
  theme_classic() + ggtitle("(F) Spring recruitment") +
  theme(legend.position = "none", axis.title.x=element_blank(),
        plot.title = element_text(size = 10), aspect.ratio = 0.5) +
  ylab("Seedlings / Quadrat") + xlab("Treatment combination") +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05), limits = c(-0.01, 0.25))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Seedling survival
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_apr18_sl
dat_nov18_sl
new_pl <- aggregate(demo_data3$new_plant, 
                    by = list(demo_data3$plot, demo_data3$subplot, 
                              demo_data3$climate, demo_data3$management), 
                    FUN = sum, na.rm = T)

colnames(new_pl) <- c("plot", "subplot", "climate", "management", "new_pl")

sl_apr18_nov18 <- right_join(sl_per_plot_nov18, sl_per_plot_apr18, by = c("plot", "subplot", "climate", "management"))
sl_apr18_nov18$sl_sum <- sl_apr18_nov18$seedling + sl_apr18_nov18$seedling_apr18

sl_surv <- right_join(new_pl, sl_apr18_nov18, by = c("plot", "subplot", "climate", "management"))

sl_surv$sl_surv <- sl_surv$new_pl/sl_surv$sl_sum
sl_surv$sl_surv <- ifelse(sl_surv$sl_surv > 1, 1, sl_surv$sl_surv)
sl_surv[sl_surv == Inf] <- NA
sl_surv[is.na(sl_surv)] <- NA

## use 15 for the se because there are 15 samples per treatment (5*3, cause of subplot)
sl_surv$treatment <- paste(sl_surv$climate, sl_surv$management, sep = "_")
table(sl_surv$treatment)

mean_sl_surv <- aggregate(sl_surv$sl_surv, by = list(sl_surv$climate, sl_surv$management), 
                          FUN = mean, na.rm = T)
mean_sl_surv$sd <- aggregate(sl_surv$sl_surv, by = list(sl_surv$climate, sl_surv$management),
                             FUN = sd, na.rm = T)[,"x"]

mean_sl_surv$se <- mean_sl_surv$sd/sqrt(15)
colnames(mean_sl_surv) <- c("climate", "management", "mean", "sd", "se")

sl_surv_pl <- ggplot(mean_sl_surv, aes(x = management, y = mean, fill = climate)) +  
  geom_bar(stat = "identity", position = position_dodge(0.94))  + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.94), 
                width = 0) + 
  geom_point(sl_surv, mapping = aes(y = sl_surv, x = management, shape = management), 
             position = position_dodge(0.94)) +  theme_classic() + 
  scale_fill_manual(values = rbPalette) + scale_y_continuous(breaks=c(0, 0.5, 1), limits = c(0, 1)) +
  theme(legend.position = "none", axis.title.x=element_blank(), 
        plot.title = element_text(size = 10), aspect.ratio = 0.5) +
  ylab("Survival rate") + xlab("Treatment combination") +
  ggtitle("(G) Establishment")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      New plants size distribution 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_dat_new <- subset(demo_data, new_plant == 1)

new_size_pl <- ggplot(demo_dat_new, aes(x = management, y = logsizet1, fill = climate)) +  geom_boxplot() + 
  geom_point(demo_dat_new, mapping = aes(x = management, y = logsizet1, shape = management),
             position = position_dodge(0.75)) + ggtitle("(H) Size distribution of new plants") +  
  scale_fill_manual(values = rbPalette) + theme_classic() +  
  theme(legend.position = "none", axis.title.x=element_blank(), 
        plot.title = element_text(size = 10), aspect.ratio = 0.5) + 
  scale_y_continuous(breaks = c(-1.38629 ,0 , 1.6094, 2.70805, 3.9120, 5.0106), 
                     labels = c(0.25, 1, 5, 15, 50, 150)) + 
  ylab(expression(paste("log (size ", italic("t"), " + 1)"))) 

Fig_2 <- ggpubr::ggarrange(surv_pl, grow_pl, flpr_pl, nrse_pl, sl_seed_nov18,
                           sl_seed_apr19 , sl_surv_pl, new_size_pl, ncol = 2, nrow = 4, common.legend = T,
                           legend = "bottom", align = "hv", widths = 3, heights = 5)
Fig_2

#save the figure as pdf
ggsave("C:\\Users/ma22buky/Documents/Julia_Paper/FigVR_binned.tiff", 
       plot = Fig_2, 
       device = "tiff", 
       dpi = res, 
       width = 18,
       height = 15, 
       units = "cm")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          LTRE, sensitivty, models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LTRE_fun <- function(demo_data)
{
  # sample the data
  data <- demo_data
  data_s <- seed_data
  
  # calculate the mean seeds 
  mean_seed <- mean(data_s$number_of_seeds, na.rm = T)
  
  # add the mean seeds to the demography data
  data$mean_seed <- ifelse(data$sample_year == 2018, mean_seed, NA)
  data$seed_per_ind <- data$mean_seed * data$number_of_flowers
  
  #--------------------------------------
  # survival
  #--------------------------------------
  sr_mod <- glm(survival ~ logsizet0, data = data, family = binomial())
  summary(sr_mod)
  
  #--------------------------------------
  # growth
  #--------------------------------------
  gr_mod <- lm(logsizet1 ~ logsizet0, data = data)
  summary(gr_mod)
  
  #--------------------------------------
  # flower probability
  #--------------------------------------
  fl_mod <- glm(flower ~ logsizet0, data = data, family = binomial())
  summary(fl_mod)
  
  #--------------------------------------
  # number of seeds per flowering plant
  #--------------------------------------
  data2 <- subset(data, flower == "1")
  sefl_mod <- glm(round(seed_per_ind) ~ logsizet0, data = data2, family = 'poisson')
  summary(sefl_mod)
  
  #--------------------------------------
  # seedling per seed
  #--------------------------------------
  data3 <- subset(data, subplot <= 3)
  
  # get the seedling count per subplot
  Sl_per_subplot_apr18 <- unique(data3[, c("plot", "subplot", "apr18SL")])
  Sl_per_subplot_nov18 <- unique(data3[, c("plot", "subplot", "nov18SL")])
  Sl_per_subplot_apr19 <- unique(data3[, c("plot", "subplot", "apr19SL")])
  
  # calculate the sum of seedling per plot
  Sl_per_plot_apr18 <- aggregate(Sl_per_subplot_apr18$apr18SL, 
                                 by = list(Sl_per_subplot_apr18$plot), FUN = sum, na.rm = T)
  Sl_per_plot_nov18 <- aggregate(Sl_per_subplot_nov18$nov18SL, 
                                 by = list(Sl_per_subplot_nov18$plot), FUN = sum, na.rm = T)
  Sl_per_plot_apr19 <- aggregate(Sl_per_subplot_apr19$apr19SL, 
                                 by = list(Sl_per_subplot_apr19$plot), FUN = sum, na.rm = T)
  
  # calculate the sum of seeds per plot
  seedling_per_seed <- aggregate(data3$seed_per_ind, by = list(data3$plot), FUN = sum, na.rm = T)
  
  # calculate the turnover of seeds to seedlings
  seedling_per_seed$seed_SL_nov18 <- Sl_per_plot_nov18$x / seedling_per_seed$x
  seedling_per_seed$seed_SL_apr19 <- Sl_per_plot_apr19$x / seedling_per_seed$x
  seedling_per_seed[seedling_per_seed == Inf] <- NA
  
  #--------------------------------
  # Seedling survival
  #--------------------------------
  sl_surv <- sum(data3$new_plant, na.rm = T) / (sum(Sl_per_plot_apr18$x) + sum(Sl_per_plot_nov18$x))
  
  #--------------------------------
  # new Individual size
  #--------------------------------
  data_new_ind <- subset(data, new_plant == 1)
  new_size_mean <- mean(data_new_ind$logsizet1, na.rm = T)
  new_size_sd <- sd(data_new_ind$logsizet1, na.rm = T)
  
  # store the parameters of the models
  pars <- data.frame(  surv_b0       = coef(sr_mod)[1],
                       surv_b1       = coef(sr_mod)[2],
                       grow_b0       = coef(gr_mod)[1],
                       grow_b1       = coef(gr_mod)[2],
                       grow_sd       = summary(gr_mod)$sigma,
                       seednum_b0    = coef(sefl_mod)[1],
                       seednum_b1    = coef(sefl_mod)[2],
                       flowprob_b0   = coef(fl_mod)[1],
                       flowprob_b1   = coef(fl_mod)[2],
                       new_size_mean  = new_size_mean,
                       new_size_sd    = new_size_sd,
                       novSL_per_fl   = mean(seedling_per_seed$seed_SL_nov18, na.rm = T),
                       aprilSL_per_fl = mean(seedling_per_seed$seed_SL_apr19, na.rm = T),
                       SL_surv        = sl_surv,
                       L = min(c(data$logsizet0,
                                 data$logsizet1),
                               na.rm=T),
                       U = max(c(data$logsizet0,
                                 data$logsizet1),
                               na.rm=T),
                       L_obs = min(c(data$logsizet0,
                                     data$logsizet1),
                                   na.rm=T),
                       U_obs = max(c(data$logsizet0,
                                     data$logsizet1),
                                   na.rm=T),
                       n = 200
  )
  
  # calculate Lambda
  return(pars)
}

## Transforms all values below/above limits in min/max size
#this is to deal with unintentional eviction
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
# Growth (transition) from size x to size y
gxy <- function(y,x,pars){
  xb <- x_range(x, pars)
  # returns a *probability density distribution* for each x value
  return(dnorm(y, 
               mean = pars$grow_b0 + pars$grow_b1*xb, 
               sd   = pars$grow_sd) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  4. Survival at size x function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sx<-function(x,pars){
  xb <- x_range(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + pars$surv_b1 * xb) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  5. transition function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# transition: Survival * growth
pxy<-function(y,x,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(y,xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Production of flowers
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fx1<-function(x,pars){
  # probability of flowering of each x size class 
  return( inv_logit(pars$flowprob_b0 + 
                      pars$flowprob_b1 * x) )
}

fx2<-function(x,pars){
  #number of flowers per flowering plant based on size class
  return( exp( pars$seednum_b0 + 
                 pars$seednum_b1 * x) )
}

#put them together- probability of flowering * flowers per flowering plant
fx3<-function(x,pars){
  xb <- x_range(x, pars)
  return( fx1(xb,pars) * fx2(xb,pars) )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Size distribution of new plants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
newpl <- function(y, pars){
  dnorm( x    = y,
         mean = pars$new_size_mean,
         sd   = pars$new_size_sd )
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Unintentional eviction using Merow's MEE article suggestion -----------
#             Kernel function
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
  
  
  # November seedlings go into 50 by 50 plant matrix
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##              Difference in vital rates
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
demo_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Demography.csv")
seed_data <- read.csv("https://raw.githubusercontent.com/Martin19910130/Bromus_IPM_publication/master/Bro_Seed_data.csv")

## Delete the new plants which grow bigger then 30
demo_data <- demo_data[-(which(demo_data$sizet1 >= 30 & demo_data$new_plant == 1)),]

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###     Add seed data to demography data
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
amb_mow <- subset(demo_data, treatment == "ambient_mowing")
amb_gra <- subset(demo_data, treatment == "ambient_grazing")
fut_mow <- subset(demo_data, treatment == "future_mowing")
fut_gra <- subset(demo_data, treatment == "future_grazing")

par_am <- LTRE_fun(demo_data = amb_mow)
par_ag <- LTRE_fun(demo_data = amb_gra)
par_fm <- LTRE_fun(demo_data = fut_mow)
par_fg <- LTRE_fun(demo_data = fut_gra)

mow_dif <- as.data.frame(t(par_am[, 1:14] - par_fm[, 1:14]))
gra_dif <- as.data.frame(t(par_ag[, 1:14] - par_fg[, 1:14]))
fut_dif <- as.data.frame(t(par_fm[, 1:14] - par_fg[, 1:14]))
amb_dif <- as.data.frame(t(par_am[, 1:14] - par_ag[, 1:14]))

# assign rownames and establish own columns for the variables, the order and the color of the graphs
colnames(mow_dif) <- "values"
mow_dif$pars      <- row.names(mow_dif)
mow_dif$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
mow_dif$mycol     <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

colnames(gra_dif) <- "values"
gra_dif$pars      <- row.names(gra_dif)
gra_dif$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
gra_dif$mycol     <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

colnames(amb_dif) <- "values"
amb_dif$pars      <- row.names(amb_dif)
amb_dif$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
amb_dif$mycol     <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

colnames(fut_dif) <- "values"
fut_dif$pars      <- row.names(fut_dif)
fut_dif$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
fut_dif$mycol     <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

mypalette <- c("#DEDEDE","#BABABA",  "#858585", "#525252", "black")
mylabels_VR <- c("S Int", "S Slope", "G Int", "G Slope", 
                 "G SD", "P Int", "P Slope", "F Int", "F Slope", " \u03B8f", "\u03B8s", "B", "\u03B7", "\u03B7 SD")
mylabels_LTRE  <- c("Survival", "Growth", "Reproduction", "Recruitment", "Establishment")

VR_amb <- ggplot(as.data.frame(amb_dif[1:14, ]), aes(y = values, x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ylim (-3.1, 3.1) +
  ggtitle("a)  Ambient") + 
  ylab("grazing       -      mowing ") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12)) +
  scale_fill_manual(values =  mypalette) 

VR_fut <- ggplot(as.data.frame(fut_dif[1:14, ]), aes(y = values, x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ylim (-3.1, 3.1) +
  ggtitle("b) Future") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_manual(values =  mypalette) 


VR_mow <- ggplot(as.data.frame(mow_dif[1:14, ]), aes(y = values, x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ylim (-3.1, 3.1) +
  ggtitle("c) Mowing") + ylab("ambient      -       future") +  theme_bw() + 
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values =  mypalette)  +
  scale_x_discrete(labels = mylabels_VR)

VR_gra <- ggplot(as.data.frame(gra_dif[1:14, ]), aes(y = values, x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ylim (-3.1, 3.1) +
  ggtitle("d) Grazing") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values =  mypalette)  +
  scale_x_discrete(labels = mylabels_VR)

VR_all <- grid.arrange(VR_amb, VR_fut, VR_mow, VR_gra,
                       ncol = 2, nrow = 2, 
                       layout_matrix = rbind(c(1, 2), c(3, 4)), 
                       widths = c(5.5, 5), heights = c(5, 5.8))

ggsave(filename = "C:\\Users/ma22buky/Documents/Julia_Paper/Fig_Dif_VR.pdf", 
       plot = VR_all, 
       device = "tiff", 
       dpi = res, 
       width = 20,
       height =15, 
       units = "cm")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##           Sensitivity & LTREs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ambient
final_out <- data.frame()
pars_man <- c()

for(i in 1:14)
{
  final_out[i, "parameter"] <- colnames(par_am)[i]  
  
  # mean pars am ag:
  pars_man <- (par_am + par_ag) / 2
  
  # "original" Lambda from combined pars
  final_out[i, "lambda_ori"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
  
  # calculate Lambda from manipulated combined pars
  pars_man[1,i] <- pars_man[1,i] + 0.001
  final_out[i, "lambda_fal"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
}
final_out

# calculate sensitivity 
final_out[, "sensitivity"] <- as.numeric((final_out$lambda_fal - final_out$lambda_ori)/0.001)

# calculate LTRE values
final_out[, "LTRE"] <- as.numeric(final_out$sensitivity * amb_dif[1:14, "values"])

# establish order and clour for plots
final_out$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
final_out$mycol       <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

# for ecological interpretation summerize vital rates into life cycle stages
final_out5 <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Sl survival"),
                         LTRE = c(final_out[1,5] + final_out[2,5], 
                                  final_out[3,5] + final_out[4,5] + final_out[5,5], 
                                  final_out[6,5] + final_out[7,5] + final_out[8,5] + final_out[9,5], 
                                  final_out[12,5] + final_out[13,5], 
                                  final_out[10,5] + final_out[11,5] + final_out[14,5]),
                         myorder = c("1", "2", "3", "4", "5"), 
                         mycol = as.factor(c(1:5)))

final_out5$LTRE_scaled = as.numeric(final_out5$LTRE / sum(abs(final_out5$LTRE)))

# Create Sensitivity figure
sensitivity_amb <- ggplot(as.data.frame(final_out), aes(y = sensitivity, x = myorder, fill = mycol)) + geom_bar(stat = "identity") + 
  ylim (0, 70) +
  ggtitle("(A) Mowing ambient * grazing ambient") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +  
  scale_x_discrete(labels = mylabels_VR) +
  scale_fill_manual(values =  mypalette) + 
  labs(y = expression(paste("Sensitivity of ", lambda)))

# Create scaled LTRE figure for life cycle stages
LTRE_amb <- ggplot(as.data.frame(final_out5), aes(y = LTRE_scaled , x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ggtitle("(A) Ambient") + 
  ylim (-0.82, 0.82) +
  labs(y = expression(paste(Delta, lambda, " (mowing - grazing)"))) +
  theme_bw() + theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12)) +
  scale_fill_manual(values =  mypalette) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete(labels = mylabels_LTRE)


## future
final_out <- data.frame()
pars_man <- c()

for(i in 1:14)
{
  final_out[i, "parameter"] <- colnames(par_fm)[i]  
  
  # mean pars fm fg:
  pars_man <- (par_fm + par_fg) / 2
  
  # "original" Lambda from combined pars
  final_out[i, "lambda_ori"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
  
  # calculate Lambda from manipulated combined pars
  pars_man[1,i] <- pars_man[1,i] + 0.001
  final_out[i, "lambda_fal"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
}
final_out

# calculate sensitivity 
final_out[, "sensitivity"] <- as.numeric((final_out$lambda_fal - final_out$lambda_ori)/0.001)

# calculate LTRE values
final_out[, "LTRE"] <- as.numeric(final_out$sensitivity * fut_dif[1:14, "values"])

# establish order and clour for plots
final_out$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
final_out$mycol       <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

# for ecological interpretation summerize vital rates into life cycle stages
final_out5 <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Sl survival"),
                         LTRE = c(final_out[1,5] + final_out[2,5], 
                                  final_out[3,5] + final_out[4,5] + final_out[5,5], 
                                  final_out[6,5] + final_out[7,5] + final_out[8,5] + final_out[9,5], 
                                  final_out[12,5] + final_out[13,5], 
                                  final_out[10,5] + final_out[11,5] + final_out[14,5]),
                         myorder = c("1", "2", "3", "4", "5"), 
                         mycol = as.factor(c(1:5)))

final_out5$LTRE_scaled = as.numeric(final_out5$LTRE / sum(abs(final_out5$LTRE)))

# Create Sensitivity figure
sensitivity_fut <- ggplot(as.data.frame(final_out), aes(y = sensitivity, x = myorder, fill = mycol)) + geom_bar(stat = "identity") + 
  ylim (0, 70) +
  ggtitle("(B) Grazing future * mowing future ") + 
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank()) +  
  scale_x_discrete(labels = mylabels_VR) +
  scale_fill_manual(values =  mypalette) 

# Create scaled LTRE figure for life cycle stages
LTRE_fut <- ggplot(as.data.frame(final_out5), aes(y = LTRE_scaled , x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ggtitle("(B) Future") + 
  ylim (-0.82, 0.82) +
  theme_bw() + theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_fill_manual(values =  mypalette) + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_x_discrete(labels = mylabels_LTRE)


## mowing
final_out <- data.frame()
pars_man <- c()

for(i in 1:14)
{
  final_out[i, "parameter"] <- colnames(par_fm)[i]  
  
  # mean pars am fm:
  pars_man <- (par_am + par_fm) / 2
  
  # "original" Lambda from combined pars
  final_out[i, "lambda_ori"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
  
  # calculate Lambda from manipulated combined pars
  pars_man[1, i] <- pars_man[1,i] + 0.001
  final_out[i, "lambda_fal"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
}
final_out

# calculate sensitivity 
final_out[, "sensitivity"] <- as.numeric((final_out$lambda_fal - final_out$lambda_ori)/0.001)

# calculate LTRE values
final_out[, "LTRE"] <- as.numeric(final_out$sensitivity * mow_dif[1:14, "values"])

# establish order and clour for plots
final_out$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
final_out$mycol       <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

# for ecological interpretation summerize vital rates into life cycle stages
final_out5 <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Sl survival"),
                         LTRE = c(final_out[1,5] + final_out[2,5], 
                                  final_out[3,5] + final_out[4,5] + final_out[5,5], 
                                  final_out[6,5] + final_out[7,5] + final_out[8,5] + final_out[9,5], 
                                  final_out[12,5] + final_out[13,5], 
                                  final_out[10,5] + final_out[11,5] + final_out[14,5]),
                         myorder = c("1", "2", "3", "4", "5"), 
                         mycol = as.factor(c(1:5)))

final_out5$LTRE_scaled = as.numeric(final_out5$LTRE / sum(abs(final_out5$LTRE)))

# Create Sensitivity figure
sensitivity_mow <- ggplot(as.data.frame(final_out), 
                          aes(y = sensitivity, x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") + 
  ylim (0, 70) +
  ggtitle("(C) Mowing ambient * mowing future") + 
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_x_discrete(labels = mylabels_VR) +
  scale_fill_manual(values =  mypalette) + 
  labs(y = expression(paste("Sensitivity of ", lambda)))

# Create scaled LTRE figure for life cycle stages
LTRE_mow <- ggplot(as.data.frame(final_out5), aes(y = LTRE_scaled , x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ggtitle("(C) Mowing") + 
  ylim (-0.82, 0.82) +
  labs(y = expression(paste(Delta, lambda, " (ambient - future)"))) +
  theme_bw() + theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12)) +
  scale_fill_manual(values =  mypalette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = mylabels_LTRE)

## grazing

final_out <- data.frame()
pars_man <- c()

for(i in 1:14)
{
  final_out[i, "parameter"] <- colnames(par_fm)[i]  
  
  # mean pars ag fg:
  pars_man <- (par_ag + par_fg) / 2
  
  # "original" Lambda from combined pars
  final_out[i, "lambda_ori"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
  
  # calculate Lambda from manipulated combined pars
  pars_man[1,i] <- pars_man[1,i] + 0.001
  final_out[i, "lambda_fal"] <- Re(eigen(ker_mee(pars_man)$k_yx)$value[1])
}
final_out

# calculate sensitivity 
final_out[, "sensitivity"] <- as.numeric((final_out$lambda_fal - final_out$lambda_ori)/0.001)

# calculate LTRE values
final_out[, "LTRE"] <- as.numeric(final_out$sensitivity * gra_dif[1:14, "values"])

# establish order and clour for plots
final_out$myorder   <- as.factor(c(1, 2, 3, 4, 5, 8, 9, 6, 7, 13, 14, 10, 11, 12))
final_out$mycol       <- as.factor(c(1, 1, 2, 2, 2, 3, 3, 3, 3, 5, 5, 4, 4, 5))

# for ecological interpretation summerize vital rates into life cycle stages
final_out5 <- data.frame(pars = c("Survival", "Growth", "Reproduction", "Recruitment", "Sl survival"),
                         LTRE = c(final_out[1,5] + final_out[2,5], 
                                  final_out[3,5] + final_out[4,5] + final_out[5,5], 
                                  final_out[6,5] + final_out[7,5] + final_out[8,5] + final_out[9,5], 
                                  final_out[12,5] + final_out[13,5], 
                                  final_out[10,5] + final_out[11,5] + final_out[14,5]),
                         myorder = c("1", "2", "3", "4", "5"), 
                         mycol = as.factor(c(1:5)))

final_out5$LTRE_scaled = as.numeric(final_out5$LTRE / sum(abs(final_out5$LTRE)))

# Create Sensitivity figure
sensitivity_gra <- ggplot(as.data.frame(final_out), aes(y = sensitivity, x = myorder, fill = mycol)) + geom_bar(stat = "identity") + 
  ylim (0, 70) +
  ggtitle("(D) Grazing ambient * grazing future") + #ylab("Sensitivity of  \u03BB") +
  theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_x_discrete(labels = mylabels_VR) +
  scale_fill_manual(values =  mypalette) 

# Create scaled LTRE figure for life cycle stages
LTRE_gra <- ggplot(as.data.frame(final_out5), aes(y = LTRE_scaled , x = myorder, fill = mycol)) + 
  geom_bar(stat = "identity") +
  ggtitle("(D) Grazing") + 
  ylim (-0.82, 0.82) +
  theme_bw() + theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_fill_manual(values =  mypalette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = mylabels_LTRE)

Sensitivity_all <- grid.arrange(sensitivity_amb, sensitivity_fut, sensitivity_mow, sensitivity_gra,
                                ncol = 2, nrow = 2, 
                                layout_matrix = rbind(c(1, 2), c(3, 4)), 
                                widths = c(5.2, 5.2), heights = c(5, 5.8))
Sensitivity_all
ggsave(filename = "C:\\Users/ma22buky/Documents/Julia_Paper/Fig_Sensitivity.jpeg", 
       plot = Sensitivity_all, 
       device = "jpeg", 
       dpi = res, 
       width = 20,
       height =15, 
       units = "cm")

LTRE_all <- grid.arrange(LTRE_amb, LTRE_fut, LTRE_mow, LTRE_gra,
                         ncol = 2, nrow = 2, 
                         layout_matrix = rbind(c(1, 2), c(3, 4)), 
                         widths = c(5.5, 5), heights = c(4.8, 5.9))

LTRE_all
ggsave(filename = "C:\\Users/ma22buky/Documents/Julia_Paper/Fig_LTRE.tiff", 
       plot = LTRE_all, 
      device = "tiff", 
       dpi = res, 
       width = 20,
       height =15, 
       units = "cm")

## sensitivity just treatments
sens_amb_gra <- data.frame()
sens_fut_gra <- data.frame()
sens_amb_mow <- data.frame()
sens_fut_mow <- data.frame()



for(i in 1:14)
{
  ## ambient grazing
  sens_amb_gra[i, "ori_lambda"] <- Re(eigen(ker_mee(par_ag)$k_yx)$value[1])
  par_ag_new <- par_ag
  par_ag_new[1,i] <- par_ag[1, i] + 0.001
  sens_amb_gra[i, "new_lambda"] <- Re(eigen(ker_mee(par_ag_new)$k_yx)$value[1])
  sens_amb_gra[i, "changed_par"] <- names(par_ag_new)[i]
  
  ## ambient mowing
  sens_amb_mow[i, "ori_lambda"] <- Re(eigen(ker_mee(par_am)$k_yx)$value[1])
  par_am_new <- par_am
  par_am_new[1,i] <- par_am[1, i] + 0.001
  sens_amb_mow[i, "new_lambda"] <- Re(eigen(ker_mee(par_am_new)$k_yx)$value[1])
  sens_amb_mow[i, "changed_par"] <- names(par_am_new)[i]
  
  ## future grazing
  sens_fut_gra[i, "ori_lambda"] <- Re(eigen(ker_mee(par_fg)$k_yx)$value[1])
  par_fg_new <- par_fg
  par_fg_new[1,i] <- par_fg[1, i] + 0.001
  sens_fut_gra[i, "new_lambda"] <- Re(eigen(ker_mee(par_fg_new)$k_yx)$value[1])
  sens_fut_gra[i, "changed_par"] <- names(par_fg_new)[i]
  
  ## future mowing
  sens_fut_mow[i, "ori_lambda"] <- Re(eigen(ker_mee(par_fm)$k_yx)$value[1])
  par_fm_new <- par_fm
  par_fm_new[1,i] <- par_fm[1, i] + 0.001
  sens_fut_mow[i, "new_lambda"] <- Re(eigen(ker_mee(par_fm_new)$k_yx)$value[1])
  sens_fut_mow[i, "changed_par"] <- names(par_fm_new)[i]
  
}

sens_amb_gra_fin <- data.frame(sensitivity_amb = (sens_amb_gra$new_lambda - sens_amb_gra$ori_lambda)/0.001,
                               par = sens_amb_gra$changed_par)
sens_amb_mow_fin <- data.frame(sensitivity_amb = (sens_amb_mow$new_lambda - sens_amb_mow$ori_lambda)/0.001,
                               par = sens_amb_mow$changed_par)
sens_fut_gra_fin <- data.frame(sensitivity_fut = (sens_fut_gra$new_lambda - sens_fut_gra$ori_lambda)/0.001,
                               par = sens_fut_gra$changed_par)
sens_fut_mow_fin <- data.frame(sensitivity_fut = (sens_fut_mow$new_lambda - sens_fut_mow$ori_lambda)/0.001,
                               par = sens_fut_mow$changed_par)

## change order to desired order
sens_amb_gra_fin$par <- factor(sens_amb_gra_fin$par, 
                               levels = c("surv_b0", "surv_b1", "grow_b0", "grow_b1", 
                                          "grow_sd", "seednum_b0", "seednum_b1", 
                                          "flowprob_b0", "flowprob_b1", "novSL_per_fl",
                                          "aprilSL_per_fl", "SL_surv", "new_size_mean",
                                          "new_size_sd"))
sens_amb_gra_fin <- sens_amb_gra_fin[order(sens_amb_gra_fin$par),]
sens_amb_gra_fin$par_names <- factor(mylabels_VR, 
                                     levels = c("S Int", "S Slope", "G Int", "G Slope", 
                                                "G SD", "P Int", "P Slope", "F Int", "F Slope", 
                                                " \u03B8f", "\u03B8s", "B", "\u03B7", 
                                                "\u03B7 SD"))

sens_amb_mow_fin$par <- factor(sens_amb_mow_fin$par, 
                               levels = c("surv_b0", "surv_b1", "grow_b0", "grow_b1", 
                                          "grow_sd", "seednum_b0", "seednum_b1", 
                                          "flowprob_b0", "flowprob_b1", "novSL_per_fl",
                                          "aprilSL_per_fl", "SL_surv", "new_size_mean",
                                          "new_size_sd"))
sens_amb_mow_fin <- sens_amb_mow_fin[order(sens_amb_mow_fin$par),]
sens_amb_mow_fin$par_names <- factor(mylabels_VR, 
                                     levels = c("S Int", "S Slope", "G Int", "G Slope", 
                                                "G SD", "P Int", "P Slope", "F Int", "F Slope", 
                                                " \u03B8f", "\u03B8s", "B", "\u03B7", 
                                                "\u03B7 SD"))

sens_fut_gra_fin$par <- factor(sens_fut_gra_fin$par, 
                               levels = c("surv_b0", "surv_b1", "grow_b0", "grow_b1", 
                                          "grow_sd", "seednum_b0", "seednum_b1", 
                                          "flowprob_b0", "flowprob_b1", "novSL_per_fl",
                                          "aprilSL_per_fl", "SL_surv", "new_size_mean",
                                          "new_size_sd"))
sens_fut_gra_fin <- sens_fut_gra_fin[order(sens_fut_gra_fin$par),]
sens_fut_gra_fin$par_names <- factor(mylabels_VR, 
                                     levels = c("S Int", "S Slope", "G Int", "G Slope", 
                                                "G SD", "P Int", "P Slope", "F Int", "F Slope", 
                                                " \u03B8f", "\u03B8s", "B", "\u03B7", 
                                                "\u03B7 SD"))

sens_fut_mow_fin$par <- factor(sens_fut_mow_fin$par, 
                               levels = c("surv_b0", "surv_b1", "grow_b0", "grow_b1", 
                                          "grow_sd", "seednum_b0", "seednum_b1", 
                                          "flowprob_b0", "flowprob_b1", "novSL_per_fl",
                                          "aprilSL_per_fl", "SL_surv", "new_size_mean",
                                          "new_size_sd"))
sens_fut_mow_fin <- sens_fut_mow_fin[order(sens_fut_mow_fin$par),]
sens_fut_mow_fin$par_names <- factor(mylabels_VR, 
                                     levels = c("S Int", "S Slope", "G Int", "G Slope", 
                                                "G SD", "P Int", "P Slope", "F Int", "F Slope", 
                                                " \u03B8f", "\u03B8s", "B", "\u03B7", 
                                                "\u03B7 SD"))

amb_mow_pl <- ggplot(sens_amb_mow_fin, aes(x = par_names, y = sensitivity_amb)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(title = "(B) Ambient mowing", x = "", y = "") +
  ylim(0, 60) 

amb_gra_pl <- ggplot(sens_amb_gra_fin, aes(x = par_names, y = sensitivity_amb)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(title = "(A) Ambient grazing", y = expression(paste("Sensitivity of ", 
                                                           lambda)), x = "") + ylim(0, 60)

fut_mow_pl <- ggplot(sens_fut_mow_fin, aes(x = par_names, y = sensitivity_fut)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "(D) Future mowing", y = "", x = "") + ylim(0, 60)

fut_gra_pl <- ggplot(sens_fut_gra_fin, aes(x = par_names, y = sensitivity_fut)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "(C) Future grazing", y = expression(paste("Sensitivity of ", 
                                                          lambda)), x = "") + ylim(0, 60)

sensitivity_all <- grid.arrange(amb_gra_pl, amb_mow_pl, fut_gra_pl, fut_mow_pl)
ggsave(filename = "C:\\Users/ma22buky/Documents/Julia_Paper/Fig_sens.tiff", 
       plot = sensitivity_all, 
       device = "tiff", 
       dpi = res, 
       width = 20,
       height =15, 
       units = "cm")

## sensitivity + elasticity
heat_colors <- sort(heat.colors(10), decreasing = T)

## ambient mowing
par_am$n <-50

k_yx <- ker_mee(par_am)$k_yx
lam   <- Re(eigen(k_yx)$values[1])
ev    <- eigen(k_yx) 

W <- ev$vectors 
w <- abs(Re(W[, 1])) #w is the right eigenvector, this also describes the reproductive value of each stage class.  Reproductive value depends on reproduction, survival and timing (i.e., an individual has to live to reach a reproductive stage in order to have reproductive value.  Low values could indicate a small proportion of individuals in a stage class survive to reproductive age).  
V <- try(Conj(solve(W)), silent = TRUE)
v <- abs(Re(V[1, ])) #v is the left eigenvector, this also describes he proportion of individuals in each stage class at stable stage distribution
s <- v %o% w #sensitivity of lambda to changes in each matrix element (aij) is proportional to the product of the ith element of the left eigenvector and the jth element of the right eigenvector 
rotate <- function(x) t(apply(x, 2, rev))
e <- s * k_yx/lam #elasticity is the proportional sensitivity of lambda to changes in each matrix element 
e <- rotate(e)

image(e)
image(rotate(s))

elast_amb_mow <- reshape::melt(e) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(D) Elasticity - ambient mowing", y = "", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25) 
  
sens_amb_mow <- reshape::melt(rotate(s)) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(C) Sensitivity - ambient mowing", y = "Size (t + 1)", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25) 

## ambient grazing
par_ag$n <- 50
k_yx <- ker_mee(par_ag)$k_yx
lam   <- Re(eigen(k_yx)$values[1])
ev    <- eigen(k_yx) 

W <- ev$vectors 
w <- abs(Re(W[, 1])) #w is the right eigenvector, this also describes the reproductive value of each stage class.  Reproductive value depends on reproduction, survival and timing (i.e., an individual has to live to reach a reproductive stage in order to have reproductive value.  Low values could indicate a small proportion of individuals in a stage class survive to reproductive age).  
V <- try(Conj(solve(W)), silent = TRUE)
v <- abs(Re(V[1, ])) #v is the left eigenvector, this also describes he proportion of individuals in each stage class at stable stage distribution
s <- v %o% w #sensitivity of lambda to changes in each matrix element (aij) is proportional to the product of the ith element of the left eigenvector and the jth element of the right eigenvector 
rotate <- function(x) t(apply(x, 2, rev))
e <- s * k_yx/lam #elasticity is the proportional sensitivity of lambda to changes in each matrix element 
e <- rotate(e)

image(e)
image(rotate(s))

elast_amb_gra <- reshape::melt(e) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(B) Elasticity - ambient grazing", y = "", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25) 

sens_amb_gra <- reshape::melt(rotate(s)) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(A) Sensitivty - ambient grazing", y = "Size (t + 1)", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25) 

##future mowing
par_fm$n <- 50
k_yx <- ker_mee(par_fm)$k_yx
lam   <- Re(eigen(k_yx)$values[1])
ev    <- eigen(k_yx) 

W <- ev$vectors 
w <- abs(Re(W[, 1])) #w is the right eigenvector, this also describes the reproductive value of each stage class.  Reproductive value depends on reproduction, survival and timing (i.e., an individual has to live to reach a reproductive stage in order to have reproductive value.  Low values could indicate a small proportion of individuals in a stage class survive to reproductive age).  
V <- try(Conj(solve(W)), silent = TRUE)
v <- abs(Re(V[1, ])) #v is the left eigenvector, this also describes he proportion of individuals in each stage class at stable stage distribution
s <- v %o% w #sensitivity of lambda to changes in each matrix element (aij) is proportional to the product of the ith element of the left eigenvector and the jth element of the right eigenvector 
rotate <- function(x) t(apply(x, 2, rev))
e <- s * k_yx/lam #elasticity is the proportional sensitivity of lambda to changes in each matrix element 
e <- rotate(e)

image(e)
image(rotate(s))

elast_fut_mow <- reshape::melt(e) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(H) Elasticity - future mowing", y = "", 
       x = "Size (t)") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25)  

sens_fut_mow <- reshape::melt(rotate(s)) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(G) Sensitivity - future mowing", y = "Size (t + 1)", 
       x = "Size (t)") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25) 

## future grazing
par_fg$n <- 50
k_yx <- ker_mee(par_fg)$k_yx
lam   <- Re(eigen(k_yx)$values[1])
ev    <- eigen(k_yx) 

W <- ev$vectors 
w <- abs(Re(W[, 1])) #w is the right eigenvector, this also describes the reproductive value of each stage class.  Reproductive value depends on reproduction, survival and timing (i.e., an individual has to live to reach a reproductive stage in order to have reproductive value.  Low values could indicate a small proportion of individuals in a stage class survive to reproductive age).  
V <- try(Conj(solve(W)), silent = TRUE)
v <- abs(Re(V[1, ])) #v is the left eigenvector, this also describes he proportion of individuals in each stage class at stable stage distribution
s <- v %o% w #sensitivity of lambda to changes in each matrix element (aij) is proportional to the product of the ith element of the left eigenvector and the jth element of the right eigenvector 
rotate <- function(x) t(apply(x, 2, rev))

e <- s * k_yx/lam #elasticity is the proportional sensitivity of lambda to changes in each matrix element 
e <- rotate(e)

image(e)
image(rotate(s))

elast_fut_gra <- reshape::melt(e) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(F) Elasticity - future grazing", y = "", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25)

sens_fut_gra <- reshape::melt(rotate(s)) %>% ggplot(aes(X1, X2, fill = value)) + geom_tile() + 
  theme_bw() + scale_fill_gradientn(colors = heat_colors) +
  labs(title = "(E) Sensitivty - future grazing", y = "Size (t + 1)", 
       x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_blank(), aspect.ratio = 0.25)

heat_maps <- ggpubr::ggarrange(sens_amb_gra, elast_amb_gra, sens_amb_mow, elast_amb_mow, sens_fut_gra,
             elast_fut_gra, sens_fut_mow, elast_fut_mow, ncol = 2,nrow = 4, align = "hv")

heat_maps
ggsave(filename = "C:\\Users/ma22buky/Documents/Julia_Paper/heat.tiff", 
       plot = heat_maps, 
       device = "tiff", 
       dpi = res, 
       width = 23.11,
       height = 15, 
       units = "cm")
