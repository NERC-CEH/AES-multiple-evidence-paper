

#### Various plotting functions

#### fit distribution
# function to determine which distribution is right for the data
# negative binomial or poisson
# data - data containing the variable of interest
# var - column of interest

# outputs
# mean and variance of column
# histgram of data with the two distributions overlayed
fit_distrib <- function(dat, var) {
  
  require(fitdistrplus)
  
  print(data.frame('mean' = mean(dat[, var]),
                   'variance' = var(dat[, var])))
  
  hist(dat[, var], prob = TRUE, breaks = 20, main = var)
  
  fit <- fitdist(dat[, var], "nbinom")
  fitD <- dnbinom(0:max(dat[, var]), size = fit$estimate[1], mu = fit$estimate[2])
  lines(fitD, col = "blue")
  
  fit <- fitdist(dat[, var], "pois")
  fitD <- dpois(0:max(dat[, var]), lambda = fit$estimate[1])
  lines(fitD, col = "red")
  
  legend("topright",legend = c('neg_binom', 'poisson'),
         col = c('blue', 'red'), lty=1)
  
  par(mfrow = c(2,2))
  #scatterplots
  plot(dat[, var] ~ dat[, 'AES1KM'], ylab = var)
  plot(dat[, var] ~ dat[, 'AES3KM'], ylab = var)
  
  #boxplot with year
  
  boxplot(dat[, var] ~ dat[, 'SURVEY_YEAR'], ylab = var)
  if(any(grepl("ROUND_NUMBER", names(dat)))){
    boxplot(dat[, var] ~ dat[, 'ROUND_NUMBER'], ylab = var)
  }
  if(any(grepl("NUMBER_TRAPS", names(dat)))){
    boxplot(dat[, var] ~ dat[, 'NUMBER_TRAPS'], ylab = var)
  }
  par(mfrow = c(1,1))
  
}



#### get predictictions and confidence intervals
# takes a model and new dataset to predict from 
# and returns the line of best fit 
# and cofidence intervals

# mod - the model, works with class = lme or glmer.nb
# pred data - the 'newx' data to predict from
# col_id - the name to call the line of best fit column
# cmult - approximation of normal distribution value
# re.form - value for random effects in predict function, default = no REs
# level - which level to predict for
# re_name - the name of the random effects column

pred_conf <- function(mod, pred_data, col_id, cmult = 1.96, re.form = NA, level = 0, re_name = 'NCA') {
  
  # predict from the model
  if(class(mod)=='lme') {
    pred_data[, col_id] <- predict(mod, pred_data, level = level)
  } else {
    pred_data[, col_id] <- predict(mod, pred_data, re.form=re.form)
  }
  
  # get the confidence intervals
  mm <- model.matrix(terms(mod), pred_data)
  ## or newdat$distance <- mm %*% fixef(fm1)
  pvar1 <- diag(mm %*% tcrossprod(vcov(mod),mm))
  
  if(class(mod)=='lme') {
    tvar1 <- pvar1+as.numeric(VarCorr(mod)[rownames(VarCorr(mod))=="(Intercept)"][1])  ## must be adapted for more complex models
  } else {
    tvar1 <- pvar1+as.numeric(VarCorr(mod)[re_name][1])  ## must be adapted for more complex models
  }
  pred_data_out <- data.frame(
    pred_data
    , plo = pred_data[, col_id]-cmult*sqrt(pvar1)
    , phi = pred_data[, col_id]+cmult*sqrt(pvar1)
    , tlo = pred_data[, col_id]-cmult*sqrt(tvar1)
    , thi = pred_data[, col_id]+cmult*sqrt(tvar1)
  )
  
  return(pred_data_out)
  
}




## function for rescaling data and generating + saving plots
# predict_data - dataset containing the predicted data
# xvars - the xvariable for plotting (same as the variable in the original data)
# yvar - the name of the column containing the predicted data
# original data - the original dataset (not the predicted one) used to unscale the data
# unscale - whether or not to unscale the data
# log_trans - whether or not to exponential transform the fitted lines - used after log transforming data for analysis
# x_var_names - name for the x axis and what to call the saved file
# saving name - name to give the saved file.

plot_results <- function(predict_data, xvars, yvar, original_data, scaled_data = NULL, unscale = TRUE, log_trans = FALSE,
                         x_var_names, ylim = c(1,1200), boxplot_x, plot.name = plot_name, saving_name = analysis_name, 
                         write = TRUE, lowland_only = FALSE, plot_dir) {
  
  require(tidyverse)
  
  if((lowland_only & !grepl('lowland', plot_dir))|
     (!lowland_only & grepl('lowland', plot_dir))) stop('Running lowland only analysis and not saving in lowland folder, or vice versa. Set lowland_only to TRUE or edit directory')
  
  if((lowland_only & is.null(scaled_data))|
     (!lowland_only & !is.null(scaled_data))) stop("Running lowland only analysis without providing a scaled dataset to unscale variables, or vice versa.\nIf running lowland only, provide a scaled dataset (data from all NCA analysis).")
  
  predicted_data <- predict_data
  original_data$boxplot_y <- original_data[, yvar]
  original_data$boxplot_x <- original_data[, boxplot_x]
  
  #set y limits
  if(!any(grepl('DIV', colnames(predict_data))) | log_trans) {
    ylim <- c(floor(min(c(exp(predicted_data$plo), original_data[, yvar]))), ceiling(max(c(exp(predicted_data$phi),original_data[, yvar]))))
  } else if(any(grepl('DIV', colnames(predict_data)))) {
    ylim <- c(floor(min(c(predicted_data$plo, original_data[, yvar]))), ceiling(max(c(predicted_data$phi, original_data[, yvar]))))
  }
  
  # if there's no scaled data, then scale using original data
  if(is.null(scaled_data)) scaled_data <- original_data
  
  if(unscale) {
    
    predicted_data[ , xvars] <- predicted_data[ , xvars] * attr(scaled_data[ , xvars], 'scaled:scale') + attr(scaled_data[ , xvars], 'scaled:center')
    unscaled_x <- c(original_data[ , xvars] * attr(scaled_data[ , xvars], 'scaled:scale') + attr(scaled_data[ , xvars], 'scaled:center'))
    
  
    ## plotting
    if(write) {
      
      png(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", x_var_names, " AES scatterplot.png"), width = 8.4, height = 8.4, pointsize = 8, units = "cm", res = 300)
      plot(original_data[, yvar] ~ unscaled_x, 
           xlab = paste0(strsplit(x_var_names, " ")[[1]][1], " AES gradient score"), 
           ylab = ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name), 
           pch = 20, col = "grey20",
           ylim = ylim, bty = 'l')
      
      if(!any(grepl('DIV', colnames(predict_data))) | log_trans) { # if not diversity model then exp transform
        
        lines(exp(predicted_data[, yvar]) ~ predicted_data[, xvars])
        lines(exp(predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
        lines(exp(predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
        
      } else if(any(grepl('DIV', colnames(predict_data)))) { # if there is diversity then don't exp transform
        
        lines((predicted_data[, yvar]) ~ predicted_data[, xvars])
        lines((predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
        lines((predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
        
      }
      
      
      dev.off()
      
      if(!lowland_only){
        
        boxp <- ggplot(original_data, aes(x = boxplot_x, y = boxplot_y, fill = SURVEY_YEAR)) +
          ylab(ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name)) +
          # stat_boxplot(geom = "errorbar")  +
          geom_boxplot() + 
          scale_fill_grey(start=0.5, end = 1) +
          # lims(y = ylim) + 
          theme_classic() +
          theme(legend.title = element_blank(),
                axis.title.x = element_blank())
        
        ggsave(boxp, filename = paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " per ", boxplot_x, " and year.png"),
               width = 18, height = 10, units = "cm")
      }
      
    } 
    
    plot(original_data[, yvar] ~ unscaled_x, 
         xlab = paste0(strsplit(x_var_names, " ")[[1]][1], " AES gradient score"), 
         ylab = ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name), 
         pch = 20, col = "grey20",
         ylim = ylim, bty = 'l')
    
    if(!any(grepl('DIV', colnames(predict_data))) | log_trans) { # if not diversity model then exp transform
      
      lines(exp(predicted_data[, yvar]) ~ predicted_data[, xvars])
      lines(exp(predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
      lines(exp(predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
      
    } else if(any(grepl('DIV', colnames(predict_data)))) { # if there is diversity then don't exp transform
      
      lines((predicted_data[, yvar]) ~ predicted_data[, xvars])
      lines((predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
      lines((predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
      
    }
    
    
    if(!lowland_only){
      
      ggplot(original_data, aes(x = boxplot_x, y = boxplot_y, fill = SURVEY_YEAR)) +
        ylab(ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name)) +
        # stat_boxplot(geom = "errorbar")  +
        geom_boxplot() + 
        scale_fill_grey(start=0.5, end = 1) +
        # lims(y = ylim) + 
        theme_classic() +
        theme(legend.title = element_blank(),
              axis.title.x = element_blank())
      
    }
  }
  
}



###Interaction plot###

# predict_data - dataset containing the predicted data
# xvars - the xvariable for plotting (same as the variable in the original data)
# yvar - the name of the column containing the predicted data
# original data - the original dataset (not the predicted one) used to unscale the data
# unscale - whether or not to unscale the data
# x_var_names - name for the x axis and what to call the saved file
# saving name - name to give the saved file.

interaction_plot <- function(model, yvar, original_data, scaled_data = NULL, unscale = TRUE,
                             ylim = c(0,6), plot.name = plot_name, saving_name = analysis_name, 
                             write = TRUE, lowland_only = FALSE, plot_dir){
  
  require(interactions)
  require(ggplot2)
  
  if((lowland_only & !grepl('lowland', plot_dir))|
     (!lowland_only & grepl('lowland', plot_dir))) stop('Running lowland only analysis and not saving in lowland folder, or vice versa. Set lowland_only to TRUE or edit directory')
  
  if((lowland_only & is.null(scaled_data))|
     (!lowland_only & !is.null(scaled_data))) stop("Running lowland only analysis without providing a scaled dataset to unscale variables, or vice versa.\nIf running lowland only, provide a scaled dataset (data from all NCA analysis).")
  
  if(class(model)[1] != "glmerMod"){ # if it's already in lme4 - leave it, will work with interact_plot function
    if(length(grep("lme.formula", model$call)) >0){
      #fit same model in lme4 to allow interaction plot
      model <- lmer(formula(paste(paste(model$call$fixed)[2], "~", model$call$fixed[3], "+",paste0("(", model$call$random,")")[2])), data = original_data)
    }
  }
  
  # if there's no scaled data, then scale using original data
  if(is.null(scaled_data)) scaled_data <- original_data
  
  
  xlab_local <- c(0,10000,20000,30000)
  
  breaks_local <- (c(0,10000,20000,30000) - attr(scaled_data[ , "AES1KM"], 'scaled:center')) / attr(scaled_data[ , "AES1KM"], 'scaled:scale')
  
  xlab_landscape <- c(0,2000,4000,6000,8000,10000,12000,14000)
  
  breaks_landscape <- (c(0,2000,4000,6000,8000,10000,12000,14000) - attr(scaled_data[ , "AES3KM"], 'scaled:center')) / attr(scaled_data[ , "AES3KM"], 'scaled:scale')
  
  
  #png(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", " interaction plot 1.png"), height = 12, width = 16, units = "cm", res = 300)
  int1 <- interact_plot(model, pred = "AES1KM", modx = "AES3KM", plot.points = TRUE, modx.values = c(-0.946, -0.311, 0.994), modx.labels = c("Low", "Medium", "High"), cex.lab = 4, cex.axis = 4, x.label = "Local AES gradient score", y.label = plot.name, legend.main = "Landscape AES \ngradient score", colors = "Greys", interval = TRUE, int.width = 0.95)
  int1 + scale_y_continuous(limits = ylim) + scale_x_continuous(breaks = breaks_local, labels = xlab_local)
  if(write) ggsave(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", "interaction plot 1.png"))
  
  
  #png(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", " interaction plot 2.png"), height = 12, width = 16, units = "cm", res = 300)
  int2 <- interact_plot(model, pred = "AES3KM", modx = "AES1KM", plot.points = TRUE, modx.values = c(-0.831, -0.221, 0.84), modx.labels = c("Low", "Medium", "High"), cex.lab = 4, cex.axis = 4, x.label = "Landscape AES gradient score", y.label = plot.name, legend.main = "Local AES \ngradient score", colors = "Greys", interval = TRUE, int.width = 0.95)
  int2 + scale_y_continuous(limits = ylim) + scale_x_continuous(breaks = breaks_landscape, labels = xlab_landscape)
  if(write) ggsave(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", "interaction plot 2.png"))
  
}


## write model outputs to file in a nice format
# model - model to save
# saving_name - name of the analysis, as in other functions
# NCA - lowland or upland? one of 'Lowland only' or 'All NCA'
# outlier - is this an outlier model? One of 'Outlier removed' or 'All data'
# write - whether to write, TRUE / FALSE
# save_dir - where to save. Setting it to the main folder (not split into lowland only)
model_table <- function(model, saving_name = analysis_name, NCA, outlier,
                        write = TRUE, 
                        save_dir = paste0(dir$directories$outpath, analysis_group, "/", analysis_name, "/")) {
  
  require(broom.mixed)
  require(tidyverse)
  
  if(class(model) == 'lme'){
    
    require(lmerTest)
    
    #fit same model in lme4 to allow interaction plot
    model <- lmer(formula(paste(paste(model$call$fixed)[2], "~", model$call$fixed[3], "+",paste0("(", model$call$random,")")[2])), data = model$data)
    
  }
  
  mod_tab <- data.frame(Analysis = saving_name,
                        NCA = NCA,
                        Outlier = outlier,
                        tidy(model)) %>%
    mutate_if(is.numeric, round, 3) %>% 
    replace(is.na(.), '/') %>% 
    dplyr::rename(Effect = effect, 
                  Group = group, 
                  Term = term,
                  Estimate = estimate,
                  'Std error' = std.error,
                  Statistic = statistic, 
                  'P value' = p.value)
  
  if(write) write.csv(mod_tab, file = paste0(save_dir, "/", NCA, " ", saving_name, "_", outlier, ".csv"))
  
  return(mod_tab)
  
}


## plot covariate effects

plot_covariates <- function(predict_data, xvars, yvar, original_data, scaled_data = NULL, unscale = TRUE, log_trans = FALSE,
                            x_var_names, ylim = c(1,1200), plot.name = plot_name, saving_name = analysis_name, 
                            write = TRUE, lowland_only = FALSE, plot_dir) {
  
  require(tidyverse)
  
  if((lowland_only & !grepl('lowland', plot_dir))|
     (!lowland_only & grepl('lowland', plot_dir))) stop('Running lowland only analysis and not saving in lowland folder, or vice versa. Set lowland_only to TRUE or edit directory')
  
  if((lowland_only & is.null(scaled_data))|
     (!lowland_only & !is.null(scaled_data))) warning("Running lowland only analysis without providing a scaled dataset to unscale variables, or vice versa.\nIf running lowland only, provide a scaled dataset (data from all NCA analysis).")
  
  predicted_data <- predict_data
  
  
  #set y limits
  if(!any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"])) | log_trans) {
    ylim <- c(floor(min(c(exp(predicted_data$plo), original_data[, yvar]))), ceiling(max(c(exp(predicted_data$phi),original_data[, yvar]))))
  } else if(any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"]))) {
    ylim <- c(floor(min(c(predicted_data$plo, original_data[, yvar]))), ceiling(max(c(predicted_data$phi, original_data[, yvar]))))
  }
  
  # if there's no scaled data, then scale using original data
  if(is.null(scaled_data)) scaled_data <- original_data
  
  if(unscale) {
    
    predicted_data[ , xvars] <- predicted_data[ , xvars] * attr(scaled_data[ , xvars], 'scaled:scale') + attr(scaled_data[ , xvars], 'scaled:center')
    unscaled_x <- c(original_data[ , xvars] * attr(scaled_data[ , xvars], 'scaled:scale') + attr(scaled_data[ , xvars], 'scaled:center'))
    
  } else {unscaled_x <- original_data[,xvars]}
  
  ## plotting
  if(write) {
    
    png(paste0(plot_dir, "/", ifelse(lowland_only, 'Lowland only ', ''), saving_name, " ", x_var_names, " plot.png"), width = 8.4, height = 8.4, pointsize = 8, units = "cm", res = 300)
    plot(original_data[, yvar] ~ unscaled_x, 
         xlab = x_var_names, 
         ylab = ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name), 
         pch = 20, col = "grey20",
         ylim = ylim, bty = 'l', cex.lab = 1.3)
    
    if(!any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"])) | log_trans) { # if not diversity model then exp transform
      
      lines(exp(predicted_data[, yvar]) ~ predicted_data[, xvars])
      lines(exp(predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
      lines(exp(predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
      
    } else if(any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"]))) { # if there is diversity then don't exp transform
      
      lines((predicted_data[, yvar]) ~ predicted_data[, xvars])
      lines((predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
      lines((predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
      
    }
    
    
    dev.off()
    
    
  } 
  
  plot(original_data[, yvar] ~ unscaled_x, 
       xlab =x_var_names, 
       ylab = ifelse(lowland_only, paste0(plot.name, ' - lowland only'), plot.name), 
       pch = 20, col = "grey20",
       ylim = ylim, bty = 'l')
  
  if(!any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"])) | log_trans) { # if not diversity model then exp transform
    
    lines(exp(predicted_data[, yvar]) ~ predicted_data[, xvars])
    lines(exp(predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
    lines(exp(predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
    
  } else if(any(grepl('DIV', colnames(predict_data)[colnames(predict_data) != "BOT_DIV"]))) { # if there is diversity then don't exp transform
    
    lines((predicted_data[, yvar]) ~ predicted_data[, xvars])
    lines((predicted_data$plo) ~ predicted_data[, xvars], lty = 2)
    lines((predicted_data$phi) ~ predicted_data[, xvars], lty = 2)
    
  }
  
  
  
  
}

