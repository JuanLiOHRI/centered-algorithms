library(ggpubr)

plot.calibration <- function(
  df,
  color = "black",
  titlestr = "Calibration Plot",
  plot.o = NULL
) {
  # -------------------------------------------------------------------
  # recreate the calibration plot from `rms::val.prob` for exploration
  # -------------------------------------------------------------------
  # https://stats.stackexchange.com/questions/563867/walk-through-rmsval-prob
  # View(rms::val.prob)
  # 2026-01-15 by Juan Li

  # Logistic calibration
  fit <- glm(y ~ logit, data = df, family = binomial)
  df$lc <- predict(fit, newdata = df, type = "response")
  # Nonparametric (loess)
  df$Sm.x <- lowess(df$p, df$y, iter = 0)$x
  df$Sm.y <- lowess(df$p, df$y, iter = 0)$y

  # Plot
  df_Diag <- data.frame(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1))
  if (!is.null(plot.o)) {
    p1 <- plot.o +
      geom_line(data = df, aes(p, lc), color = color) +
      geom_line(
        data = df,
        aes(Sm.x, Sm.y),
        linetype = "dashed",
        color = color
      ) +
      theme_bw() +
      labs(title = titlestr)
  } else {
    p1 <- ggplot(df, aes(p, lc)) +
      geom_line(color = color) +
      geom_line(aes(Sm.x, Sm.y), linetype = "dashed", color = color) +
      geom_line(
        data = df_Diag,
        aes(x, y),
        color = color,
        linewidth = 2,
        alpha = 0.1
      ) +
      theme_bw() +
      labs(
        x = "Predicted Probability",
        y = "Observed Probability"
      )
  }

  return(list(df = df, plot = p1))
}

# ==================

calibration <- function(prediction, outcome, package = "Steyerberg", ...)
{
  # input:
  #   outcome:       actual outcome
  #   prediction:    model predicted values, should be the same length as outcome.
  #   group:         a vector of the same length as outcome, subgroup categories for strong calibration
  #   package:       can be "Steyerberg" or "Harrell"
  #   ...:      optional arguments for insider functions
  
  # output:          a named vector
  
  if (length(outcome) != length(prediction)) 
    stop("lengths of outcome and prediction do not agree")
  
  if (length(unique(outcome)) == 2)
  {
    # Logistic regression: based on rms::val.prob
    if (min(prediction, na.rm = TRUE) < 0 | max(prediction, na.rm = TRUE) > 1) 
      stop("Logistic regression: predicted probability should be between 0 and 1.")
    
    # remove probability 0 or 1
    outcome    <- outcome[prediction > 0 & prediction < 1]
    prediction <- prediction[prediction > 0 & prediction < 1]
    if (package == "Steyerberg")
    {
      stats <- CalibrationCurves::val.prob.ci.2(prediction, outcome, ...)
    } else
    {
      stats <- rms::val.prob(prediction, outcome, ...)
    }
  } else
  {
    stats <- calibrationCon(prediction, outcome, ...)
  }
  
  return(stats)
}

# ===========================

calibrationCon <- function(prediction, outcome, group = NULL, group_name = NULL, g = 1, violin = FALSE, marginPlt = FALSE, alpha = 0.5, nrow = 1,
labelPos.x = NULL, labelPos.y = NULL) {
  # prediction:  vector of the predicted outcome
  # outcome:     vector of the observed outcome
  # group:       vector of group labels if want to assess calibration within subgroups
  # group_name:  a string vector, levels of group
  # g:           number of bins if want to bin the data instead of showing all data points
  # violin:      set to TRUE if want to using violin plot instead of boxplot when binning the data
  # marginPlt:   if want to show distributions at the margin
  # alpha:       alpha value of data points
  # nrow:        number of rows in the facet
  # labelPos.x:  x position of the label 
  # labelPos.y:  y position of the label 
  
  if (!is.null(group))
  {
    if (length(outcome) != length(group)) 
      stop("lengths of outcome and group do not agree")
    group <- as.character(group)
  }
  
  df <- data.frame(prediction = prediction, outcome = outcome)
  if (!is.null(group)) {
    if (!is.null(group_name)) {
      df$group_fct <- group
      df$group_fct <- factor(df$group_fct, levels = group_name)
    } else {
      df$group_fct <- factor(group)
      group_name <- levels(df$group_fct)
    }
  } else {
    df <- df %>% mutate(group_fct = " ")
    df$group_fct <- factor(df$group_fct)
    group_name <- levels(df$group_fct)
  }
  df <- df %>% mutate(er = abs(outcome - prediction))
  
  min <- min(min(outcome, na.rm = TRUE), min(prediction, na.rm = TRUE))
  max <- max(max(outcome, na.rm = TRUE), max(prediction, na.rm = TRUE))

  if (is.null(labelPos.x)) labelPos.x <- min*1.1
  if (is.null(labelPos.y)) labelPos.y <- max*0.9
  
  # ------- stats -----------
  # whole data
  lr <- lm(outcome ~ prediction, data = df)
  eavg <- mean(df$er, na.rm = TRUE)
  emax <- max(df$er, na.rm = TRUE)
  e90 <- unname(quantile(df$er, 0.9, na.rm = TRUE))
  
  stats <- data.frame(group     = " ",
                      Intercept = unname(lr$coefficients[1]),
                      Slope     = unname(lr$coefficients[2]),
                      Emax      = emax,
                      E90       = e90,
                      Eavg      = eavg)
  
  # each group
  if (!is.null(group))
  {
    for (i in 1:length(group_name))
    {
      temp <- df %>% filter(group_fct == group_name[i])
      lr <- lm(outcome ~ prediction, data = temp)
      eavg <- mean(temp$er, na.rm = TRUE)
      emax <- max(temp$er, na.rm = TRUE)
      e90 <- unname(quantile(temp$er, 0.9, na.rm = TRUE))
      stats <- rbind(stats,c(group_name[i], unname(lr$coefficients[1]), unname(lr$coefficients[2]), emax, e90, eavg))
    }
  }
  
  # ------- calibration plot ------- 
  
  if (g == 1) # single group, scatter plot
  {
    pList <- list()
    
    for (i in seq_len(length(levels(df$group_fct)))) {
      res <- unlist(stats %>% filter(group == group_name[i]) %>% select(-group))
      res <- as.character(paste(paste(names(res), round(as.numeric(res), 2), sep = ": "), collapse = "\n"))
      p <- ggplot(df %>% filter(group_fct == group_name[i]) %>% droplevels(), aes(prediction, outcome))+
        geom_point(alpha = alpha)+
        geom_smooth(method = 'loess', formula = 'y ~ x', color = "red")+
        labs(title = group_name[i],
             x = "Predicted value",
             y = "Observed value",
             color = NULL)+
        xlim(min, max)+
        ylim(min, max)+
        annotate("segment", x = min, y = min, xend = max, yend = max, linetype = "dotted", color = "black")+
        annotate("text", x = labelPos.x, y = labelPos.y, label = res)+
        coord_fixed()+
        theme_bw()
      
      if (marginPlt)
      {
        p <- ggExtra::ggMarginal(p, type = "histogram")
      }
      
      pList[[i]] <- p
    }
    
    if(is.null(group)) {
      p <- pList[[1]]
    } else {
      p <- ggarrange(plotlist = pList, nrow = nrow)
    }
    
    print(p)
  } else # bin the data to box plot or violin plot
  {
    pList <- list()
    
    for (i in seq_len(length(levels(df$group_fct)))) {
      res <- unlist(stats %>% filter(group == group_name[i]) %>% select(-group))
      res <- as.character(paste(paste(names(res), round(as.numeric(res), 2), sep = ": "), collapse = "\n"))
      
      dfi <- df %>% 
        filter(group_fct == group_name[i]) %>% droplevels() %>% 
        mutate(prediction_bin = factor(ntile(prediction, n=g)))
      
      data.bin <- plyr::ddply(dfi, "prediction_bin", function(DF) {
        data.frame(median=plyr::numcolwise(median)(DF))
      })
      levels(dfi$prediction_bin) <- data.bin$median.prediction
      
      p <- ggplot(dfi,aes(x=prediction, y=outcome, group = prediction_bin))
      if (violin)
      {
        p <- p + 
          geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) 
      } else
      {
        p <- p + 
          geom_boxplot()
      }
      
      p <- p +
        geom_smooth(method = "loess", formula = 'y ~ x', color = "red", group = 1)+
        labs(title = group_name[i],
             x = "Predicted value",
             y = "Observed value",
             color = NULL)+
        xlim(min, max)+
        ylim(min, max)+
        annotate("segment", x = min, y = min, xend = max, yend = max, linetype = "dotted", color = "black")+
        annotate("text", x = min*1.1, y = max*0.9, label = res)+
        coord_fixed()+
        theme_bw()
      
      if (marginPlt)
      {
        p <- p + geom_point(alpha = 0)
        p <- ggExtra::ggMarginal(p, type = "histogram")
      }
      
      pList[[i]] <- p
    }
    
    if(is.null(group)) {
      p <- pList[[1]]
    } else {
      p <- ggarrange(plotlist = pList, nrow = nrow)
    }
    
    print(p)
  }
  
  return(stats)
}



