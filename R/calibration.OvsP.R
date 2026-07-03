require(dplyr)
require(tibble)
require(stringr)
require(ggplot2)
require(ggpubr) # `ggarrange`

#' Subgroup OvsP plots
#'
#' @param res Returned result from `calibration`
#' @param group The vector of group labels if want to assess calibration within subgroups
#' @param group_con Flag if the group variable is continuous
#' @param group_name A string vector, levels of group
#' @param xlab Title of x-axis in the subgroup OvsP plot
#' @export
calibration.OvsP <- function(
  res,
  group,
  group_con,
  group_name = NULL,
  xlab = "Group"
) {
  res <- as.data.frame(res)

  if (is.null(res$group)) {
    res <- rownames_to_column(res, var = "group")
  }
  res <- res %>% filter(group != "Overall")
  
  if (!group_con) {
    res$group <- factor(res$group)
    if (!is.null(group_name)) {
      levels(res$group) <- group_name
    }

    p <- ggplot(res, aes(group, Obs)) +
      geom_col(fill = "white", color = "black") +
      geom_point(
        aes(group, Pavg, fill = "Predicted average"),
        shape = 23,
        size = 3
      ) +
      scale_fill_manual(values = c("Predicted average" = "blue")) +
      theme_bw() +
      theme(
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "top"
      ) +
      labs(x = xlab, y = "Observed average")
    print(p)
  } else {
    # extract bin boundaries and centers
    bin_bry <- as.numeric(unique(unlist(str_extract_all(
      res$group,
      "\\d+\\.*\\d*"
    ))))
    bin_ctr <- rep(0, nrow(res))
    for (i in 1:nrow(res)) {
      bin_ctr[i] <- (bin_bry[i] + bin_bry[i + 1]) / 2
    }
    res$bin_ctr <- bin_ctr

    # bar plot
    p1 <- ggplot(res, aes(bin_ctr, Obs)) +
      geom_col(fill = "white", color = "black") +
      geom_point(
        aes(bin_ctr, Pavg, fill = "Predicted average"),
        shape = 23,
        size = 3
      ) +
      scale_x_continuous(breaks = res$bin_ctr, labels = res$group) +
      scale_fill_manual(values = c("Predicted average" = "blue")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "top"
      ) +
      labs(x = xlab, y = "Observed average")

    # density plot
    p2 <- ggplot(data.frame(y = group), aes(y)) +
      geom_density() +
      theme_bw()
    print(ggarrange(p1, p2, nrow = 2, heights = c(1, 0.4)))
  }
}