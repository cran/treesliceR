## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.show = "hold"
)

## ----message = FALSE, warning = FALSE-----------------------------------------
# Loading the packages that we'll use
library(ape)
library(ggplot2)
library(ggpubr)

## ----message = FALSE, warning = FALSE-----------------------------------------
# Loading it
library(treesliceR)

## ----eval = FALSE-------------------------------------------------------------
#  tree <- pass_trees[[1]]

## ----eval = FALSE-------------------------------------------------------------
#  recent <- squeeze_root(tree = tree, time = 33, dropNodes = T)
#  old <- squeeze_tips(tree = tree, time = 33, dropNodes = T)

## ----eval = FALSE-------------------------------------------------------------
#  oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
#  plot(tree, main = "Complete tree", show.tip.label = F); axisPhylo()
#  plot(old, main = "Old tree", show.tip.label = F); axisPhylo()
#  plot(recent, main = "Recent tree", show.tip.label = F); axisPhylo()
#  par(oldpar) # Returning to the original display

## ----echo = FALSE, fig.align = 'center'---------------------------------------
knitr::include_graphics("Sliced_pas_tree.png")

## ----eval = FALSE-------------------------------------------------------------
#  DR_diff <- DR(tree = recent)[,2] - DR(tree = old)[,2]

## ----eval = FALSE-------------------------------------------------------------
#  tree$Species_info$DR_diff <- DR_diff

## ----eval = FALSE-------------------------------------------------------------
#  # tapply() for means
#  fam_DR <- tapply(tree$Species_info$DR_diff, tree$Species_info$Family, mean)
#  # tapply() for standard deviations
#  fam_DR_sd <- tapply(tree$Species_info$DR_diff, tree$Species_info$Family, sd)

## ----eval = FALSE-------------------------------------------------------------
#  # Creating the families DR data frame
#  fam_df <- data.frame(Family = names(fam_DR), DR_diff = fam_DR, DR_sd = fam_DR_sd)
#  # Sorting them based on DR's value
#  fam_df <- fam_df[order(fam_df$DR_diff),]
#  # Turning this order into a factor to plot it
#  fam_df$Family <- factor(fam_df$Family, levels = fam_df$Family)

## ----eval = FALSE-------------------------------------------------------------
#  ggplot(fam_df, aes(x = Family, y = DR_diff,
#                      ymin = DR_diff - DR_sd,
#                      ymax = DR_diff + DR_sd)) +
#    geom_pointrange(color = "#d90429") +
#    geom_hline(yintercept = 0, linetype="dashed", color = "black") +
#    coord_flip() + theme_minimal() +
#    theme(axis.title = element_text(size = 13),
#          axis.text = element_text(size = 10),
#          axis.line = element_line(colour = "black"),
#          panel.grid.major.x = element_blank(),
#          panel.grid.minor.x = element_blank()) +
#    ylab(expression(paste(DR["recent"]-DR["past"]))) + xlab(NULL)

## ----echo = FALSE, fig.align='center', out.width="60%"------------------------
knitr::include_graphics("OneTree_DR.png")

## ----eval = TRUE--------------------------------------------------------------
head(pass_asb[[1]][, 1:4])

## ----eval = FALSE-------------------------------------------------------------
#  vec <- c(250, 500, 750, 1000, 1250, 1500, 1750, 2000)

## ----eval = FALSE-------------------------------------------------------------
#  sens_turn <- CpR_sensitivity(tree = tree, vec = vec, samp = 100,
#                               asb = pass_asb, rate = "CpB", comp = "turnover")
#  sens_nest <- CpR_sensitivity(tree = tree, vec = vec, samp = 100,
#                               asb = pass_asb, rate = "CpB", comp = "nestedness")

## ----eval = FALSE-------------------------------------------------------------
#  # Store each graph within a respective object
#  turn_sens_plot <- CpR_sensitivity_plot(sens_turn, rate = "CpB", stc = "mean") +
#    geom_vline(xintercept = 1000, linetype="dashed", color = "black")
#  nest_sens_plot <- CpR_sensitivity_plot(sens_nest, rate = "CpB", stc = "mean") +
#    geom_vline(xintercept = 1000, linetype="dashed", color = "black")
#  # To plot them together
#  ggarrange(turn_sens_plot, nest_sens_plot,
#                      labels = c("a)", "b)"), ncol = 2, nrow = 1)
#  

## ----echo = FALSE, fig.align='center', out.width="80%"------------------------
knitr::include_graphics("OneTree_sensitivity.png")

## ----eval = FALSE-------------------------------------------------------------
#  # For turnover component
#  turn <- CpB(tree = tree, n = 1000, asb = pass_asb, comp = "turnover")
#  # For nestedness component
#  nest <- CpB(tree = tree, n = 1000, asb = pass_asb, comp = "nestedness")

## ----ggplot2, eval = FALSE----------------------------------------------------
#  turn_1 <- CpR_graph(data = turn, rate = "CpB", qtl = TRUE)
#  turn_2 <- CpR_graph(data = turn, rate = "CpB", qtl = TRUE, map = AU_grid)
#  # To plot them together
#  ggarrange(turn_1, turn_2,
#                      labels = c("a)", "b)"), ncol = 2, nrow = 1)

## ----echo = FALSE, fig.align='center'-----------------------------------------
knitr::include_graphics("OneTree_turn.png")

## ----eval = FALSE-------------------------------------------------------------
#  nest1 <- CpR_graph(data = nest, rate = "CpB", qtl = TRUE)
#  nest2 <- CpR_graph(data = nest, rate = "CpB", qtl = TRUE, map = AU_grid)
#  # To plot them together
#  ggarrange(nest1, nest2,
#                      labels = c("a)", "b)"), ncol = 2, nrow = 1)

## ----echo = FALSE, fig.align='center'-----------------------------------------
knitr::include_graphics("OneTree_nest.png")

## ----eval = FALSE-------------------------------------------------------------
#  recent <- lapply(tree, function(x){return(squeeze_root(x, 33, dropNodes = T))})
#  old <- lapply(tree, function(x){return(squeeze_tips(x, 33, dropNodes = T))})

## ----eval = FALSE-------------------------------------------------------------
#  DR_rec <- lapply(recent, function(x){DR(x)})
#  DR_old <- lapply(old, function(x){DR(x)})

## ----eval = FALSE-------------------------------------------------------------
#  # Recent
#  f_DRrec <- DR_rec[[1]]
#  colnames(f_DRrec)[2] <- 1
#  # Old
#  f_DRold <- DR_old[[1]]
#  colnames(f_DRold)[2] <- 1
#  
#  # Looping
#  for (i in 2:length(DR_rec)) {
#    # Recent
#    f_DRrec <- merge(f_DRrec, DR_rec[[i]], by = "Species", sort = FALSE)
#    colnames(f_DRrec)[i + 1] <- i
#    # Old
#    f_DRold <- merge(f_DRold, DR_old[[i]], by = "Species", sort = FALSE)
#    colnames(f_DRold)[i + 1] <- i
#  }

## ----eval = FALSE-------------------------------------------------------------
#  # Recent
#  DR_rec_mean <- apply(f_DRrec[, -1], 1, mean)
#  DR_rec_sd <- apply(f_DRrec[, -1], 1, sd)
#  
#  # Old
#  DR_old_mean <- apply(f_DRold[, -1], 1, mean)
#  DR_old_sd <- apply(f_DRold[, -1], 1, sd)
#  
#  # Capturing their mean difference
#  DR_diff <- DR_rec_mean - DR_old_mean

## ----eval = FALSE-------------------------------------------------------------
#  df <- data.frame(Specie = f_DRrec[, 1], DR_diff = DR_diff)
#  # Capturing families information
#  df <- merge(df, tree[[1]]$Species_info, by = "Specie", sort = FALSE)
#  
#  ## Displaying it graphically for families
#  fam_DR <- tapply(df$DR_diff, df$Family, mean)
#  fam_DR_sd <- tapply(df$DR_diff, df$Family, sd)

## ----eval = FALSE-------------------------------------------------------------
#  fam_df <- data.frame(Family = names(fam_DR), DR_diff = fam_DR, DR_sd = fam_DR_sd)
#  # Ordering to plot
#  fam_df <- fam_df[order(fam_df$DR_diff),]
#  fam_df$Family <- factor(fam_df$Family, levels = fam_df$Family)

## ----eval = FALSE-------------------------------------------------------------
#  ggplot(fam_df, aes(x = Family, y = DR_diff,
#                      ymin = DR_diff - DR_sd,
#                      ymax = DR_diff + DR_sd)) +
#    geom_pointrange(color = "#d90429") +
#    geom_hline(yintercept = 0, linetype="dashed", color = "black") +
#    coord_flip() + theme_minimal() +
#    theme(axis.title = element_text(size = 13),
#          axis.text = element_text(size = 10),
#          axis.line = element_line(colour = "black"),
#          panel.grid.major.x = element_blank(),
#          panel.grid.minor.x = element_blank()) +
#    ylab(expression(paste(DR["recent"]-DR["past"]))) + xlab(NULL)

## ----echo = FALSE, fig.align='center', out.width="60%"------------------------
knitr::include_graphics("AllTrees_DR.png")

## ----eval = FALSE-------------------------------------------------------------
#  CpB_turn <- lapply(pass_trees, function(x){
#    return(CpB(tree = x, n = 1000, asb = pass_asb, comp = "turnover", ncor = 5))
#  })

## ----eval = FALSE-------------------------------------------------------------
#  CpB_val <- sapply(CpB_turn, function(x){return(x[,1])})
#  pB_val <- sapply(CpB_turn, function(x){return(x[,2])})
#  pBO_val <- sapply(CpB_turn, function(x){return(x[,3])})
#  # Creating the new data frame
#  turn_100trees <- data.frame(CpB = apply(CpB_val, 1, mean),
#                              pB = apply(pB_val, 1, mean),
#                              pBO = apply(pBO_val, 1, mean))

## ----eval = FALSE-------------------------------------------------------------
#  turn_1 <- CpR_graph(data = turn_100trees, rate = "CpB", qtl = TRUE)
#  turn_2 <- CpR_graph(data = turn_100trees, rate = "CpB", qtl = TRUE, map = AU_grid)
#  # To plot them together
#  ggarrange(turn_1, turn_2,
#            labels = c("a)", "b)"), ncol = 2, nrow = 1)

## ----echo = FALSE, fig.align='center'-----------------------------------------
knitr::include_graphics("AllTrees_turn.png")

## ----eval = FALSE-------------------------------------------------------------
#  CpB_nest <- lapply(pass_trees, function(x){
#    return(CpB(tree = x, n = 1000, asb = pass_asb, comp = "nestedness", ncor = 5))
#  })

## ----eval = FALSE-------------------------------------------------------------
#  CpB_val <- sapply(CpB_nest, function(x){return(x[,1])})
#  pB_val <- sapply(CpB_nest, function(x){return(x[,2])})
#  pBO_val <- sapply(CpB_nest, function(x){return(x[,3])})
#  # DF
#  nest_100trees <- data.frame(CpB = apply(CpB_val, 1, mean),
#                              pB = apply(pB_val, 1, mean),
#                              pBO = apply(pBO_val, 1, mean))

## ----eval = FALSE-------------------------------------------------------------
#  nest_1 <- CpR_graph(data = nest_100trees, rate = "CpB", qtl = TRUE)
#  nest_2 <- CpR_graph(data = nest_100trees, rate = "CpB", qtl = TRUE, map = AU_grid)
#  # To plot them together
#  ggarrange(nest_1, nest_2,
#            labels = c("a)", "b)"), ncol = 2, nrow = 1)

## ----echo = FALSE, fig.align='center'-----------------------------------------
knitr::include_graphics("AllTrees_nest.png")

