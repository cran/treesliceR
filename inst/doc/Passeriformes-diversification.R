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
library(sf)

## ----message = FALSE, warning = FALSE-----------------------------------------
# Loading it
library(treesliceR)

## ----eval = TRUE--------------------------------------------------------------
tree <- pass_trees[[1]]

## ----eval = TRUE--------------------------------------------------------------
ancient <- prune_tips(tree, 0.75, qtl = TRUE)
recent <- prune_tips(tree, 0.25, qtl = TRUE, method = 2)

## ----eval = TRUE--------------------------------------------------------------
oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
plot(tree, main = "Complete tree", show.tip.label = F); axisPhylo()
plot(ancient, main = "Ancient tree", show.tip.label = F); axisPhylo()
plot(recent, main = "Recent tree", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

## ----eval = TRUE--------------------------------------------------------------
anc_mat <- pass_mat[, which(colnames(pass_mat) %in% ancient$tip.label)]

## ----eval = TRUE--------------------------------------------------------------
anc_end_mat <- anc_mat[which(colSums(anc_mat) <= quantile(colSums(anc_mat), 0.30))]

## ----eval = TRUE--------------------------------------------------------------
# All recent
rec_mat <- pass_mat[, which(colnames(pass_mat) %in% recent$tip.label)]
# Endemics
rec_end_mat <- rec_mat[which(colSums(rec_mat) <= quantile(colSums(rec_mat), 0.30))]

## ----eval = TRUE--------------------------------------------------------------
AU_grid <- cbind(AU_grid, 
                 SR_anc = rowSums(anc_mat), SR_end_anc = rowSums(anc_end_mat), # Ancient
                 SR_rec = rowSums(rec_mat), SR_end_rec = rowSums(rec_end_mat)) # Recent


## ----eval = TRUE, warning = FALSE, fig.align='center', out.width="80%"--------
# Before plotting, let's set a cool palette:
pal <- colorRampPalette(c("#30009B", "#000083", "#009BFE", "#00BC00", "#FEF600", "#FE6230", "#DD0000" ))

# All species
fig_a <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_anc)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_anc)),
                       na.value="white") +
  labs(fill = expression(SR["ancient"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# Only endemics
fig_b <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_end_anc)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_end_anc)),
                       na.value="white") +
  labs(fill = expression(SR["ancient"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# To plot them together
ggarrange(fig_a, fig_b,
          labels = c("a)", "b)"), ncol = 2, nrow = 1)


## ----eval = TRUE, warning = FALSE, fig.align='center', out.width="80%"--------
# All species
fig_c <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_rec)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(1, max(AU_grid$SR_rec)),
                       na.value="white") +
  labs(fill = expression(SR["recent"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# Only endemics
fig_d <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_end_rec)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(1, max(AU_grid$SR_end_rec)),
                       na.value="white") +
  labs(fill = expression(SR["recent"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# To plot them together
ggarrange(fig_c, fig_d,
          labels = c("c)", "d)"), ncol = 2, nrow = 1)


## ----eval = TRUE--------------------------------------------------------------
head(pass_mat[, 1:4])

## ----eval = FALSE-------------------------------------------------------------
#  vec <- c(250, 500, 750, 1000, 1250, 1500, 1750, 2000)

## ----eval = FALSE-------------------------------------------------------------
#  sens_turn <- CpR_sensitivity(tree = tree, vec = vec, samp = 100,
#                               mat = pass_mat, adj = AU_adj, rate = "CpB", comp = "turnover")
#  sens_nest <- CpR_sensitivity(tree = tree, vec = vec, samp = 100,
#                               mat = pass_mat, adj = AU_adj, rate = "CpB", comp = "nestedness")

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
#  turn <- CpB(tree = tree, n = 1000, mat = pass_mat, adj = AU_adj, comp = "turnover")
#  # For nestedness component
#  nest <- CpB(tree = tree, n = 1000, mat = pass_mat, adj = AU_adj, comp = "nestedness")

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

## ----echo=TRUE----------------------------------------------------------------
# Data frames for ancient lineages:
mat_anc_rich <- as.data.frame(matrix(nrow = nrow(AU_grid), ncol = 100))
mat_anc_endemics <- as.data.frame(matrix(nrow = nrow(AU_grid), ncol = 100))

# Data frames for recent ones:
mat_rec_rich <- as.data.frame(matrix(nrow = nrow(AU_grid), ncol = 100))
mat_rec_endemics <- as.data.frame(matrix(nrow = nrow(AU_grid), ncol = 100))

## ----echo=TRUE----------------------------------------------------------------
# Create a for loop that iterate along all trees
for (i in 1:100) {
  # Select phylogeny "i" available within the package:
  tree <- pass_trees[[i]]
  # Prune the "i" phylogeny based on quantiles:
  ancient <- prune_tips(tree, 0.75, qtl = T)
  recent <- prune_tips(tree, 0.25, qtl = T, method = 2)
  
  # Capture the presence-absence matrix for ancient species:
  anc_mat <- pass_mat[, which(colnames(pass_mat) %in% ancient$tip.label)]
  mat_anc_rich[, i] <- rowSums(anc_mat) # Save their species richness
  # Capture ancient endemic richness in Australia:
  anc_mat <- anc_mat[which(colSums(anc_mat) <= quantile(colSums(anc_mat), 0.30))]
  mat_anc_endemics[, i] <- rowSums(anc_mat) # Save their species richness
  
  # Capturing the presence-absence matrix for recent species:
  rec_mat <- pass_mat[, which(colnames(pass_mat) %in% recent$tip.label)]
  mat_rec_rich[, i] <- rowSums(rec_mat) # Save their species richness
  # Capturing recent endemic richness in Australia:
  rec_mat <- rec_mat[which(colSums(rec_mat) <= quantile(colSums(rec_mat), 0.30))]
  mat_rec_endemics[, i] <- rowSums(rec_mat) # Save their species richness
}

## ----echo = TRUE--------------------------------------------------------------
# Assigning ancient species richness
AU_grid$SR_anc <- rowMeans(mat_anc_rich)
AU_grid$SR_end_anc <- rowMeans(mat_anc_endemics)
# Assigning recent species richness
AU_grid$SR_rec <- rowMeans(mat_rec_rich)
AU_grid$SR_end_rec <- rowMeans(mat_rec_endemics)


## ----eval = TRUE, warning = FALSE, fig.align='center', out.width="80%"--------
# All ancient species
fig_a <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_anc)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_anc)),
                       na.value="white") +
  labs(fill = expression(SR["ancient"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# Only ancient endemics
fig_b <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_end_anc)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_end_anc)),
                       na.value="white") +
  labs(fill = expression(SR["ancient"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# All recent species
fig_c <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_rec)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_rec)),
                       na.value="white") +
  labs(fill = expression(SR["recent"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# Only recent endemics
fig_d <- ggplot() +
  geom_sf(data = AU_grid, aes(fill = SR_end_rec)) +
  scale_fill_gradientn(colours = pal(100),
                       limits = c(0.5, max(AU_grid$SR_end_rec)),
                       na.value="white") +
  labs(fill = expression(SR["recent"])) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(.15, .14),
        legend.title = element_text(size = 13),
        legend.key.size = unit(0.75, "cm"),
        legend.direction = "horizontal")

# To plot them together
ggarrange(fig_a, fig_b, fig_c, fig_d,
          labels = c("a)", "b)", "c)", "d)"), 
          ncol = 2, nrow = 2)


## ----eval = FALSE-------------------------------------------------------------
#  CpB_turn <- lapply(pass_trees, function(x){
#    return(CpB(tree = x, n = 1000, mat = pass_mat, adj = AU_adj, comp = "turnover", ncor = 5))
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
#    return(CpB(tree = x, n = 1000, mat = pass_mat, adj = AU_adj, comp = "nestedness", ncor = 5))
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

