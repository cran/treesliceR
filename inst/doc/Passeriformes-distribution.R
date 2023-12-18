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

## -----------------------------------------------------------------------------
tree <- pass_trees[[1]]

## -----------------------------------------------------------------------------
older <- prune_tips(tree, 0.75, qtl = T)

## -----------------------------------------------------------------------------
younger <- prune_tips(tree, 0.25, qtl = T, method = 2)

## -----------------------------------------------------------------------------
oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
plot(older, main = "All species", show.tip.label = F); axisPhylo()
plot(older, main = "Older species", show.tip.label = F); axisPhylo()
plot(younger, main = "Younger species", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

## -----------------------------------------------------------------------------
older <- older$tip.label
younger <- younger$tip.label

## -----------------------------------------------------------------------------
head(pass_mat[, 1:3])

## -----------------------------------------------------------------------------
positions <- which(colnames(pass_mat) %in% older)
old_mat <- pass_mat[, positions]

## -----------------------------------------------------------------------------
positions <- which(colnames(pass_mat) %in% younger)
yng_mat <- pass_mat[, positions]

## -----------------------------------------------------------------------------
grid <- AU_grid
grid$SR_old <- rowSums(old_mat)
grid$SR_yng <- rowSums(yng_mat)

## -----------------------------------------------------------------------------
# Younger
AU_yng <- ggplot() +
          geom_sf(data = grid, aes(fill = SR_yng)) +
          scale_fill_gradient(low = "white", high = "red") +          
          labs(fill = "SR (Young)") +
          theme_void()
# Older
AU_old <- ggplot() +
          geom_sf(data = grid, aes(fill = SR_old)) +
          scale_fill_gradient(low = "white", high = "red") +
          labs(fill = "SR (Old)") +
          theme_void()
# Plotting
ggarrange(AU_yng, AU_old, labels = c("a)", "b)"), ncol = 2, nrow = 1)

## -----------------------------------------------------------------------------
model <- lm(SR_yng ~ SR_old, data = grid)

## -----------------------------------------------------------------------------
summary(model)

## -----------------------------------------------------------------------------
plot(SR_yng ~ SR_old, data = grid,
     ylab = "Species richness (younger)",
     xlab = "Species richness (older)")
# Regression line
abline(model, col = "red")

## -----------------------------------------------------------------------------
grid$res <- residuals(model)
# Mapping it
ggplot() + 
  geom_sf(data = grid, aes(fill = res)) +
  scale_fill_gradient2(low = "red", mid = "grey80", high = "blue") +
  labs(fill = "Residuals") +
  theme_void()

