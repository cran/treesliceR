## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.show = "hold"
)

## ----message = FALSE, warning = FALSE, eval=TRUE------------------------------
# Loading it
library(ape)

## ----message = FALSE, warning = FALSE, eval=TRUE------------------------------
# Loading it
library(treesliceR)

## ----readData, echo=TRUE, eval=TRUE-------------------------------------------
tree <- pass_trees[[1]]

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_pruned10 <- prune_tips(tree = tree, time = 10, qtl = F) # keep species older than 10my
tree_pruned30 <- prune_tips(tree = tree, time = 30, qtl = F) # keep species older than 30my

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
plot(tree, main = "All species", show.tip.label = F); axisPhylo()
plot(tree_pruned10, main = "Species older than 10my", show.tip.label = F); axisPhylo()
plot(tree_pruned30, main = "Species older than 30my", show.tip.label = T); axisPhylo()
par(oldpar) # Returning to the original display


## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_pruned10_inverse <- prune_tips(tree = tree, time = 10, qtl = F, method = 2) # keep species younger than 10my
tree_pruned30_inverse <- prune_tips(tree = tree, time = 30, qtl = F, method = 2) # keep species younger than 30my

# plotting phylogenies
oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
plot(tree, main = "All species", show.tip.label = F); axisPhylo()
plot(tree_pruned10_inverse, main = "Species younger than 10my", show.tip.label = F); axisPhylo()
plot(tree_pruned30_inverse, main = "Species younger than 30my", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_pruned25q <- prune_tips(tree, 0.25, qtl = T, method = 2)
plot(tree_pruned25q, main = "Species with ages younger than 25th quantile", show.tip.label = F); axisPhylo()

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_squeeze10 <- squeeze_tips(tree = tree, time = 10)
tree_squeeze30 <- squeeze_tips(tree = tree, time = 30)
tree_squeeze50 <- squeeze_tips(tree = tree, time = 50)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
oldpar <- par(mfrow = c(1, 3)) # Setting an 1x3 graphical display
plot(tree_squeeze10, main = "squeezed at 10my", show.tip.label = F); axisPhylo()
plot(tree_squeeze30, main = "Squeezed at 30my", show.tip.label = F); axisPhylo()
plot(tree_squeeze50, main = "Squeezed at 50my", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_squeeze30_drop <- squeeze_tips(tree = tree, time = 30, criteria = "my", dropNodes = TRUE)
tree_squeeze30 # full binary tree
tree_squeeze30_drop # tree with nodes dropped

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
PD_total <- sum(tree$edge.length)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
tree_squeeze10 <- squeeze_tips(tree = tree, time = PD_total/10, criteria = "pd")
tree_squeeze50 <- squeeze_tips(tree = tree, time = PD_total/2, criteria = "pd")
oldpar <- par(mfrow = c(1, 2)) # Setting an 1x2 graphical display
plot(tree_squeeze10, main = "Tree squeezed at 10% PD", show.tip.label = F); axisPhylo()
plot(tree_squeeze50, main = "Tree squeezed at 50% PD", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

## ----squeezeRoot, echo=TRUE, eval=TRUE----------------------------------------
tree_root50 <- squeeze_root(tree = tree, time = 50)
plot(tree_root50, main = "Tree sliced rootwardly in 50my", show.tip.label = F); axisPhylo()

## ----echo=TRUE,eval=TRUE------------------------------------------------------
tree_int <- squeeze_int(tree = tree, from = 30, to = 10,)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
oldpar <- par(mfrow = c(1, 2)) # Setting an 1x2 graphical display
plot(tree, main = "Original tree", show.tip.label = F); axisPhylo()
plot(tree_int, main = "Tree slice 10-30my", show.tip.label = F); axisPhylo()
par(oldpar) # Returning to the original display

