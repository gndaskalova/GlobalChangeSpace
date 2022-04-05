load("combined_pop_marine_simple.RData")

pop_mar_pca <- prcomp(combined_pop_marine_simple[, -6], scale = TRUE,
                      center = TRUE)

library(ggbiplot)
devtools::install_github("andrewljackson/SIBER@master",
                         build_vignettes = FALSE)
library(SIBER)

(g <- ggbiplot(pop_mar_pca, obs.scale = 1, var.scale = 1, 
               groups = combined_pop_marine_simple$sampling2, ellipse = TRUE,
               circle = TRUE))

# create the siber object to calculate overlap
# this package was made for stable isotope analysis and it wants the 
# object to have specific column names, like iso1 and iso2
ind <- get_pca_ind(pop_mar_pca)
test <- as.data.frame(ind$coord[,1:2])
test$group <- NA
test$group <- combined_pop_marine_simple$sampling2
test$community <- "1"

colnames(test) <- c("iso1", "iso2", "group", "community")
siber.example <- createSiberObject(test)

group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
test$group <- factor(test$group, levels = c("Marine Random sampling", "Marine Living Planet Database"),
                     labels = c("2", "3"))
test <- test %>% dplyr::select(iso1, iso2, group, community)
colnames(test) <- c("iso1", "iso2", "group", "community")

ellipse1 <- "1.2"

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "1.3"

siber.example <- createSiberObject(test)

# The overlap of the maximum likelihood fitted standard ellipses are 
# estimated using
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                                   p.interval = 0.95, n = 100)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
prop.95.over
# 69.9% overlap between LPD marine and global sampling


library(hagis)
library(vegan)
library(ape)
pathotype.anosim <- anosim(combined_pop_marine_simple, combined_pop_marine_simple$sampling2)


# PCA
iris_c <- scale(combined_pop_marine_simple[, -6])
pca <- rda(iris_c)

# plot
plot(pca, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green')
points(pca, display='sites', col = cols[combined_pop_marine_simple$sampling2], pch = 16)
ordihull(pca, groups=combined_pop_marine_simple$sampling2)
ordispider(pca, groups = combined_pop_marine_simple$sampling2, label = TRUE)

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis(iris_c ~ sampling2, data = combined_pop_marine_simple, method='eu')


install.packages("Momocs")
library(Momocs)

pca <- PCA(efourier(combined_pop_marine_simple, 12))
get_chull_area(pca)


pc <- prcomp(combined_pop_marine_simple[, -6], scale = TRUE,
            center = TRUE)
# Plot PC points for PC1 and PC2 and overlay convex hull:
plot(pc$x[,1:2])
hpts <- chull(pc$x[,1:2]) # Identify vertices of convex hull
test3 <- polygon(pc$scores[hpts, 1:2])

test <- chull(pc$x[,1:2])
test2 <- chull(pc$x[,2])

library(sp)

chull.poly <- Polygon(hpts, hole=F)
chull.area <- chull.poly@area