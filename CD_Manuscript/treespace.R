#https://thibautjombart.github.io/treespace/articles/introduction.html

library("treespace")
library("phangorn")
library("adegenet")
library("adegraphics")
library("rgl")

# generate list of trees
suppressWarnings(RNGversion("3.5.0"))
set.seed(1)
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")

# use treespace
res <- treespace(x, nf=3)
names(res)

# table.image
table.image(res$D, nclass=30)

# table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))

# 2-dimensional space is given by the first 2 PCs of the MDS.
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)


data(woodmiceTrees)
wm.res <- treespace(woodmiceTrees,nf=3)

# PCs are stored in:
head(wm.res$pco$li)
