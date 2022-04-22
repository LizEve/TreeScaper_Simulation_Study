#https://thibautjombart.github.io/treespace/articles/introduction.html

library("treespace")
library("phangorn")
library("adegenet")
library("adegraphics")
library("rgl")

# generate list of trees
suppressWarnings(RNGversion("3.5.0"))
#set.seed(1)
#x <- rmtree(10, 20)
#names(x) <- paste("tree", 1:10, sep = "")
setwd('/Users/chatnoir/Projects/TreeScaper_ExampleData/Parks_2017_gene_trees/gene_trees')
#x<-read.tree('trees.bootstrap._1.inclade1.ortho1.occ_0.2/RAxML_bootstrap.OG0005948_1.inclade1.ortho1.occ_0.2.all_bootstrap.tre', keep.multi=TRUE)
#y<-read.tree('trees.bootstrap._1.inclade1.ortho1.occ_0.2/RAxML_bootstrap.OG0006169_1.inclade1.ortho1.occ_0.2.all_bootstrap.tre', keep.multi=TRUE)
#z <- c(x, y)
z2 <- read.tree('occ0.2sp100/all.0.2.sub10.tre') #320
z5 <- read.tree('occ0.5sp100/all.0.5.sub10.tre') #640
z8 <- read.tree('occ0.8sp100/all.0.8.sub10.tre') #960
#r2 <- sample(z2,size=200)
#r5 <- sample(z5,size=200)
#r8 <- sample(z8,size=200)

z <- c(z2,z5,z8)


rt <- root(z,"hemi",resolve.root = TRUE)
rt
#rand <- sample(rt,size=500)

#write.tree(rt, "Park.sub10.tree")
#write.tree(z, "Park.sub10.unroot.tree")

# Label where trees came from 

names(rt)[1:320] <- paste0("oh2",1:320)
names(rt)[321:640] <- paste0("oh5",1:320)
names(rt)[641:960] <- paste0("oh8",1:320)

#names(rt)[1:100] <- paste0("one",1:100)
#names(rt)[101:200] <- paste0("two",1:100)
Dtype <- c(rep("oh2",320),rep("oh5",320),rep("oh8",320))


# create vector corresponding to tree inference method:
#Dtype <- c(rep("oh2",200),rep("oh5",200),rep("oh8",200))

# use treespace to find and project the distances: Uses kenji colin automatically
Dscape <- treespace(rt, nf=3, method="RF")


plotGrovesD3(Dscape$pco, groups=Dtype)



# use treespace
res <- treespace(rt, nf=3)
names(res)

# table.image
#table.image(res$D, nclass=30)

# table.value with some customization
#table.value(res$D, nclass=6, method="color", 
#            symbol="circle", col=redpal(6))

# 2-dimensional space is given by the first 2 PCs of the MDS.
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)


# Find groups
wm.groves <- findGroves(res, nclust=3)
#names(wm.groves)
plotGrovesD3(wm.groves)

plotGrovesD3(wm.groves, xax=2, yax=3)

colours <- fac2col(wm.groves$groups, col.pal=funky)
plot3d(wm.groves$treespace$pco$li[,1],
       wm.groves$treespace$pco$li[,2],
       wm.groves$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")

