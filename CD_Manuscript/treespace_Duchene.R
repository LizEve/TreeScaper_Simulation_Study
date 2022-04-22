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
setwd('/Users/chatnoir/Projects/TreeScaper_ExampleData/Duchene_2017_doi_10.5061_dryad.353q5__v1/Phylogenomics')
z <- read.tree('all.loci.trees')

# Keep names http://blog.phytools.org/2011/12/dropping-same-tip-from-many-trees-in.html

keep <- c("Peramelidae_Isoodon_obesulus", "Burramyidae_Burramys_parvus", "Thylacomyidae_Macrotis_lagotis", "Notoryctidae_Notoryctes_typhlops", "Vombatidae_Lasiorhinus_latifrons", "Potoroidae_Bettongia_penicillata", "Potoroidae_Bettongia_lesueur", "Dasyuridae_Dasyurus_maculatus", "Phalangeridae_Trichosurus_vulpecula", "Phalangeridae_Spilocuscus_maculatus", "Potoroidae_Aepyprymnus_rufescens", "Pseudochiridae_Pseudocheirus_peregrinus_north", "Peramelidae_Perameles_nasuta", "Myrmecobiidae_Myrmecobius_fasciatus", "Phalangeridae_Strigocuscus_pelengensis", "Dasyuridae_Phascogale_tapoatafa", "Vombatidae_Vombatus_ursinus", "Dasyuridae_Sminthopsis_macroura", "Phalangeridae_Phalanger_orientalis", "Acrobatidae_Acrobates_pygmaeus", "Potoroidae_Potorous_tridactylus_north", "Dasyuridae_Phascogale_calura", "Dasyuridae_Sarcophilus_harrisii", "Dasyuridae_Planigale_tenuirostris", "Peramelidae_Echymipera_kalubu", "Macropodidae_Dorcopsulus_vanheurni", "Pseudochiridae_Pseudochirulus_forbesi", "Didelphidae_Didelphis_virginiana", "Marmosidae_Monodelphis_domestica", "Dasyuridae_Antechinus_stuartii", "Macropodidae_Onychogalea_fraenata", "Hypsiprymnodontidae_Hypsiprymnodon_moschatus", "Microbiotheriidae_Dromiciops_gliroides", "Pseudochiridae_Petauroides_volans", "Macropodidae_Wallabia_bicolor", "Pseudochiridae_Pseudocheirus_occidentalis", "Tarsipedidae_Tarsipes_rostratus", "Acrobatidae_Distoechurus_pennatus", "Burramyidae_Cercartetus_nanus_northern", "Phalangeridae_Phalanger_carmelitae", "Peramelidae_Peroryctes_raffrayana", "Pseudochiridae_Pseodochirops_corinnae", "Dasyuridae_Dasyurus_geffroii", "Peramelidae_Perameles_gunnii", "Phascolarctidae_Phascolarctos_cinereus")
  
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))


ntrees<-lapply(z,keep.tip,tip=keep)
"multiPhylo"->class(ntrees)

# use treespace
res <- treespace(ntrees, nf=3)
names(res)

# table.image
#table.image(res$D, nclass=30)

# table.value with some customization
#table.value(res$D, nclass=6, method="color", 
#            symbol="circle", col=redpal(6))

# 2-dimensional space is given by the first 2 PCs of the MDS.
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)


# Find groups
wm.groves <- findGroves(res, nclust=4)
names(wm.groves)

plotGrovesD3(wm.groves)

plotGrovesD3(wm.groves, xax=2, yax=3)

colours <- fac2col(wm.groves$groups, col.pal=funky)
plot3d(wm.groves$treespace$pco$li[,1],
       wm.groves$treespace$pco$li[,2],
       wm.groves$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")

