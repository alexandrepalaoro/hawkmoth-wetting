library(treeio)
library(ape)
library(phytools)
library(geiger)
library(phylolm)
library(scales)


## Loading data
wet = read.csv("data/hawkmoth-avgs.csv", h=T,sep = ';')

# Sometimes my excel adds weird characters. Just fixing it
colnames(wet)[1] = "species"

## Loading phylogeny
phy = read.beast("phylo/Hawkmoth.202tx.BEAST.tre")

# Getting just the species that are present in the phylo
wet.phy1 = wet[!is.na(wet$sp.in.tree),]

# pruning
pruned = keep.tip(phy@phylo,wet.phy1$sp.in.tree)

plot(pruned);axisPhylo()
nodelabels();edgelabels();tiplabels()

## Adding the missing species manually
pruned1 = bind.tip(pruned,"E fasciatus", where = 4, position = 10)

plot(pruned1)

wet$species2 = gsub("_", " ", wet$species)
pruned1$tip.label = c("Agrius cingulata", "Paratrea plebeja","Dolba hyloeus",
                      "Manduca rustica", "Manduca quinquemaculata",
                      "Manduca sexta", "Enyo lugubris", "Eumorpha pandorus",
                      "Eumorpha fasciatus", "Hemaris thysbe",  "Darapsa myron",  
                      "Hyles lineata", "Xylophanes tersa")


rownames(wet) = wet$species2

name.check(pruned1,wet)


### Adding a distinction between the two subfamilies in our dataset

cols<-c("grey","orange"); names(cols)<-1:2

plot(pruned1)
nodelabels()

pruned2 = paintSubTree(pruned1, node=15, state='2')
plot(pruned2)

dev.off()

#### ADVANCING CONTACT ANGLE MODELS OF EVOLUTION ####


m1.BM = phylolm(avg.ca ~ prob.wing , data = wet, phy = pruned2, model = "BM", measurement_error = T)
summary(m1.BM)

m1.L = phylolm(avg.ca ~ prob.wing , data = wet, phy = pruned2, model = "lambda", lower.bound = 1e-10, upper.bound = 1)
summary(m1.L)

m1.K = phylolm(avg.ca ~ prob.wing , data = wet, phy = pruned2, model = "kappa", measurement_error = T)
summary(m1.K)

m1.d = phylolm(avg.ca ~ prob.wing , data = wet, phy = pruned2, model = "delta", upper.bound = 180, measurement_error =  T)
summary(m1.d)

AIC(m1.BM)
AIC(m1.L)
AIC(m1.K)
AIC(m1.d)


#### BM IS THE BEST, LET's BOOTSTRAP IT ####

m1.BM = phylolm(avg.ca ~ prob.wing , data = wet, phy = pruned2, model = "BM", measurement_error = T, boot = 10000)
summary(m1.BM)

##### SAME MODELS BUT NOW USING THE SLOPE OF CHANGE #### 

m2.BM = phylolm(slope.ca ~ prob.wing , data = wet, phy = pruned2, model = "BM", measurement_error = T)
summary(m2.BM)

m2.L = phylolm(slope.ca ~ prob.wing , data = wet, phy = pruned2, model = "lambda", lower.bound = 1e-10, upper.bound = 1)
summary(m2.L)

m2.K = phylolm(slope.ca ~ prob.wing , data = wet, phy = pruned2, model = "kappa", measurement_error = T)
summary(m2.K)

m2.d = phylolm(slope.ca ~ prob.wing , data = wet, phy = pruned2, model = "delta", upper.bound = 180, measurement_error =  T)
summary(m2.d)


AIC(m2.BM)
AIC(m2.L)
AIC(m2.K)
AIC(m2.d)

m2.BM$aic - m2.d$aic

##### BOOTSTRAPING ####

m2.BM = phylolm(slope.ca ~ prob.wing , data = wet, phy = pruned2, model = "BM", measurement_error = T, boot = 10000)
summary(m2.BM)


##### PHENOGRAMS ####

avg.ca = setNames(wet$avg.ca, rownames(wet))
slope.ca = setNames(wet$slope.ca, rownames(wet))

phenogram(pruned2,avg.ca,colors=cols, ftype = 'i', ylim = c(40,80),
          ylab = "Average advancing contact angle (degrees)", xlab = "Time since root (my)")

phenogram(pruned2,slope.ca,colors=cols, ftype = 'i', 
          ylab = "Slope of change in contact angle along the proboscis", xlab = "Time since root (my)")


### FIGURE ####

## I added the numbers of the phenograms in photoshop. That's why I exported it to PDF.
## PDFs are better to work in photoshop (because they are vector graphs).

pdf("figures/pheno-ca-probwing.pdf",w=14,h=12)

par(mfrow=c(2,2), las = 1, bty = 'l')

plot(avg.ca ~ prob.wing, data = wet, cex = 2, bty = 'l', las = 1,
     pch = 21, bg = c("grey","orange")[as.numeric(as.factor(wet$group))],
     xlab = "Ratio between proboscis length and forewing length",
     ylab = "Average advancing contact angles (degrees)",
     ylim = c(40,90))
abline(m1.BM, lwd = 5, lty = 2, col = alpha('black',0.7))

arrows(x0=wet$prob.wing, x1=wet$prob.wing,
       y0=wet$avg.ca-wet$sd.ca, y1=wet$avg.ca+wet$sd.ca, 
       code=3, angle=90, length=0.1, lwd = 2)

points( x = wet$prob.wing, y = wet$avg.ca, cex = 2.5,
        pch = 21, bg = c("grey50","orange")[as.numeric(as.factor(wet$group))])

text(avg.ca ~ prob.wing, data = wet, labels = lab, font = 2, pos = 4.5, cex = 0.9)

legend("topright",legend="(A)",bty='n',cex=1.2)

phenogram(pruned2,avg.ca,colors=cols, ftype = 'i', ylim = c(40,80),
          ylab = "Average advancing contact angle (degrees)", xlab = "Time since root (my)")

legend("topright",legend="(B)",bty='n',cex=1.2)


plot(slope.ca ~ prob.wing, data = wet, cex = 2.5, bty = 'l', las = 1,
     pch = 21, bg = c("grey","orange")[as.numeric(as.factor(wet$group))],
     xlab = "Ratio between proboscis length and forewing length",
     ylab = "Slope of change in contact angle along the proboscis")

abline(m2.BM, lwd = 5, lty = 2, col = alpha('black',0.7))

text(slope.ca ~ prob.wing, data = wet, labels = lab, font = 2, pos = 4.5, cex = 0.9)

legend("topright",legend="(C)",bty='n',cex=1.2)

phenogram(pruned2,slope.ca,colors=cols, ftype = 'i', 
          ylab = "Slope of change in contact angle along the proboscis", xlab = "Time since root (my)")
legend("topright",legend="(D)",bty='n',cex=1.2)

dev.off()


#### ANALYSES OF THE CURVATURE ####

m3.BM = phylolm(curv.food ~ avg.ca , data = wet, phy = pruned2, model = "BM", measurement_error = T)
summary(m3.BM)

m3.L = phylolm(curv.food ~ avg.ca , data = wet, phy = pruned2, model = "lambda", lower.bound = 1e-10, upper.bound = 1)
summary(m3.L)

m3.K = phylolm(curv.food ~ avg.ca , data = wet, phy = pruned2, model = "kappa", measurement_error = T)
summary(m3.K)

m3.d = phylolm(curv.food ~ avg.ca , data = wet, phy = pruned2, model = "delta", upper.bound = 180, measurement_error =  T)
summary(m3.d)


AIC(m3.BM)
AIC(m3.L)
AIC(m3.K)
AIC(m3.d)

##### BM IS BEST, BOOTSTRAPPING ####

m3.BM = phylolm(curv.food ~ avg.ca , data = wet, phy = pruned2, model = "BM", measurement_error = T, boot = 10000)
summary(m3.BM)

### FIGURE ####

pdf("figures/curv-food-graphs.pdf",w=14,h=5)
par(mfrow = c(1,3), las = 1, bty = 'l')

plot(radius.curv ~ radius.food, data = wet, bty = 'l', las = 1, cex = 2.5,
     pch = 21, bg = c("grey","orange")[as.numeric(as.factor(wet$group))],
     ylab = "Radius of curvature of the galea",
     xlab = "Radius of curvature of the food canal",
     ylim = c(0,0.14), xlim = c(0,0.14))
abline(0,1,lwd = 2, lty = 2)

text(radius.curv ~ radius.food, data = wet, labels = lab, font = 2, pos = 3, cex = 0.9)

legend("topleft",legend="(A)",bty='n',cex=1.2)

plot(curv.food ~ avg.ca, data = wet, cex = 2.5, bty = 'l', las = 1,
     pch = 21, bg = c("grey","orange")[as.numeric(as.factor(wet$group))],
     ylab = "Ratio of radii of curvature (outer/inner)",
     xlab = "Average contact angle (degrees)")

text(curv.food ~ avg.ca, data = wet, labels = lab, font = 2, pos = 4.5, cex = 0.9)

legend("topright",legend="(B)",bty='n',cex=1.2)

curv.food = setNames(wet$curv.food, rownames(wet))

phenogram(pruned2,curv.food,colors=cols, ftype = 'i',
          ylab = "Ratio of radii of curvature (outer/inner)", xlab = "Time since root (my)")

legend("topleft",legend="(C)",bty='n',cex=1.2)

dev.off()
