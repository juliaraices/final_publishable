# January 2016
# Final R program to make:
# - make controls
# - make graphs 1, 2, 3, 4, suplementals et al

# reads the output
total <- read.table("output/final.output", header=T)
total$age <- factor(total$age, levels=c("old","new"))
total$dnds <- as.numeric(total$dnds)
# creates a subgroup of autosomal genes
total.a <- subset(total, total$XorA=="A")

# creates a subset for genes with dN data (positive selection signature)
al <- subset(total, total$dn != "NA")

### Groups are related to expression during spermatogenesis
# creates a subset with genes from Equal expression category
equal <- subset(al, al$Group=="Equal")

# subset for PostMeiotic genes (and the sub-subsets for autosomal and X-linked genes)
postmeiotic <- subset(al, al$Group=="PostMeiotic")
postmeiotic.a <- subset(postmeiotic, postmeiotic$XorA=="A")
postmeiotic.x <- subset(postmeiotic, postmeiotic$XorA=="X")

# subset for Meiotic genes (and the sub-subsets for autosomal and X-linked genes)
meiotic <- subset(al, al$Group=="Meiotic")
meiotic.a <- subset(meiotic, meiotic$XorA=="A")
meiotic.x <- subset(meiotic, meiotic$XorA=="X")

# subset for Mitotic genes (and the sub-subsets for autosomal and X-linked genes)
mitotic <- subset(al, al$Group=="Mitotic")
mitotic.a <- subset(mitotic, mitotic$XorA=="A")
mitotic.x <- subset(mitotic, mitotic$XorA=="X")

# subset for Meiotic + MeioticPostMeiotic genes (and the sub-subsets for autosomal and X-linked genes)
meiotic.meioticpostmeiotic <- subset(al, al$Group=="Meiotic" | al$Group=="MeioticPostmeiotic")
meiotic.meioticpostmeiotic.a <- subset(meiotic.meioticpostmeiotic, meiotic.meioticpostmeiotic$XorA=="A")
meiotic.meioticpostmeiotic.x <- subset(meiotic.meioticpostmeiotic, meiotic.meioticpostmeiotic$XorA=="X")

# subset for Meiotic + MeioticPostMeiotic + PostMeiotic genes (and the sub-subsets for autosomal and X-linked genes)
meiotic.meioticpostmeiotic.postmeiotic <- subset(al, al$Group=="Meiotic" | al$Group=="MeioticPostmeiotic" | al$Group=="PostMeiotic")
meiotic.meioticpostmeiotic.postmeiotic.a <- subset(meiotic.meioticpostmeiotic.postmeiotic, meiotic.meioticpostmeiotic.postmeiotic$XorA=="A")
meiotic.meioticpostmeiotic.postmeiotic.x <- subset(meiotic.meioticpostmeiotic.postmeiotic, meiotic.meioticpostmeiotic.postmeiotic$XorA=="X")

# subset for MeioticPostMeiotic + PostMeiotic genes (and the sub-subsets for autosomal and X-linked genes)
meioticpostmeiotic.postmeiotic <- subset(al, al$Group=="MeioticPostmeiotic" | al$Group=="PostMeiotic")
meioticpostmeiotic.postmeiotic.a <- subset(meioticpostmeiotic.postmeiotic, meioticpostmeiotic.postmeiotic$XorA=="A")
meioticpostmeiotic.postmeiotic.x <- subset(meioticpostmeiotic.postmeiotic, meioticpostmeiotic.postmeiotic$XorA=="X")

# subset for MeioticPostMeiotic genes (and the sub-subsets for autosomal and X-linked genes)
meioticpostmeiotic <- subset(al, al$Group=="MeioticPostmeiotic")
meioticpostmeiotic.a <- subset(meioticpostmeiotic, meioticpostmeiotic$XorA=="A")
meioticpostmeiotic.x <- subset(meioticpostmeiotic, meioticpostmeiotic$XorA=="X")

# subset of equal genes to act as control for the PostMeiotic group, as their expression during PostMeiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.1 <- subset(equal, equal$PostMeiosis >=(median(postmeiotic$PostMeiosis)-2)) #pos
if(wilcox.test(postmeiotic$PostMeiosis, control.1$PostMeiosis)$p.value <= 0.05) print("FUCK 1")
control.1.a <- subset(control.1, control.1$XorA=="A")
if(wilcox.test(postmeiotic.a$PostMeiosis, control.1.a$PostMeiosis)$p.value <= 0.05) print("FUCK 1")
control.1.x <- subset(control.1, control.1$XorA=="X")
if(wilcox.test(postmeiotic.x$PostMeiosis, control.1.x$PostMeiosis)$p.value <= 0.05) print("FUCK 1")

# subset of equal genes to act as control for the Meiotic group, as their expression during Meiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.2 <- subset(equal, equal$Meiosis >= (median(meiotic$Meiosis)-2.2)) #mei
if(wilcox.test(meiotic$Meiosis, control.2$Meiosis)$p.value <= 0.05) print("FUCK 2")
control.2.a <- subset(control.2, control.2$XorA=="A")
if(wilcox.test(meiotic.a$Meiosis, control.2.a$Meiosis)$p.value <= 0.05) print("FUCK 2")
control.2.x <- subset(control.2, control.2$XorA=="X")
if(wilcox.test(meiotic.x$Meiosis, control.2.x$Meiosis)$p.value <= 0.05) print("FUCK 2")

# subset of equal genes to act as control for the MeioticPostMeiotic group, as their expression during PostMeiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.3 <- subset(equal, equal$PostMeiosis >= (median(meioticpostmeiotic$PostMeiosis)-2)) #meipos$Pos
if(wilcox.test(meioticpostmeiotic$PostMeiosis, control.3$PostMeiosis)$p.value <= 0.05) print("FUCK 3")
control.3.a <- subset(control.3, control.3$XorA=="A")
if(wilcox.test(meioticpostmeiotic.a$PostMeiosis, control.3.a$PostMeiosis)$p.value <= 0.05) print("FUCK 3")
control.3.x <- subset(control.3, control.3$XorA=="X")
if(wilcox.test(meioticpostmeiotic.x$PostMeiosis, control.3.x$PostMeiosis)$p.value <= 0.05) print("FUCK 3")

# subset of equal genes to act as control for the MeioticPostMeiotic group, as their expression during Meiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.4 <- subset(equal, equal$Meiosis >= (median(meioticpostmeiotic$Meiosis)-3)) #meipos$Mei
if(wilcox.test(meioticpostmeiotic$Meiosis, control.4$Meiosis)$p.value <= 0.05) print("FUCK 4")
control.4.a <- subset(control.4, control.4$XorA=="A")
if(wilcox.test(meioticpostmeiotic.a$Meiosis, control.4.a$Meiosis)$p.value <= 0.05) print("FUCK 4")
control.4.x <- subset(control.4, control.4$XorA=="X")
if(wilcox.test(meioticpostmeiotic.x$Meiosis, control.4.x$Meiosis)$p.value <= 0.05) print("FUCK 4")

# subset of equal genes to act as control for the Meiotic + MeioticPostMeiotic group, as their expression during Meiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.5 <- subset(equal, equal$Meiosis>=(median(meiotic.meioticpostmeiotic$Meiosis)-2.7)) #mei.meipos
if(wilcox.test(meiotic.meioticpostmeiotic$Meiosis, control.5$Meiosis)$p.value <= 0.05) print("FUCK 5")
control.5.a <- subset(control.5, control.5$XorA=="A")
if(wilcox.test(meiotic.meioticpostmeiotic.a$Meiosis, control.5.a$Meiosis)$p.value <= 0.05) print("FUCK 5")
control.5.x <- subset(control.5, control.5$XorA=="X")
if(wilcox.test(meiotic.meioticpostmeiotic.x$Meiosis, control.5.x$Meiosis)$p.value <= 0.05) print("FUCK 5")

# subset of equal genes to act as control for the MeioticPostMeiotic + PostMeiotic group, as their expression during PostMeiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.6 <- subset(equal, equal$PostMeiosis>=(median(meioticpostmeiotic.postmeiotic$PostMeiosis)-2.5)) #meipos.pos
if(wilcox.test(meioticpostmeiotic.postmeiotic$PostMeiosis, control.6$PostMeiosis)$p.value <= 0.05) print("FUCK 6")
control.6.a <- subset(control.6, control.6$XorA=="A")
if(wilcox.test(meioticpostmeiotic.postmeiotic.a$PostMeiosis, control.6.a$PostMeiosis)$p.value <= 0.05) print("FUCK 6")
control.6.x <- subset(control.6, control.6$XorA=="X")
if(wilcox.test(meioticpostmeiotic.postmeiotic.x$PostMeiosis, control.6.x$PostMeiosis)$p.value <= 0.05) print("FUCK 6")

# subset of equal genes to act as control for the Meiotic + MeioticPostMeiotic + PostMeiotic group, as their expression during Meiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.7 <- subset(equal, equal$Meiosis>=(median(meiotic.meioticpostmeiotic.postmeiotic$Meiosis)-2.4)) #mei.meipos.pos$Mei
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$Meiosis, control.7$Meiosis)$p.value <= 0.05) print("FUCK 7")
control.7.a <- subset(control.7, control.7$XorA=="A")
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$Meiosis, control.7.a$Meiosis)$p.value <= 0.05) print("FUCK 7")
control.7.x <- subset(control.7, control.7$XorA=="X")
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$Meiosis, control.7.x$Meiosis)$p.value <= 0.05) print("FUCK 7")

# subset of equal genes to act as control for the Meiotic + MeioticPostMeiotic + PostMeiotic group, as their expression during PostMeiosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.8 <- subset(equal, equal$PostMeiosis>=(median(meiotic.meioticpostmeiotic.postmeiotic$PostMeiosis)-2)) #mei.meipos.pos$Pos
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$PostMeiosis, control.8$PostMeiosis)$p.value <= 0.05) print("FUCK 8")
control.8.a <- subset(control.8, control.8$XorA=="A")
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$PostMeiosis, control.8.a$PostMeiosis)$p.value <= 0.05) print("FUCK 8")
control.8.x <- subset(control.8, control.8$XorA=="X")
if(wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$PostMeiosis, control.8.x$PostMeiosis)$p.value <= 0.05) print("FUCK 8")

# subset of equal genes to act as control for the Mitotic group, as their expression during Mitosis shall be the same (and the sub-subsets for autosomal and X-linked genes)
control.mit <- subset(equal, equal$Mitosis>=(median(mitotic$Mitosis)-2.3))
if(wilcox.test(mitotic$Mitosis, control.mit$Mitosis)$p.value <= 0.05) print("FUCK MIT")
control.mit.a <- subset(control.mit, control.mit$XorA=="A")
if(wilcox.test(mitotic.a$Mitosis, control.mit.a$Mitosis)$p.value <= 0.05) print("FUCK MIT")
control.mit.x <- subset(control.mit, control.mit$XorA=="X")
if(wilcox.test(mitotic.x$Mitosis, control.mit.x$Mitosis)$p.value <= 0.05) print("FUCK MIT")

# gets the proportion of X-linked genes, new X-linked genes, autosomal genes and new autosomal genes in the whole set minus the equal expression category
xall <- length(subset(total$XorA, total$XorA=="X" & total$Group !="Equal"))/length(subset(total$XorA, total$Group !="Equal"))

newxall <- length(subset(total$XorA, total$XorA=="X" & total$age=="new" & total$Group!="Equal"))/length(subset(total$XorA, total$age=="new" & total$Group!="Equal"))

aall <- length(subset(total$XorA, total$XorA=="A" & total$Group!="Equal"))/length(subset(total$XorA, total$Group!="Equal"))

newaall <- length(subset(total$XorA, total$XorA=="A" & total$age=="new" & total$Group!="Equal"))/length(subset(total$XorA, total$age=="new" & total$Group!="Equal"))


# creates a data.frame with the number of genes with a given age (new or old), chromosomal location (x-linked or autosomal), and expression category (mitotic, meiotic or post-meiotic), as well as the proportion of such genes in the data.frame
# datab <- data.frame(Age = c(rep("old",6), rep("new",6)),
#                     Chromosome = c(rep(c("A","X"), 6)), 
#                     Group = c(rep(c('Mitotic','Meiotic','PostMeiotic'), each=2), rep(c('Mitotic','Meiotic','PostMeiotic'), each=2)),
#                     Count = c(length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='X')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='X')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='X')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='X')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='X')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='A')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='X'))),
#                     Proportion=c(length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),
#                         length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')),
#                         length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new')),
#                         length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new'))))
# 
# datab$Group <- factor(datab$Group, levels=c("Mitotic", "Meiotic", "PostMeiotic"))
# datab$Age <- factor(datab$Age, levels=c("old", "new"))
# datab$Chromosome <- factor(datab$Chromosome, levels=c("A", "X"))
# 
# 

# way better way to the what was above:
library(plyr)
datab <- aggregate(total$id, by = list(total$age, total$XorA, total$Group), FUN=length)
datab <- subset(datab, datab$Group.3=="Mitotic" | datab$Group.3=="Meiotic" | datab$Group.3=="PostMeiotic")
datab$Proportion <- c(length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')), #old A mei
        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')), #new A mei
        length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')), #old X mei
        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')), #new X mei
        length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')), #old A mit
        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),#new A mit
        length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')), #old X mit
        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),#new X mit
        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')), #old A pos
        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new')), #new A pos
        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')), #old X pos
        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new'))) #new X pos

datab <- rename(datab, c("Group.1"="Age", "Group.2"="Chromosome", "Group.3"="Group", "x"="Count"))

datab$Group <- factor(datab$Group, levels=c("Mitotic", "Meiotic", "PostMeiotic"))

x.datab <- subset(datab, datab$Chromosome=="X")
a.datab <- subset(datab, datab$Chromosome=="A")

############### How to do the proportions better?

# subset of the total data.frame to create a matrix with ageXexpression category, and be able to create graphs and stats
mmp <- subset(total.a, total.a$Group=='Mitotic' | total.a$Group=='Meiotic' | total.a$Group=='PostMeiotic')

mmp$Group = factor(mmp$Group) ### "joga fora" os factors vazios
matx = table(mmp$age, mmp$Group)
matx <- matx[,c("Mitotic", "Meiotic", "PostMeiotic")]


# matx <- matrix(data=c(length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='old')),
#                       length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='new')),
#                       length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='old')),
#                       length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='new')),
#                       length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='old')),
#                       length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='new'))), nrow=2)
# rownames(matx)<-c("old", "new")
# colnames(matx)<-c("Mitotic", "Meiotic", "PostMeiotic")
matp <- prop.table(matx, 1)*100




#################### Figure 1 BW ###################
pdf("figure1_bw.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
# graph with the proportion of genes in each class, and subsequent stats represented by *
barplot(matp, beside=T, col=c("dimgray","gray"), xlab="Spermatogenesis phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "Post-Meiotic"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
legend(x=7, y=60, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 

# change the sets so that you can test 2 by 2 the groups per age
#test <- matrix(data=c(
#            length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='old')),
#            length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='new')),
#            length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='old')),
#            length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='new')),
#            length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='old')),
#            length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='new'))),
#            nrow=2)
#fisher.test(test)
# mit X mei = < 2.2e-16
# mit X pos = 1.628e-11
# mei X pos = 7.866e-07

#test1 <- matrix(data=c(length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='old')),
#               length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='new')),
#               sum(test[1,])-length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='old')),
#               sum(test[2,])-length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='new'))), nrow=2)
#fisher.test(test1)
# mit X resto = < 2.2e-16
# mei X resto = < 2.2e-16
# pos X resto = 0.01091

###################### Figure 2 BW #################################
pdf("figure2_bw.pdf")#, res=100)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
# boxplot with dn/ds and alpha values
## dN/dS ##
boxplot(subset(mitotic.a$dnds, mitotic.a$age=="old"), subset(mitotic.a$dnds, mitotic.a$age=="new"),
    subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"),
    col=c("dimgrey", "grey"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1.5,3.5), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=1, y=-0.87, inset=0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(2.6, 2), labels=c(rep("***", 2)), cex=3)
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
#wilcox.test(subset(mitotic.a$dnds, mitotic.a$age=="old"), subset(mitotic.a$dnds, mitotic.a$age=="new"))
#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"))


boxplot(mitotic.a$dnds, control.mit.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        col=c("Grey 20", "Grey 42"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=-0.01, legend=c("Experimental", "Control"), fill=c("Grey 20","Grey 42"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(rep(1.5, 2)), labels=c("***", "***"), cex=3)
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
#wilcox.test(mitotic.a$dnds, control.mit.a$dnds)
#wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds)

## alpha ##
boxplot(subset(mitotic.a$alpha, mitotic.a$age=="old"), subset(mitotic.a$alpha, mitotic.a$age=="new"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"),
        col=c("dimgrey", "grey"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
#text(x=1.5, y=-1.7, labels="p = 0.3606", cex=1.2)
text(x=c(3.5), y=c(-1.7), labels=c("***"), cex=3)
text(x=c(1.5), y=c(-2), labels=c("p = 0.289"), cex=1)
mtext("(c)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
#wilcox.test(subset(mitotic.a$alpha, mitotic.a$age=="old"), subset(mitotic.a$alpha, mitotic.a$age=="new"))
#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"))

boxplot(mitotic.a$alpha, control.mit.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha,
        col=c("Grey 20", "Grey 42"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=3.5, y=-1.7, labels="*", cex=2.5)
text(x=c(1.5), y=c(-2), labels=c("p=0.250"), cex=1)
mtext("(d)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
#wilcox.test(mitotic.a$alpha, control.mit.a$alpha)
#wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha)
dev.off()

######################## Figure 3 BW ####################################
library(ggplot2)


x.count <- c(subset(x.datab$Count, x.datab$Age=="old" & x.datab$Group=="Mitotic"), subset(x.datab$Count, x.datab$Age=="old" & x.datab$Group!="Mitotic"), subset(x.datab$Count, x.datab$Age=="new" & x.datab$Group=="Mitotic"), subset(x.datab$Count, x.datab$Age=="new" & x.datab$Group!="Mitotic"))

a.count <- c(subset(a.datab$Count, a.datab$Age=="old" & a.datab$Group=="Mitotic"), subset(a.datab$Count, a.datab$Age=="old" & a.datab$Group!="Mitotic"), subset(a.datab$Count, x.datab$Age=="new" & a.datab$Group=="Mitotic"), subset(a.datab$Count, a.datab$Age=="new" & a.datab$Group!="Mitotic"))

hline.data <- data.frame(z = c(aall, newaall), Age = c("old","new"))
labs.data <- data.frame(s=c(2, 1, 3), f=c(0.878, 0.74, 0.905), z = c("***", "*", "***"), Age = c("old","new", "new"), Chromosome=c("A", "A", "A"))
ps <- data.frame(a=c(1,3,1,2), b=c(0.834, 0.834, 0.74, 0.815), d=c("0.096", "0.493", "", "0.777"), Age=c("old", "old", "new", "new"), Chromosome=c("A", "A", "A", "A"))
texto <- data.frame(x=c(rep(c(1,2,3),4)), y=c(rep(0.99, 6), rep(0.61, 6)), lab=c(x.count, a.count), Age=c(rep(c(rep("old",3), rep("new", 3)),2)))

facet_names <- list(
    'old'="Old Genes",
    'new'="New Genes")

facet_labeller <- function(variable, value){
    return(facet_names[value])
}



# graph with proportion of X-linked and autossomal genes in each expression category for new and old genes.
pdf("figure3_bw.pdf", width=15, height=10)#, res=300) # width=22, height=10, 
ggplot(datab, aes(x=Group, y=Proportion, fill=Chromosome)) +
    geom_bar(position='stack', stat='identity') +
    facet_grid(.~Age, labeller=facet_labeller) +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme(text = element_text(size=30)) +
    scale_x_discrete(name="") +
    geom_hline(aes(yintercept =z), hline.data) +
    geom_text(aes(x=x,y=y,label=lab), texto, size=7, inherit.aes=F) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    geom_text(aes(x=a, y=b, label=d), ps, size=5) +
    scale_fill_manual(values=c("Grey 70", "Gray 50"), name="", breaks=c("A", "X"), labels=c("Autosomal gene", "X-linked gene"))
dev.off()

#x.test <- matrix(data=c(
    #subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic' & datab$Chromosome=='X'))
    
    #subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic' & datab$Chromosome=='X'))

    #subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='old'))-subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic' & datab$Chromosome=='X'))
    
    #subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic' & datab$Chromosome=='X'))
    
    #subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic' & datab$Chromosome=='X'))

    #subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic' & datab$Chromosome=='A'),
    #subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic' & datab$Chromosome=='X'),
    #(sum(subset(datab$Count, datab$Chromosome=='A' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic' & datab$Chromosome=='A')),
    #(sum(subset(datab$Count, datab$Chromosome=='X' & datab$Age=='new'))-subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic' & datab$Chromosome=='X'))
#), nrow=2)


chisq.test(x.test)
fisher.test(x.test)

### Same graphs as before, now in collors.

######################## Figure 1 Color ################################
pdf("figure1_color.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("powderblue","salmon"), xlab="Spermatogenesis phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "Post-Meiotic"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
legend(x=7, y=60, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()


######################### Figure 2 Color ############################
pdf("figure2_color.pdf")#, res=300)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
## dN/dS ##
boxplot(subset(mitotic.a$dnds, mitotic.a$age=="old"), subset(mitotic.a$dnds, mitotic.a$age=="new"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"), 
        col=c("powderblue", "salmon"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1.5,3.5), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=1, y=-0.87, inset=0.01, legend=c("new", "old"), fill=c("salmon","powderblue"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(2.6, 2), labels=c(rep("***", 2)), cex=2.5)
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

boxplot(mitotic.a$dnds, control.mit.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        col=c("pink", "plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=-0.01, legend=c("Experimental", "Control"), fill=c("pink","plum"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(rep(1.5, 2)), labels=c("***", "***"), cex=2.5)
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
## alpha ##
boxplot(subset(mitotic.a$alpha, mitotic.a$age=="old"), subset(mitotic.a$alpha, mitotic.a$age=="new"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"),
        col=c("powderblue", "salmon"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=3.5, y=-1.7, labels="***", cex=2.5)
mtext("(c)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

boxplot(mitotic.a$alpha, control.mit.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha,
        col=c("pink", "plum"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=3.5, y=-1.7, labels="*", cex=2.5)
text(x=c(1.5), y=c(-2), labels=c("p=0.250"), cex=1)
mtext("(d)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
#wilcox.test(mitotic.a$alpha, control.mit.a$alpha)
#wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha)
dev.off()

#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=="new"), subset(control.7$alpha, control.7$age=="new"))

#summary(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=="new"))
summary(subset(control.7$alpha, control.7$age=="new"))

####################### Figure 3 Color #########################
pdf("figure3_color.pdf", width=15, height=10)#, res=300) # width=22, height=10, 
ggplot(datab, aes(x=Group, y=Proportion, fill=Chromosome)) +
    facet_grid(.~Age, labeller=facet_labeller) +
    geom_bar(position='stack', stat='identity') +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme(text = element_text(size=30)) +
    scale_x_discrete(name="") +
    geom_hline(aes(yintercept =z), hline.data) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    geom_text(aes(x=x,y=y,label=lab), texto, size=7, inherit.aes=F) +
    geom_text(aes(x=a, y=b, label=d), ps, size=5) +
    scale_fill_manual(values=c("Khaki", "Thistle"), name="", breaks=c("A", "X"), labels=c("Autosomal gene", "X-linked gene"))
dev.off()


########################### SUPPLEMENTAL FIGURES ##################################
############## BW ###############
pdf("sup1_bw.pdf", width=30, height=20)#, res=300)#)
par(mfrow=c(3,1), mar=c(8,6,4,2)+0.1)
## dN/dS ##
boxplot(subset(postmeiotic.a$dnds, postmeiotic.a$age=="new"), subset(postmeiotic.a$dnds, postmeiotic.a$age=="old"),
        subset(meiotic.a$dnds, meiotic.a$age=="new"), subset(meiotic.a$dnds, meiotic.a$age=="old"),
        subset(meioticpostmeiotic.a$dnds, meioticpostmeiotic.a$age=="new"), subset(meioticpostmeiotic.a$dnds, meioticpostmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.a$dnds, meiotic.meioticpostmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.a$dnds, meiotic.meioticpostmeiotic.a$age=="old"),
        subset(meioticpostmeiotic.postmeiotic.a$dnds, meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meioticpostmeiotic.postmeiotic.a$dnds, meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(mitotic.a$dnds, mitotic.a$age=="new"), subset(mitotic.a$dnds, mitotic.a$age=="old"),
        col=c("grey", "dimgrey"), ylim=c(0,4), ylab="dN/dS", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for New and Old genes from Experimental Groups", cex.main=3)
axis(1, at=c(seq(1.5,13.5, by=2)), labels=c("Experimental I", "Experimental II", "Experimental III/IV", "Experimental V", "Experimental VI", "Experimental VII/VIII", "Experimental IX"), cex.axis=2)
legend(x=13, y=4.6, inset=-0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), cex=2, horiz = T, xpd = T, bty = "n")
text(x=c(seq(1.5,13.5, by=2)), y=c(rep(2, 7)), labels=c(rep("***", 7)), cex=3)

boxplot(subset(control.1.a$dnds, control.1.a$age=="new"), subset(control.1.a$dnds, control.1.a$age=="old"),
        subset(control.2.a$dnds, control.2.a$age=="new"), subset(control.2.a$dnds, control.2.a$age=="old"),
        subset(control.3.a$dnds, control.3.a$age=="new"), subset(control.3.a$dnds, control.3.a$age=="old"),
        subset(control.4.a$dnds, control.4.a$age=="new"), subset(control.4.a$dnds, control.4.a$age=="old"),
        subset(control.5.a$dnds, control.5.a$age=="new"), subset(control.5.a$dnds, control.5.a$age=="old"),
        subset(control.6.a$dnds, control.6.a$age=="new"), subset(control.6.a$dnds, control.6.a$age=="old"),
        subset(control.7.a$dnds, control.7.a$age=="new"), subset(control.7.a$dnds, control.7.a$age=="old"),
        subset(control.8.a$dnds, control.8.a$age=="new"), subset(control.8.a$dnds, control.8.a$age=="old"),
        subset(control.mit.a$dnds, control.mit.a$age=="new"), subset(control.mit.a$dnds, control.mit.a$age=="old"),
        col=c("grey", "dimgrey"), ylim=c(0,4), ylab="dN/dS", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for New and Old genes from Control Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("Control I", "Control II", "Control III", "Control IV", "Control V", "Control VI", "Control VII", "Control VIII", "Control IX"), cex.axis=2)
#legend("topright", inset=0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), title= "Age", cex=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(2.7, 9)), labels=c("","**","", "", "", "", "*", "", "**"), cex=3)

boxplot(postmeiotic.a$dnds, control.1.a$dnds,
        meiotic.a$dnds, control.2.a$dnds,
        meioticpostmeiotic.a$dnds, control.3.a$dnds,
        meioticpostmeiotic.a$dnds, control.4.a$dnds,
        meiotic.meioticpostmeiotic.a$dnds, control.5.a$dnds,
        meioticpostmeiotic.postmeiotic.a$dnds, control.6.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.7.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        mitotic.a$dnds, control.mit.a$dnds,
        col=c("Grey 20", "Grey 42"), ylim=c(0,4), ylab="dN/dS", xlab="Group Number", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for Experimental and Control genes for each Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"), cex.axis=2)
legend(y=4.7, x=16, inset=-0.01, legend=c("Experimental", "Control"), fill=c("Grey 20","Grey 42"), cex=2, xpd = T, bty ="n", horiz = T)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(1.5, 9)), labels=c("","***","*", "**", "**", "", "***", "***", "**"), cex=3)

## alpha ##
boxplot(subset(postmeiotic.a$alpha, postmeiotic.a$age=="new"), subset(postmeiotic.a$alpha, postmeiotic.a$age=="old"),
        subset(meiotic.a$alpha, meiotic.a$age=="new"), subset(meiotic.a$alpha, meiotic.a$age=="old"),
        subset(meioticpostmeiotic.a$alpha, meioticpostmeiotic.a$age=="new"), subset(meioticpostmeiotic.a$alpha, meioticpostmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.a$alpha, meiotic.meioticpostmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.a$alpha, meiotic.meioticpostmeiotic.a$age=="old"),
        subset(meioticpostmeiotic.postmeiotic.a$alpha, meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meioticpostmeiotic.postmeiotic.a$alpha, meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(mitotic.a$alpha, mitotic.a$age=="new"), subset(mitotic.a$alpha, mitotic.a$age=="old"),
        col=c("grey", "dimgrey"), ylim=c(-6,1), ylab="alpha", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for New and Old genes from Experimental Groups", cex.main=3)
axis(1, at=c(seq(1.5,13.5, by=2)), labels=c("Experimental I", "Experimental II", "Experimental III/IV", "Experimental V", "Experimental VI", "Experimental VII/VIII", "Experimental IX"), cex.axis=2)
#legend("topright", inset=0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), title= "Age", cex=2)
text(x=c(seq(1.5,13.5, by=2)), y=c(rep(0.9, 7)), labels=c("**","","*", "", "", "", "**", "***", ""), cex=3)

boxplot(subset(control.1.a$alpha, control.1.a$age=="new"), subset(control.1.a$alpha, control.1.a$age=="old"),
        subset(control.2.a$alpha, control.2.a$age=="new"), subset(control.2.a$alpha, control.2.a$age=="old"),
        subset(control.3.a$alpha, control.3.a$age=="new"), subset(control.3.a$alpha, control.3.a$age=="old"),
        subset(control.4.a$alpha, control.4.a$age=="new"), subset(control.4.a$alpha, control.4.a$age=="old"),
        subset(control.5.a$alpha, control.5.a$age=="new"), subset(control.5.a$alpha, control.5.a$age=="old"),
        subset(control.6.a$alpha, control.6.a$age=="new"), subset(control.6.a$alpha, control.6.a$age=="old"),
        subset(control.7.a$alpha, control.7.a$age=="new"), subset(control.7.a$alpha, control.7.a$age=="old"),
        subset(control.8.a$alpha, control.8.a$age=="new"), subset(control.8.a$alpha, control.8.a$age=="old"),
        subset(control.mit.a$alpha, control.mit.a$age=="new"), subset(control.mit.a$alpha, control.mit.a$age=="old"),
        col=c("grey", "dimgrey"), ylim=c(-6,1), ylab="alpha", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for New and Old genes from Control Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("Control I", "Control II", "Control III", "Control IV", "Control V", "Control VI", "Control VII", "Control VIII", "Control IX"), cex.axis=2)
#legend("topright", inset=0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), title= "Age", cex=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(0.9, 9)), labels=c("","","", "", "", "", "*", "", ""), cex=3)

boxplot(postmeiotic.a$alpha, control.1.a$alpha,
        meiotic.a$alpha, control.2.a$alpha,
        meioticpostmeiotic.a$alpha, control.3.a$alpha,
        meioticpostmeiotic.a$alpha, control.4.a$alpha,
        meiotic.meioticpostmeiotic.a$alpha, control.5.a$alpha,
        meioticpostmeiotic.postmeiotic.a$alpha, control.6.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.7.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha,
        mitotic.a$alpha, control.mit.a$alpha,
        col=c("Grey 20", "Grey 42"), ylim=c(-6,1), ylab="alpha", xlab="Group Number", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for Experimental and Control genes for each Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"), cex.axis=3)
#legend("topright", inset=0.007, legend=c("Experimental", "Control"), fill=c("Grey 20","Grey 42"), cex=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(0.9, 9)), labels=c("","","", "", "", "", "", "", ""), cex=3)
dev.off()


################ Color ##################
pdf("sup1_color.pdf", width=30, height=20)#, res=300)#, width=30, height=20)
par(mfrow=c(3,1), mar=c(8,6,4,2)+0.1)
## dN/dS ##
boxplot(subset(postmeiotic.a$dnds, postmeiotic.a$age=="new"), subset(postmeiotic.a$dnds, postmeiotic.a$age=="old"),
        subset(meiotic.a$dnds, meiotic.a$age=="new"), subset(meiotic.a$dnds, meiotic.a$age=="old"),
        subset(meioticpostmeiotic.a$dnds, meioticpostmeiotic.a$age=="new"), subset(meioticpostmeiotic.a$dnds, meioticpostmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.a$dnds, meiotic.meioticpostmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.a$dnds, meiotic.meioticpostmeiotic.a$age=="old"),
        subset(meioticpostmeiotic.postmeiotic.a$dnds, meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meioticpostmeiotic.postmeiotic.a$dnds, meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(mitotic.a$dnds, mitotic.a$age=="new"), subset(mitotic.a$dnds, mitotic.a$age=="old"),
        col=c("salmon", "powderblue"), ylim=c(0,4), ylab="dN/dS", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for New and Old genes from Experimental Groups", cex.main=3)
axis(1, at=c(seq(1.5,13.5, by=2)), labels=c("Experimental I", "Experimental II", "Experimental III/IV", "Experimental V", "Experimental VI", "Experimental VII/VIII", "Experimental IX"), cex.axis=2)
legend(x=13, y=4.6, inset=-0.01, legend=c("new", "old"), fill=c("salmon","powderblue"), cex=2, horiz = T, xpd=T, bty="n")
text(x=c(seq(1.5,13.5, by=2)), y=c(rep(2, 7)), labels=c(rep("***", 7)), cex=3)

boxplot(subset(control.1.a$dnds, control.1.a$age=="new"), subset(control.1.a$dnds, control.1.a$age=="old"),
        subset(control.2.a$dnds, control.2.a$age=="new"), subset(control.2.a$dnds, control.2.a$age=="old"),
        subset(control.3.a$dnds, control.3.a$age=="new"), subset(control.3.a$dnds, control.3.a$age=="old"),
        subset(control.4.a$dnds, control.4.a$age=="new"), subset(control.4.a$dnds, control.4.a$age=="old"),
        subset(control.5.a$dnds, control.5.a$age=="new"), subset(control.5.a$dnds, control.5.a$age=="old"),
        subset(control.6.a$dnds, control.6.a$age=="new"), subset(control.6.a$dnds, control.6.a$age=="old"),
        subset(control.7.a$dnds, control.7.a$age=="new"), subset(control.7.a$dnds, control.7.a$age=="old"),
        subset(control.8.a$dnds, control.8.a$age=="new"), subset(control.8.a$dnds, control.8.a$age=="old"),
        subset(control.mit.a$dnds, control.mit.a$age=="new"), subset(control.mit.a$dnds, control.mit.a$age=="old"),
        col=c("salmon", "powderblue"), ylim=c(0,4), ylab="dN/dS", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for New and Old genes from Control Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("Control I", "Control II", "Control III", "Control IV", "Control V", "Control VI", "Control VII", "Control VIII", "Control IX"), cex.axis=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(2.7, 9)), labels=c("","**","", "", "", "", "*", "", "**"), cex=3)

boxplot(postmeiotic.a$dnds, control.1.a$dnds,
        meiotic.a$dnds, control.2.a$dnds,
        meioticpostmeiotic.a$dnds, control.3.a$dnds,
        meioticpostmeiotic.a$dnds, control.4.a$dnds,
        meiotic.meioticpostmeiotic.a$dnds, control.5.a$dnds,
        meioticpostmeiotic.postmeiotic.a$dnds, control.6.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.7.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        mitotic.a$dnds, control.mit.a$dnds,
        col=c("pink", "plum"), ylim=c(0,4), ylab="dN/dS", xlab="Group Number", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("dN/dS for Experimental and Control genes for each Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"), cex.axis=2)
legend(x=16, y=4.6, inset=-0.01, legend=c("Experimental", "Control"), fill=c("pink","plum"), cex=2, horiz = T, xpd=T, bty="n")
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(1.5, 9)), labels=c("","***","*", "**", "**", "", "***", "***", "**"), cex=3)

## alpha ##
boxplot(subset(postmeiotic.a$alpha, postmeiotic.a$age=="new"), subset(postmeiotic.a$alpha, postmeiotic.a$age=="old"),
        subset(meiotic.a$alpha, meiotic.a$age=="new"), subset(meiotic.a$alpha, meiotic.a$age=="old"),
        subset(meioticpostmeiotic.a$alpha, meioticpostmeiotic.a$age=="new"), subset(meioticpostmeiotic.a$alpha, meioticpostmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.a$alpha, meiotic.meioticpostmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.a$alpha, meiotic.meioticpostmeiotic.a$age=="old"),
        subset(meioticpostmeiotic.postmeiotic.a$alpha, meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meioticpostmeiotic.postmeiotic.a$alpha, meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"),
        subset(mitotic.a$alpha, mitotic.a$age=="new"), subset(mitotic.a$alpha, mitotic.a$age=="old"),
        col=c("salmon", "powderblue"), ylim=c(-6,1), ylab="alpha", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for New and Old genes from Experimental Groups", cex.main=3)
axis(1, at=c(seq(1.5,13.5, by=2)), labels=c("Experimental I", "Experimental II", "Experimental III/IV", "Experimental V", "Experimental VI", "Experimental VII/VIII", "Experimental IX"), cex.axis=2)
text(x=c(seq(1.5,13.5, by=2)), y=c(rep(0.9, 7)), labels=c("**","","*", "", "", "", "**", "***", ""), cex=3)

boxplot(subset(control.1.a$alpha, control.1.a$age=="new"), subset(control.1.a$alpha, control.1.a$age=="old"),
        subset(control.2.a$alpha, control.2.a$age=="new"), subset(control.2.a$alpha, control.2.a$age=="old"),
        subset(control.3.a$alpha, control.3.a$age=="new"), subset(control.3.a$alpha, control.3.a$age=="old"),
        subset(control.4.a$alpha, control.4.a$age=="new"), subset(control.4.a$alpha, control.4.a$age=="old"),
        subset(control.5.a$alpha, control.5.a$age=="new"), subset(control.5.a$alpha, control.5.a$age=="old"),
        subset(control.6.a$alpha, control.6.a$age=="new"), subset(control.6.a$alpha, control.6.a$age=="old"),
        subset(control.7.a$alpha, control.7.a$age=="new"), subset(control.7.a$alpha, control.7.a$age=="old"),
        subset(control.8.a$alpha, control.8.a$age=="new"), subset(control.8.a$alpha, control.8.a$age=="old"),
        subset(control.mit.a$alpha, control.mit.a$age=="new"), subset(control.mit.a$alpha, control.mit.a$age=="old"),
        col=c("salmon", "powderblue"), ylim=c(-6,1), ylab="alpha", xlab="Group", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for New and Old genes from Control Groups", cex=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("Control I", "Control II", "Control III", "Control IV", "Control V", "Control VI", "Control VII", "Control VIII", "Control IX"), cex.axis=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(0.9, 9)), labels=c("","","", "", "", "", "*", "", ""), cex=3)

boxplot(postmeiotic.a$alpha, control.1.a$alpha,
        meiotic.a$alpha, control.2.a$alpha,
        meioticpostmeiotic.a$alpha, control.3.a$alpha,
        meioticpostmeiotic.a$alpha, control.4.a$alpha,
        meiotic.meioticpostmeiotic.a$alpha, control.5.a$alpha,
        meioticpostmeiotic.postmeiotic.a$alpha, control.6.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.7.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha,
        mitotic.a$alpha, control.mit.a$alpha,
        col=c("pink", "plum"), ylim=c(-6,1), ylab="alpha", xlab="Group Number", xaxt="n", outline=F, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
title("Alpha for Experimental and Control genes for each Groups", cex.main=3)
axis(1, at=c(seq(1.5,17.5, by=2)), labels=c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"), cex.axis=2)
text(x=c(seq(1.5,17.5, by=2)), y=c(rep(0.9, 9)), labels=c("","","", "", "", "", "", "", ""), cex=3)
dev.off()

dev.off()

######################## Last analysis ##############################
# wilcox new genes X vs new genes A (mit, mei, pm)

wilcox.test(subset(total$Mitosis, total$XorA=='X'), subset(total$Mitosis, total$XorA=='A'))
ks.test(subset(total$Mitosis, total$XorA=='X'), subset(total$Mitosis, total$XorA=='A'))

wilcox.test(subset(total$Meiosis, total$XorA=='X'), subset(total$Meiosis, total$XorA=='A'))
ks.test(subset(total$Meiosis, total$XorA=='X'), subset(total$Meiosis, total$XorA=='A'))

wilcox.test(subset(total$PostMeiosis, total$XorA=='X'), subset(total$PostMeiosis, total$XorA=='A'))
ks.test(subset(total$PostMeiosis, total$XorA=='X'), subset(total$PostMeiosis, total$XorA=='A'))

# wilcox alpha new genes only (hap vs cont, mit vs cont)
#wilcox.test(subset(mitotic.a$alpha, mitotic.a$age=='new'), subset(control.mit.a$alpha, control.mit.a$age=='new'))
#ks.test(subset(mitotic.a$alpha, mitotic.a$age=='new'), subset(control.mit.a$alpha, control.mit.a$age=='new'))

#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=='new'), subset(control.8$alpha, control.8$age=='new'))
#ks.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=='new'), subset(control.8$alpha, control.8$age=='new'))
#de fato o n é muito baixo pra qlq coisa....


# wilcox alpha (hap vs mit)
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$alpha, mitotic$alpha)
ks.test(meiotic.meioticpostmeiotic.postmeiotic$alpha, mitotic$alpha)

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, mitotic.a$alpha)
ks.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, mitotic.a$alpha)

# wilcox alpha new genes only (hap vs mit) nada =/
#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=='new'), subset(mitotic$alpha, mitotic$age=='new'))
#ks.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=='new'), subset(mitotic$alpha, mitotic$age=='new'))

#wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=='new'), subset(mitotic.a$alpha, mitotic.a$age=='new'))
#ks.test(subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=='new'), subset(mitotic.a$alpha, mitotic.a$age=='new'))

