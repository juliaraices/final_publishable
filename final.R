# January 2016
# Final R program to make:
# - make controls
# - make graphs 1, 2, 3, 4, suplementals et al

total <- read.table("final.output", header=T)
#total <- subset(totaal, totaal$fetValue<=0.05)
total$age <- factor(total$age, levels=c("old","new"))
total$dnds <- as.numeric(total$dnds)
total.a <- subset(total, total$XorA=="A")


summary(subset(total$Meiosis, total$age=="old"))
summary(total$Mitosis)
#max(total$Mitosis)

al <- subset(total, total$dn != "NA")

equal <- subset(al, al$Group=="Equal")

postmeiotic <- subset(al, al$Group=="PostMeiotic")
postmeiotic.a <- subset(postmeiotic, postmeiotic$XorA=="A")
postmeiotic.x <- subset(postmeiotic, postmeiotic$XorA=="X")

meiotic <- subset(al, al$Group=="Meiotic")
meiotic.a <- subset(meiotic, meiotic$XorA=="A")
meiotic.x <- subset(meiotic, meiotic$XorA=="X")

mitotic <- subset(al, al$Group=="Mitotic")
mitotic.a <- subset(mitotic, mitotic$XorA=="A")
mitotic.x <- subset(mitotic, mitotic$XorA=="X")

meiotic.meioticpostmeiotic <- subset(al, al$Group=="Meiotic" | al$Group=="MeioticPostmeiotic")
meiotic.meioticpostmeiotic.a <- subset(meiotic.meioticpostmeiotic, meiotic.meioticpostmeiotic$XorA=="A")
meiotic.meioticpostmeiotic.x <- subset(meiotic.meioticpostmeiotic, meiotic.meioticpostmeiotic$XorA=="X")

meiotic.meioticpostmeiotic.postmeiotic <- subset(al, al$Group=="Meiotic" | al$Group=="MeioticPostmeiotic" | al$Group=="PostMeiotic")
meiotic.meioticpostmeiotic.postmeiotic.a <- subset(meiotic.meioticpostmeiotic.postmeiotic, meiotic.meioticpostmeiotic.postmeiotic$XorA=="A")
meiotic.meioticpostmeiotic.postmeiotic.x <- subset(meiotic.meioticpostmeiotic.postmeiotic, meiotic.meioticpostmeiotic.postmeiotic$XorA=="X")

meioticpostmeiotic.postmeiotic <- subset(al, al$Group=="MeioticPostmeiotic" | al$Group=="PostMeiotic")
meioticpostmeiotic.postmeiotic.a <- subset(meioticpostmeiotic.postmeiotic, meioticpostmeiotic.postmeiotic$XorA=="A")
meioticpostmeiotic.postmeiotic.x <- subset(meioticpostmeiotic.postmeiotic, meioticpostmeiotic.postmeiotic$XorA=="X")

meioticpostmeiotic <- subset(al, al$Group=="MeioticPostmeiotic")
meioticpostmeiotic.a <- subset(meioticpostmeiotic, meioticpostmeiotic$XorA=="A")
meioticpostmeiotic.x <- subset(meioticpostmeiotic, meioticpostmeiotic$XorA=="X")


control.1 <- subset(equal, equal$PostMeiosis >=7.253) #pos
control.1.a <- subset(control.1, control.1$XorA=="A")
control.1.x <- subset(control.1, control.1$XorA=="X")

control.2 <- subset(equal, equal$Meiosis >= 7.423) #mei
control.2.a <- subset(control.2, control.2$XorA=="A")
control.2.x <- subset(control.2, control.2$XorA=="X")

control.3 <- subset(equal, equal$PostMeiosis >= 9.102) #meipos$Pos
control.3.a <- subset(control.3, control.3$XorA=="A")
control.3.x <- subset(control.3, control.3$XorA=="X")

control.4 <- subset(equal, equal$Meiosis >= 8.346) #meipos$Mei
control.4.a <- subset(control.4, control.4$XorA=="A")
control.4.x <- subset(control.4, control.4$XorA=="X")

control.5 <- subset(equal, equal$Meiosis>=7.708) #mei.meipos
control.5.a <- subset(control.5, control.5$XorA=="A")
control.5.x <- subset(control.5, control.5$XorA=="X")

control.6 <- subset(equal, equal$PostMeiosis>=7.402) #meipos.pos
control.6.a <- subset(control.6, control.6$XorA=="A")
control.6.x <- subset(control.6, control.6$XorA=="X")

control.7 <- subset(equal, equal$Meiosis>=6.603) #mei.meipos.pos$Mei
control.7.a <- subset(control.7, control.7$XorA=="A")
control.7.x <- subset(control.7, control.7$XorA=="X")

control.8 <- subset(equal, equal$PostMeiosis>=7.208) #mei.meipos.pos$Pos
control.8.a <- subset(control.8, control.8$XorA=="A")
control.8.x <- subset(control.8, control.8$XorA=="X")

control.mit <- subset(equal, equal$Mitosis>=7.262)
control.mit.a <- subset(control.mit, control.mit$XorA=="A")
control.mit.x <- subset(control.mit, control.mit$XorA=="X")


xall <- length(subset(total$XorA, total$XorA=="X" & total$Group !="Equal"))/length(subset(total$XorA, total$Group !="Equal"))

newxall <- length(subset(total$XorA, total$XorA=="X" & total$age=="new" & total$Group!="Equal"))/length(subset(total$XorA, total$age=="new" & total$Group!="Equal"))

aall <- length(subset(total$XorA, total$XorA=="A" & total$Group!="Equal"))/length(subset(total$XorA, total$Group!="Equal"))

newaall <- length(subset(total$XorA, total$XorA=="A" & total$age=="new" & total$Group!="Equal"))/length(subset(total$XorA, total$age=="new" & total$Group!="Equal"))

datab <- data.frame(Age = c(rep("old",6), rep("new",6)),
                    Chromosome = c(rep(c("A","X"), 6)), 
                    Group = c(rep(c('Mitotic','Meiotic','PostMeiotic'), each=2), rep(c('Mitotic','Meiotic','PostMeiotic'), each=2)),
                    Count = c(length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='X')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='X')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='X')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='X')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='X')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='A')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='X'))),
                    Proportion=c(length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='old')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),
                        length(subset(total$id, total$Group=='Mitotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Mitotic' & total$age=='new')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')),
                        length(subset(total$id, total$Group=='Meiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='Meiotic' & total$age=='new')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='A'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new')),
                        length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new' & total$XorA=='X'))/length(subset(total$id, total$Group=='PostMeiotic' & total$age=='new'))))

datab$Group <- factor(datab$Group, levels=c("Mitotic", "Meiotic", "PostMeiotic"))
datab$Age <- factor(datab$Age, levels=c("old", "new"))
datab$Chromosome <- factor(datab$Chromosome, levels=c("A", "X"))

x.datab <- subset(datab, datab$Chromosome=="X")
a.datab <- subset(datab, datab$Chromosome=="A")

mm <- subset(total, total$Group=='Mitotic' | total$Group=='Meiotic' | total$Group=='PostMeiotic')
mmp <- subset(mm, mm$XorA=='A')

matx <- matrix(data=c(length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='old')),
                      length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='new')),
                      length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='old')),
                      length(subset(mmp$id, mmp$Group=='Meiotic' & mmp$age=='new')),
                      length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='old')),
                      length(subset(mmp$id, mmp$Group=='PostMeiotic' & mmp$age=='new'))), nrow=2)
rownames(matx)<-c("old", "new")
colnames(matx)<-c("Mitotic", "Meiotic", "PostMeiotic")
matp <- prop.table(matx, 1)*100




#################### Figure 1 BW ###################
pdf("figure1_bw.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("dimgray","gray"), xlab="Spermatogenesis phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "Post-Meiotic"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
legend(x=7, y=60, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 

# change the sets so that you can test 2 by 2 the groups per age
#test <- matrix(data=c(length(subset(mmp$id, mmp$Group=='Mitotic' & mmp$age=='old')),
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
## dN/dS ##
boxplot(subset(mitotic.a$dnds, mitotic.a$age=="old"), subset(mitotic.a$dnds, mitotic.a$age=="new"),
    subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"),
    col=c("dimgrey", "grey"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1.5,3.5), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=1, y=-1, inset=0.01, legend=c("new", "old"), fill=c("grey","dimgrey"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(2.6, 2), labels=c(rep("***", 2)), cex=3)
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

boxplot(mitotic.a$dnds, control.mit.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        col=c("Grey 20", "Grey 42"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=0, y=-1, inset=-0.01, legend=c("Experimental", "Control"), fill=c("Grey 20","Grey 42"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(rep(1.5, 2)), labels=c("***", "***"), cex=3)
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

## alpha ##
boxplot(subset(mitotic.a$alpha, mitotic.a$age=="old"), subset(mitotic.a$alpha, mitotic.a$age=="new"),
        subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="old"), subset(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, meiotic.meioticpostmeiotic.postmeiotic.a$age=="new"),
        col=c("dimgrey", "grey"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
#text(x=1.5, y=-1.7, labels="p = 0.3606", cex=1.2)
text(x=3.5, y=-1.7, labels="***", cex=3)
mtext("(c)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

boxplot(mitotic.a$alpha, control.mit.a$alpha,
        meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha,
        col=c("Grey 20", "Grey 42"), ylim=c(-5,1), ylab="alpha", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
#text(x=1.5, y=-2, labels="p = 0.4675", cex=1.2)
text(x=3.5, y=-1.7, labels="", cex=2.5)
mtext("(d)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
dev.off()

######################## Figure 3 BW ####################################
library(ggplot2)

hline.data <- data.frame(z = c(aall, newaall), Age = c("old","new"))
labs.data <- data.frame(s=c(2, 1, 3), f=c(0.878, 0.74, 0.905), z = c("***", "*", "***"), Age = c("old","new", "new"), Chromosome=c("A", "A", "A"))

facet_names <- list(
    'old'="Old Genes",
    'new'="New Genes")

facet_labeller <- function(variable, value){
    return(facet_names[value])
}

pdf("figure3_bw.pdf", width=15, height=10)#, res=300) # width=22, height=10, 
ggplot(datab, aes(x=Group, y=Proportion, fill=Chromosome)) +
    geom_bar(position='stack', stat='identity') +
    facet_grid(.~Age, labeller=facet_labeller) +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme(text = element_text(size=30)) +
    scale_x_discrete(name="") +
    geom_hline(aes(yintercept =z), hline.data) +
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.99,6)), label=c(x.datab$Count), size=7) +
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.61, 6)), label=c(a.datab$Count), size=7) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    #geom_text(x=c(2,1), y=c(0.875, 0.74), label=c("***", "*"), Age=c("old", "new")) + #comparações contra todos os genes que não são equal
    #geom_text(x=c(2, 4, 6)a, y=c(0.886, 0.747, 0.914), label=c("***", "*", "***")) + #comparações contra todos os genes que não são equal
    scale_fill_manual(values=c("Grey 70", "Gray 50"), name="", breaks=c("A", "X"), labels=c("Autosomal gene", "X-linked gene"))
dev.off()

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
legend(x=1, y=-1, inset=0.01, legend=c("new", "old"), fill=c("salmon","powderblue"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(seq(1.5,3.5, by=2)), y=c(2.6, 2), labels=c(rep("***", 2)), cex=2.5)
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)

boxplot(mitotic.a$dnds, control.mit.a$dnds,
        meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds,
        col=c("pink", "plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(seq(1.5,3.5, by=2)), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
legend(x=0, y=-1, inset=-0.01, legend=c("Experimental", "Control"), fill=c("pink","plum"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
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
mtext("(d)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
dev.off()

wilcox.test(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=="new"), subset(control.7$alpha, control.7$age=="new"))

summary(subset(meiotic.meioticpostmeiotic.postmeiotic$alpha, meiotic.meioticpostmeiotic.postmeiotic$age=="new"))
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
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.99,6)), label=c(x.datab$Count), size=7) +
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.61, 6)), label=c(a.datab$Count), size=7) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    #geom_text(x=c(2,1), y=c(0.875, 0.74), label=c("***", "*"), Age=c("old", "new")) + #comparações contra todos os genes que não são equal
    #geom_text(x=c(2, 4, 6)a, y=c(0.886, 0.747, 0.914), label=c("***", "*", "***")) + #comparações contra todos os genes que não são equal
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


################# Expression #######################
boxplot(total.a$Mitosis, )
ggplot(total.a, aes(x=age, y=Mitosis, fill=age)) +
    geom_boxplot()
ggplot(total.a, aes(x=age, y=Meiosis, fill=age)) +
    geom_boxplot()
ggplot(total.a, aes(x=age, y=PostMeiosis, fill=age)) +
    geom_boxplot()


ggplot(datab, aes(x=Group, y=Proportion, fill=Chromosome)) +
    facet_grid(.~Age, labeller=facet_labeller) +
    geom_bar(position='stack', stat='identity') +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme(text = element_text(size=30)) +
    scale_x_discrete(name="") +
    geom_hline(aes(yintercept =z), hline.data) +
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.99,6)), label=c(x.datab$Count), size=7) +
    geom_text(x=c(1, 2, 3, 4, 5, 6), y=c(rep(0.61, 6)), label=c(a.datab$Count), size=7) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    #geom_text(x=c(2,1), y=c(0.875, 0.74), label=c("***", "*"), Age=c("old", "new")) + #comparações contra todos os genes que não são equal
    #geom_text(x=c(2, 4, 6)a, y=c(0.886, 0.747, 0.914), label=c("***", "*", "***")) + #comparações contra todos os genes que não são equal
    scale_fill_manual(values=c("Khaki", "Thistle"), name="", breaks=c("A", "X"), labels=c("Autosomal gene", "X-linked gene"))

