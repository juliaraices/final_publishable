# January 2016
# Final R program to make:
# - make controls
# - make graphs 1, 2, 3, 4, suplementals et al

divas <- read.table("input/melsubgroup_analysis_results_flydivas_v1.2", header=T, sep = "\t")
# reads the output
total <- read.table("output/final.output", header=T)
total$age <- factor(total$age, levels=c("old","new"))
total$dnds <- as.numeric(total$dnds)
# cross my output with flydivas output using cg column from flydivas and id from my output 
mydivas <- merge(total, divas, by.x="id", by.y = "cg")
a.mydivas <- subset(mydivas, mydivas$XorA=="A") # only autosomal genes
# creates a subgroup of autosomal genes
total.a <- subset(total, total$XorA=="A")
# creates a subset for genes with dN data (positive selection signature)
al <- subset(total, total$dn != "NA")

#################### Figure 1 ###################
# graph with the proportion of genes in each class, and subsequent stats represented by *
# dados a serem usados: todos os genes dos grupos mitotico, meiotico e posmeiotoicos apenas
mmp <- subset(total.a, total.a$Group=='Mitotic' | total.a$Group=='Meiotic' | total.a$Group=='PostMeiotic')
mmp$Group = factor(mmp$Group) ### "joga fora" os factors vazios
matx = table(mmp$age, mmp$Group)
matx <- matx[,c("Mitotic", "Meiotic", "PostMeiotic")]
matp <- prop.table(matx, 1)*100
# estatísticas
mtmmtx <- matx[,c(1,2)]
mtpmtx <- matx[,c(1,3)]
mpmtx <- matx[,c(2,3)]
chisq.test(matx)
chisq.test(mtmmtx)
chisq.test(mtpmtx)
chisq.test(mpmtx)
fisher.test(matx)
fisher.test(mtmmtx)
fisher.test(mtpmtx)
fisher.test(mpmtx)
# qual é mais apropriado? o chi-quadrado ou o fisher? ¯\_(ツ)_/¯
########## BW ##########
pdf("figure1_bw.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("dimgray","gray"), xlab="Spermatogenesis phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "Post-Meiotic"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
text(x=c(1.5,2.5,4.5,5.5,7.5,8.5), y=c(rep(11,6)), labels=c(matx[1,1], matx[2,1], matx[1,2], matx[2,2], matx[1,3], matx[2,3]), cex=1.5)
legend(x=7, y=60, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 
pdf("figure1_bw_pt.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("dimgray","gray"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "Pós-Meiótica"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
text(x=c(1.5,2.5,4.5,5.5,7.5,8.5), y=c(rep(11,6)), labels=c(matx[1,1], matx[2,1], matx[1,2], matx[2,2], matx[1,3], matx[2,3]), cex=1.5)
legend(x=6, y=60, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 
########## Color ##########
pdf("figure1_color.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("powderblue","salmon"), xlab="Spermatogenesis phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "Post-Meiotic"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
text(x=c(1.5,2.5,4.5,5.5,7.5,8.5), y=c(rep(11,6)), labels=c(matx[1,1], matx[2,1], matx[1,2], matx[2,2], matx[1,3], matx[2,3]), cex=1.5)
legend(x=7, y=60, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()
pdf("figure1_color_pt.pdf")#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(matp, beside=T, col=c("powderblue","salmon"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(10,60), xpd=F, xaxt="n")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "Pós-Meiótica"), cex.axis=1.5, lwd=0)
text(x=c(3.6,6.6,5), y=c(30,32,52), labels=c("***", "***", "***"), cex=3)
segments(2.8,28,4.2,28, cex=2, lwd=4)
segments(5.8,30,7.2,30, cex=2, lwd=4)
segments(2.8,50,7.2,50, cex=2,lwd=4)
text(x=c(1.5,2.5,4.5,5.5,7.5,8.5), y=c(rep(11,6)), labels=c(matx[1,1], matx[2,1], matx[1,2], matx[2,2], matx[1,3], matx[2,3]), cex=1.5)
legend(x=6, y=60, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()

###################### Figure 2 #################################
# boxplot with dn/ds and alpha values
# dados a serem utilizados: genes autossomicos com valores de dn/ds e alpha, dos grupos mitótico e haplóide (Meiotic + MeioticPostMeiotic + Postmeiotic)
int <- subset(al, al$XorA=="A" & al$Group!="Impossible" & al$Group!="Equal" & al$Group!="TheV")
mit <- subset(int, int$Group=='Mitotic')
hap <- subset(int, int$Group=='Meiotic' | int$Group=='MeioticPostmeiotic' | int$Group=='PostMeiotic')
wilcox.test(subset(int$dnds, int$age=='old'), subset(int$dnds, int$age=='new'))
wilcox.test(subset(mit$dnds, mit$age=='old'), subset(hap$dnds, hap$age=='old'))
wilcox.test(mit$dnds, hap$dnds)
wilcox.test(subset(mit$alpha, mit$age=='old'), subset(hap$alpha, hap$age=='old'))
#flydivas
interest.a <- subset(a.mydivas, a.mydivas$Group!="TheV" | a.mydivas$Group!="Impossible" | a.mydivas!="Equal")
interest.a.old <- subset(interest.a, interest.a$age=='old')
interest.a.new <- subset(interest.a, interest.a$age=='new')
new.old.a.interest <- matrix(data=c(length(subset(interest.a.new$id, interest.a.new$pos78=='does' | interest.a.new$pos_12=='does' | interest.a.new$pos_88a=='does')), (length(interest.a.new$id) - length(subset(interest.a.new$id, interest.a.new$pos78=='does' | interest.a.new$pos_12=='does' | interest.a.new$pos_88a=='does'))), length(subset(interest.a.old$id, interest.a.old$pos78=='does' | interest.a.old$pos_12=='does' | interest.a.old$pos_88a=='does')), (length(interest.a.old$id) - length(subset(interest.a.old$id, interest.a.old$pos78=='does' | interest.a.old$pos_12=='does' | interest.a.old$pos_88a=='does')))), ncol=2)
new.old.a.interest.p <- prop.table(new.old.a.interest, 2)*100
fisher.test(new.old.a.interest) # between new and old genes from interest groups (all but equal, impossible e theV), there is no difference in the number of genes with positive selection markers (p = 0.2635)
# testing between haploid group and mitotic group genes
hap.a <- subset(a.mydivas, a.mydivas$HaploidGroup=='haploid_group')
mit.a <- subset(a.mydivas, a.mydivas$Group=='Mitotic')
mit.hap.a.test <- matrix(data=c(length(subset(hap.a$id, hap.a$pos78=='does' | hap.a$pos_12=='does' | hap.a$pos_88a=='does')), (length(hap.a$id) - length(subset(hap.a$id, hap.a$pos78=='does' | hap.a$pos_12=='does' | hap.a$pos_88a=='does'))), length(subset(mit.a$id, mit.a$pos78=='does' | mit.a$pos_12=='does' | mit.a$pos_88a=='does')), (length(mit.a$id) - length(subset(mit.a$id, mit.a$pos78=='does' | mit.a$pos_12=='does' | mit.a$pos_88a=='does')))), ncol=2)
mit.hap.a.test.p <- prop.table(mit.hap.a.test, 2)*100
fisher.test(mit.hap.a.test) # there are more autosomal genes with signatures of poositive selection in the haploid group then in the mitotic group (p = 2.795e-06)
bab <- apply(new.old.a.interest.p, 1, rev)
bab <- apply(bab, 1, rev)
beb <- apply(mit.hap.a.test.p, 1, rev)
beb <- apply(beb, 1, rev)
pdf("figure2_bw.pdf")#, res=100)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
## dN/dS ##
boxplot(subset(int$dnds, int$age=="old"), subset(int$dnds, int$age=="new"),
        col=c("ivory3", "ivory3"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Old", "New"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=0.01, legend=c("w/o selection", "w/ selection"), fill=c("dimgrey","grey"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(1.5), y=c(2), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(subset(int$dnds, int$age=='old'))$stats[5]+0.1), (boxplot.stats(subset(int$dnds, int$age=='new'))$stats[5] + 0.1)), labels=c(length(subset(int$dnds, int$age=="old")), length(subset(int$dnds, int$age=="new"))))
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
boxplot(mit$dnds, hap$dnds, col=c("ivory3"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=c(1.5), y=c(1.5), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(mit$dnds)$stats[5]+0.1), (boxplot.stats(hap$dnds)$stats[5]+0.1)), labels = c(length(mit$dnds), length(hap$dnds)))
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
### % ###
barplot(bab, col = c("dimgrey","grey"), xpd=F, xlab = "", ylab="Genes Percentage", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Old", "New"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(new.old.a.interest[2,2], new.old.a.interest[1,2], new.old.a.interest[2,1], new.old.a.interest[1,1]), cex=1)
mtext("(c)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
barplot(beb, col = c("dimgrey","grey"), xpd=F, xlab = "", ylab="Genes Percentage", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(mit.hap.a.test[2,2], mit.hap.a.test[1,2], mit.hap.a.test[2,1], mit.hap.a.test[1,1]), cex=1)
mtext("(d)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
dev.off()
pdf("figure2_bw_pt.pdf")#, res=100)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
## dN/dS ##
boxplot(subset(int$dnds, int$age=="old"), subset(int$dnds, int$age=="new"),
        col=c("ivory3", "ivory3"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Antigos", "Novos"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=0.01, legend=c("sem seleção", "com seleção"), fill=c("dimgrey","grey"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(1.5), y=c(2), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(subset(int$dnds, int$age=='old'))$stats[5]+0.1), (boxplot.stats(subset(int$dnds, int$age=='new'))$stats[5] + 0.1)), labels=c(length(subset(int$dnds, int$age=="old")), length(subset(int$dnds, int$age=="new"))))
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
boxplot(mit$dnds, hap$dnds,
        col=c("ivory3"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Mitótico", "Haplóide"), cex.axis=1.5, lwd=0)
text(x=c(1.5), y=c(1.5), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(mit$dnds)$stats[5]+0.1), (boxplot.stats(hap$dnds)$stats[5]+0.1)), labels = c(length(mit$dnds), length(hap$dnds)))
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
## % ##
barplot(bab, col = c("dimgrey","grey"), xpd=F, xlab = "", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Antigos", "Novos"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(new.old.a.interest[2,2], new.old.a.interest[1,2], new.old.a.interest[2,1], new.old.a.interest[1,1]), cex=1)
legend(x=6, y=60, inset=c(-5,-0.5),legend=c("com seleção", "sem seleção"), fill=c("dimgrey","grey"), bty="n", cex=1.5, xpd = T)
mtext("(c)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
barplot(beb, col = c("dimgrey","grey"), xpd=F, xlab = "", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Mitótico", "Haplóide"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(mit.hap.a.test[2,2], mit.hap.a.test[1,2], mit.hap.a.test[2,1], mit.hap.a.test[1,1]), cex=1)
mtext("(d)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
dev.off()
########## Color ##########
pdf("figure2_color.pdf")#, res=300)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
## dN/dS ##
boxplot(subset(int$dnds, int$age=="old"), subset(int$dnds, int$age=="new"),
        col=c("plum", "plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1.5,3.5), labels=c("Old", "New"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=0.01, legend=c("w/o selection", "w/ selection"), fill=c("orchid4","orchid1"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(1.5), y=c(2), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(subset(int$dnds, int$age=='old'))$stats[5]+0.1), (boxplot.stats(subset(int$dnds, int$age=='new'))$stats[5] + 0.1)), labels=c(length(subset(int$dnds, int$age=="old")), length(subset(int$dnds, int$age=="new"))))
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
boxplot(mit$dnds, hap$dnds,
        col=c("plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=c(1.5), y=c(1.5), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(mit$dnds)$stats[5]+0.1), (boxplot.stats(hap$dnds)$stats[5]+0.1)), labels = c(length(mit$dnds), length(hap$dnds)))
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
## % ##
barplot(bab, col = c("orchid4","orchid1"), xpd=F, xlab = "", ylab="Genes Percentage", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Old", "New"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(new.old.a.interest[2,2], new.old.a.interest[1,2], new.old.a.interest[2,1], new.old.a.interest[1,1]), cex=1)
legend(x=6, y=60, inset=c(-5,-0.5),legend=c("w/ selection", "w/o selection"), fill=c("orchid4","orchid1"), bty="n", cex=1.5, xpd = T)
mtext("(c)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
barplot(beb, col = c("orchid4","orchid1"), xpd=F, xlab = "", ylab="Genes Percentage", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Mitotic", "Haploid"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(mit.hap.a.test[2,2], mit.hap.a.test[1,2], mit.hap.a.test[2,1], mit.hap.a.test[1,1]), cex=1)
mtext("(d)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
dev.off()
pdf("figure2_color_pt.pdf")#, res=300)#, width=10, height=10)
par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
## dN/dS ##
boxplot(subset(int$dnds, int$age=="old"), subset(int$dnds, int$age=="new"),
        col=c("plum", "plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Antigos", "Novos"), cex.axis=1.5, lwd=0)
legend(x=0, y=-0.87, inset=0.01, legend=c("com seleção", "sem seleção"), fill=c("orchid1","orchid4"), horiz=TRUE, cex=1.3, bty="n", xpd=TRUE)
text(x=c(1.5), y=c(2), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(subset(int$dnds, int$age=='old'))$stats[5]+0.1), (boxplot.stats(subset(int$dnds, int$age=='new'))$stats[5] + 0.1)), labels=c(length(subset(int$dnds, int$age=="old")), length(subset(int$dnds, int$age=="new"))))
mtext("(a)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
boxplot(mit$dnds, hap$dnds,
        col=c("plum"), ylim=c(0,4), ylab="dN/dS", xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, outline=F)
axis(1, at=c(1,2), labels=c("Mitótico", "Haplóide"), cex.axis=1.5, lwd=0)
text(x=c(1.5), y=c(1.5), labels=c("***"), cex=3)
text(x=c(1,2), y=c((boxplot.stats(mit$dnds)$stats[5]+0.1), (boxplot.stats(hap$dnds)$stats[5]+0.1)), labels = c(length(mit$dnds), length(hap$dnds)))
mtext("(b)", 3, line=0.5, at=0, cex=1.5, xpd=TRUE)
## % ##
barplot(bab, col = c("orchid4","orchid1"), xpd=F, xlab = "", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Antigos", "Novos"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(new.old.a.interest[2,2], new.old.a.interest[1,2], new.old.a.interest[2,1], new.old.a.interest[1,1]), cex=1)
legend(x=6, y=60, inset=c(-5,-0.5),legend=c("com seleção", "sem seleção"), fill=c("orchid4","orchid1"), bty="n", cex=1.5, xpd = T)
mtext("(c)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
barplot(beb, col = c("orchid4","orchid1"), xpd=F, xlab = "", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis = 1.5, ylim=c(90,100), xaxt="n")
axis(1, at=c(0.65, 1.95), labels=c("Mitótico", "Haplóide"), cex.axis=1.5, lwd=0)
text(x=c(1.3), y=c(97), labels=c("***"), cex=3)
segments(0.65,96.5,1.95,96.5, cex=2, lwd=4)
text(x=c(0.65, 0.65, 1.95, 1.95), y=c(90.5, 99.5, 90.5, 99.5), labels=c(mit.hap.a.test[2,2], mit.hap.a.test[1,2], mit.hap.a.test[2,1], mit.hap.a.test[1,1]), cex=1)
mtext("(d)", 3, line=1.1, at=0, cex=1.5, xpd=TRUE)
dev.off()

######################## Figure 3 ####################################
# bibliotecas:
library(ggplot2)
library(plyr)
# dados a seres utilizados: dados totais X e autossomicos dos genes dos grupos mitotico, meiotico e pos-meiotico.
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
aall <- length(subset(total$XorA, total$XorA=="A" & total$Group!="Equal"))/length(subset(total$XorA, total$Group!="Equal"))
newaall <- length(subset(total$XorA, total$XorA=="A" & total$age=="new" & total$Group!="Equal"))/length(subset(total$XorA, total$age=="new" & total$Group!="Equal"))
x.count <- c(subset(x.datab$Count, x.datab$Age=="old" & x.datab$Group=="Mitotic"), subset(x.datab$Count, x.datab$Age=="old" & x.datab$Group!="Mitotic"), subset(x.datab$Count, x.datab$Age=="new" & x.datab$Group=="Mitotic"), subset(x.datab$Count, x.datab$Age=="new" & x.datab$Group!="Mitotic"))
a.count <- c(subset(a.datab$Count, a.datab$Age=="old" & a.datab$Group=="Mitotic"), subset(a.datab$Count, a.datab$Age=="old" & a.datab$Group!="Mitotic"), subset(a.datab$Count, x.datab$Age=="new" & a.datab$Group=="Mitotic"), subset(a.datab$Count, a.datab$Age=="new" & a.datab$Group!="Mitotic"))
# estatísticas:
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Mitotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Mitotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.3341 | mmpp: 0.09214 | al: 0.1246
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='Mitotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Mitotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Mitotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.3389 | mmpp: 0.09631 | al: 0.132
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Meiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Meiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.001026 | mmpp: 0.001512 | al: 0.002704 
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='Meiotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Meiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='Meiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.001365 | mmpp: 0.002063 | al: 0.003811
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='PostMeiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='PostMeiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.6651 | mmpp: 0.4619 | al: 0.4684
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='old' & datab$Group=='PostMeiotic'), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='PostMeiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='old' & datab$Group!='PostMeiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.6855 | mmpp: 0.4933 | al: 0.501
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Mitotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Mitotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.06279 | mmpp: 0.0504 | al: 0.01002
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='Mitotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Mitotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Mitotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.07479 | mmpp: 0.05897 | al: 0.01308
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Meiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Meiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.7945 | mmpp: 0.668 | al: 0.6996
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='Meiotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Meiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='Meiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.8164 | mmpp: 0.7776 | al: 0.7194
fisher.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='PostMeiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='PostMeiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.016 | mmpp: 0.01609 | al: 0.06117
chisq.test(matrix(data=c(subset(datab$Count, datab$Age=='new' & datab$Group=='PostMeiotic'), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='PostMeiotic' & datab$Chromosome=='A')), sum(subset(datab$Count, datab$Age=='new' & datab$Group!='PostMeiotic' & datab$Chromosome=='X'))), nrow = 2)) # tot: 0.02671 | mmpp: 0.02031 | al: 0.07154
# dados e funções pros gráficos:
hline.data <- data.frame(z = c(aall, newaall), Age = c("old","new"))
labs.data <- data.frame(s=c(2, 1, 3), f=c(0.878, 0.74, 0.905), z = c("***", "", "***"), Age = c("old","new", "new"), Chromosome=c("A", "A", "A"))
ps <- data.frame(a=c(1,3,1,2), b=c(0.834, 0.834, 0.74, 0.815), d=c("0.096", "0.493", "0.059", "0.777"), Age=c("old", "old", "new", "new"), Chromosome=c("A", "A", "A", "A"))
texto <- data.frame(x=c(rep(c(1,2,3),4)), y=c(rep(0.99, 6), rep(0.61, 6)), lab=c(x.count, a.count), Age=c(rep(c(rep("old",3), rep("new", 3)),2)))
facet_names <- list(
    'old'="Old Genes",
    'new'="New Genes")
facet_names_pt <- list(
    'old'="Genes Antigos",
    'new'="Genes Novos")
facet_labeller <- function(variable, value){
    return(facet_names[value])
}
facet_labeller_pt <- function(variable, value){
    return(facet_names_pt[value])
}
########## BW ##########
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
pdf("figure3_bw_pt.pdf", width=15, height=10)#, res=300) # width=22, height=10, 
ggplot(datab, aes(x=Group, y=Proportion, fill=Chromosome)) +
    geom_bar(position='stack', stat='identity') +
    facet_grid(.~Age, labeller=facet_labeller_pt) +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme(text = element_text(size=30)) +
    scale_x_discrete(name="") +
    geom_hline(aes(yintercept =z), hline.data) +
    geom_text(aes(x=x,y=y,label=lab), texto, size=7, inherit.aes=F) +
    geom_text(aes(x=s, y=f, label=z), labs.data, size=20) +
    geom_text(aes(x=a, y=b, label=d), ps, size=5) +
    scale_fill_manual(values=c("Grey 70", "Gray 50"), name="", breaks=c("A", "X"), labels=c("Genes Autossômicos", "Genes ligados ao X"))
dev.off()
########## Color ##########
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
pdf("figure3_color_pt.pdf", width=15, height=10)#, res=300) # width=22, height=10, 
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
    scale_fill_manual(values=c("Khaki", "Thistle"), name="", breaks=c("A", "X"), labels=c("Genes Autossômicos", "Genes ligados ao X"))
dev.off()


######################## Figure 4 ##################################
# dados a serem usados: genes totais autossomicos para todas as classes usadas
total.a$Class <- as.numeric(total.a$Class)
classes <- subset(total.a, total.a$Class<=13)
classes$Class = factor(classes$Class) ### "joga fora" os factors vazios
clx = table(classes$age, classes$Class)
clp <- prop.table(clx, 1)*100
########## BW ##########
pdf("figure4_bw.pdf", width = 9)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(clp, beside=T, col=c("dimgray","gray"), xlab="Numerical Class", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,38, by=3)), labels=c(seq(1,13,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 
pdf("figure4_bw_pt.pdf", width = 9)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(clp, beside=T, col=c("dimgray","gray"), xlab="Classes Numéricas", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,38, by=3)), labels=c(seq(1,13,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 
########## Color ##########
pdf("figure4_color.pdf", width = 9)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(clp, beside=T, col=c("powderblue","salmon"), xlab="Numerical Class", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,38, by=3)), labels=c(seq(1,13,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()
pdf("figure4_color_pt.pdf", width = 9)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(clp, beside=T, col=c("powderblue","salmon"), xlab="Classes Numéricas", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,38, by=3)), labels=c(seq(1,13,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()



######################## Figure 5 ##################################
# dados a serem usados: genes totais autossomicos para todas as classes
total.a$Class <- factor(total.a$Class)
ttx = table(total.a$age, total.a$Class)
ttp <- prop.table(ttx, 1)*100
########## BW ##########
pdf("figure5_bw.pdf", width = 13)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(ttp, beside=T, col=c("dimgray","gray"), xlab="Numerical Class", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,56, by=3)), labels=c(seq(1,19,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off() 
pdf("figure5_bw_pt.pdf", width = 13)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(ttp, beside=T, col=c("dimgray","gray"), xlab="Classes Numéricas", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,56, by=3)), labels=c(seq(1,19,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
dev.off()
########## Color ##########
pdf("figure5_color.pdf", width = 13)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(ttp, beside=T, col=c("powderblue","salmon"), xlab="Numerical Class", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,56, by=3)), labels=c(seq(1,19,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("new genes", "old genes"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()
pdf("figure5_color_pt.pdf", width = 13)#, res=300)
par(mar=c(6,5,4,2)+0.1)
barplot(ttp, beside=T, col=c("powderblue","salmon"), xlab="Classes Numéricas", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,40), xpd=F, xaxt="n")
axis(1, at=c(seq(2,56, by=3)), labels=c(seq(1,19,by=1)), cex.axis=1.5, lwd=0)
legend(x=3, y=35, inset=c(-5,-0.5),legend=c("genes novos", "genes antigos"), fill=c("salmon","powderblue"), bty="n", cex=1.5, xpd = T)
dev.off()


######################## Figure 6 ##################################
# dados a serem usados: genes ligados ao X e autossomicos para as classes de interesse
xatb <- subset(total, total$Group=='Mitotic' | total$Group=='Meiotic' | total$Group=='PostMeiotic')
xatb$Group = factor(xatb$Group) ### "joga fora" os factors vazios
xax <- table(xatb$XorA, xatb$Group)
xap <- prop.table(xax, 1)*100
xaxn <- table(subset(xatb$XorA, xatb$age=='new'), subset(xatb$Group, xatb$age=='new'))
xapn <- prop.table(xaxn, 1)*100
xaxo <- table(subset(xatb$XorA, xatb$age=='old'), subset(xatb$Group, xatb$age=='old'))
xapo <- prop.table(xaxo, 1)*100
########## BW ##########
pdf("figure6_bw.pdf", width = 13)#, res=300)
par(mfrow=c(1,3), mar=c(6,5,4,2)+0.1)
barplot(xap, beside=T, col=c("gray","dimgray"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main = "All genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
barplot(xapn, beside=T, col=c("gray","dimgray"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="New Genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
legend(x=1.5, y=60, inset=c(-5,-5),legend=c("autossomic genes", "X-linked genes"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
barplot(xapo, beside=T, col=c("gray","dimgray"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Old Genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
dev.off() 
pdf("figure6_bw_pt.pdf", width = 13)#, res=300)
par(mfrow=c(1,3), mar=c(6,5,4,2)+0.1)
barplot(xap, beside=T, col=c("gray","dimgray"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main = "Todos os genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
barplot(xapn, beside=T, col=c("gray","dimgray"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Genes Novos")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
legend(x=1.5, y=60, inset=c(-5,-5),legend=c("genes autossômicos", "genes ligados ao X"), fill=c("gray","dimgray"), bty="n", cex=1.5, xpd = T)
barplot(xapo, beside=T, col=c("gray","dimgray"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Genes Antigos")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
dev.off() 
########## Color ##########
pdf("figure6_color.pdf", width = 13)#, res=300)
par(mfrow=c(1,3), mar=c(6,5,4,2)+0.1)
barplot(xap, beside=T, col=c("Khaki","Thistle"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main = "All genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
barplot(xapn, beside=T, col=c("Khaki","Thistle"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="New Genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
legend(x=1.5, y=60, inset=c(-5,-5),legend=c("autossomic genes", "X-linked genes"), fill=c("Khaki","Thistle"), bty="n", cex=1.5, xpd = T)
barplot(xapo, beside=T, col=c("Khaki","Thistle"), xlab="Spermatogenesis Phase", ylab="Percentage of genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Old Genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitotic", "Meiotic", "PostMeiotic"), cex.axis=1.5, lwd=0)
dev.off() 
pdf("figure6_color_pt.pdf", width = 13)#, res=300)
par(mfrow=c(1,3), mar=c(6,5,4,2)+0.1)
barplot(xap, beside=T, col=c("Khaki","Thistle"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main = "Todos os genes")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
barplot(xapn, beside=T, col=c("Khaki","Thistle"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Genes novos")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
legend(x=1.5, y=60, inset=c(-5,-5),legend=c("genes autossômicos", "genes ligados ao X"), fill=c("Khaki","Thistle"), bty="n", cex=1.5, xpd = T)
barplot(xapo, beside=T, col=c("Khaki","Thistle"), xlab="Fase da Espermatogênese", ylab="Porcentagem de genes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, ylim=c(0,65), xpd=F, xaxt="n", main="Genes antigos")
axis(1, at=c(seq(2,8, by=3)), labels=c("Mitótica", "Meiótica", "PósMeiótica"), cex.axis=1.5, lwd=0)
dev.off() 


########## Figure 7 ##########
# grafico da expressão de genes novos ao longo da espermatogenese
# bibliotecas necessarias:
library(ggplot2)
library(gridExtra)
# dados usados: genes autossomicos das classes 1 a 10 e 12
total.a$Class <- factor(total.a$Class)
interest <- subset(total.a, total.a$Class!='13' & total.a$Class!='19' & total.a$Class!='18' & total.a$Class!='17' & total.a$Class!='16' & total.a$Class!='15' & total.a$Class!='14')
interest$Class <- factor(interest$Class)
# estatisticas
wilcox.test(interest$Mitosis~interest$age)
wilcox.test(interest$Meiosis~interest$age)
wilcox.test(interest$PostMeiosis~interest$age)
# gráficos:
########## BW Eng ##########
mit <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Mitosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Genes' age", y="Expression during Mitosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("0.1725",length(subset(interest$Mitosis, interest$age=="old")),length(subset(interest$Mitosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
mei <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Meiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Genes' age", y="Expression during Meiosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$Meiosis, interest$age=="old")), length(subset(interest$Meiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pm <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$PostMeiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Genes' age", y="Expression during PostMeiosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$PostMeiosis, interest$age=="old")), length(subset(interest$PostMeiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pdf("figure7_bw.pdf")
grid.arrange(arrangeGrob(mit, mei, pm,  nrow=1, ncol=3))
dev.off() 
########## BW Pt ##########
mit <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Mitosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Idade dos Genes", y="Expressão durante a Mitose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("0.1725",length(subset(interest$Mitosis, interest$age=="old")),length(subset(interest$Mitosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
mei <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Meiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Idade dos Genes", y="Expressão durante a Meiose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$Meiosis, interest$age=="old")), length(subset(interest$Meiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pm <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$PostMeiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("dimgrey", "grey")) +
    labs(x="Idade dos Genes", y="Expressão durante a Pós-Meiose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$PostMeiosis, interest$age=="old")), length(subset(interest$PostMeiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pdf("figure7_bw_pt.pdf")
grid.arrange(arrangeGrob(mit, mei, pm,  nrow=1, ncol=3))
dev.off() 
########## Color Eng ##########
mit <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Mitosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Genes' age", y="Expression during Mitosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("0.1725",length(subset(interest$Mitosis, interest$age=="old")),length(subset(interest$Mitosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
mei <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Meiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Genes' age", y="Expression during Meiosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$Meiosis, interest$age=="old")), length(subset(interest$Meiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pm <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$PostMeiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Genes' age", y="Expression during PostMeiosis") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$PostMeiosis, interest$age=="old")), length(subset(interest$PostMeiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pdf("figure7_color.pdf")
grid.arrange(arrangeGrob(mit, mei, pm,  nrow=1, ncol=3))
dev.off() 
########## Color Pt ##########
mit <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Mitosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Idade dos Genes", y="Expressão durante a Mitose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("0.1725",length(subset(interest$Mitosis, interest$age=="old")),length(subset(interest$Mitosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
mei <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$Meiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Idade dos Genes", y="Expressão durante a Meiose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$Meiosis, interest$age=="old")), length(subset(interest$Meiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pm <- ggplot(data.frame(interest), aes(x=interest$age, y=interest$PostMeiosis, color=age)) +
    geom_boxplot(notch = T) +
    scale_color_manual(values=c("powderblue", "salmon")) +
    labs(x="Idade dos Genes", y="Expressão durante a Pós-Meiose") + 
    theme(legend.position = "none") + annotate("text", x=c(1.5,1,2), y=c(13.5,3,3), label=c("***", length(subset(interest$PostMeiosis, interest$age=="old")), length(subset(interest$PostMeiosis, interest$age=="new")))) +ggtitle("") +
    scale_y_continuous(limits = c(3, 15))
pdf("figure7_color_pt.pdf")
grid.arrange(arrangeGrob(mit, mei, pm,  nrow=1, ncol=3))
dev.off() 

