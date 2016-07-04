# July 2016
# Julia Raices

droso <- read.table("final.output", header=T)

equal <- subset(droso, droso$Group=="Equal")

meiotic <- subset(droso, droso$Group=="Meiotic")

post <- subset(droso, droso$Group=="PostMeiotic")

haploid <- subset(droso, droso$Group=="Meiotic" | droso$Group=="PostMeiotic" | droso$Group=="MeioticPostmeiotic")

bootz <- rep("NA", 100)

a = 0
b = 0

summary(meiotic$Meiosis)

length(subset(meiotic$Meiosis, meiotic$Meiosis >= 4 & meiotic$Meiosis < 5))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 5 & meiotic$Meiosis < 6))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 6 & meiotic$Meiosis < 7))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 7 & meiotic$Meiosis < 8))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 8 & meiotic$Meiosis < 9))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 9 & meiotic$Meiosis < 10))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 10 & meiotic$Meiosis < 11))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 11 & meiotic$Meiosis < 12))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 12 & meiotic$Meiosis < 13))
length(subset(meiotic$Meiosis, meiotic$Meiosis >= 13 & meiotic$Meiosis < 14))

z <- (subset(meiotic, meiotic$Meiosis >= 12 & meiotic$Meiosis < 13))
w <- (subset(meiotic, meiotic$Meiosis >= 13 & meiotic$Meiosis < 14))

y <- rbind(z, w)

for(i in 1:2){
    eq <- equal
    ea <- subset(eq, eq$Meiosis >= 4 & eq$Meiosis < 5)
    eb <- subset(eq, eq$Meiosis >= 5 & eq$Meiosis < 6)
    ec <- subset(eq, eq$Meiosis >= 6 & eq$Meiosis < 7)
    ed <- subset(eq, eq$Meiosis >= 7 & eq$Meiosis < 8)
    ee <- subset(eq, eq$Meiosis >= 8 & eq$Meiosis < 9)
    ef <- subset(eq, eq$Meiosis >= 9 & eq$Meiosis < 10)
    eg <- subset(eq, eq$Meiosis >= 10 & eq$Meiosis < 11)
    eh <- subset(eq, eq$Meiosis >= 11 & eq$Meiosis < 12)
    ei <- subset(eq, eq$Meiosis >= 12 & eq$Meiosis < 13)
    ej <- subset(eq, eq$Meiosis >= 13 & eq$Meiosis < 14)
    ta <- ea[sample(nrow(ea), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 4 & meiotic$Meiosis < 5))), ]
    tb <- eb[sample(nrow(eb), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 5 & meiotic$Meiosis < 6))), ]
    tc <- ec[sample(nrow(ec), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 6 & meiotic$Meiosis < 7))), ]
    td <- ed[sample(nrow(ed), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 7 & meiotic$Meiosis < 8))), ]
    te <- ee[sample(nrow(ee), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 8 & meiotic$Meiosis < 9)))]
    tf <- ef[sample(nrow(ef), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 9 & meiotic$Meiosis < 10))), ]
    tg <- eg[sample(nrow(eg), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 10 & meiotic$Meiosis < 11))), ]
    th <- eh[sample(nrow(eh), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 11 & meiotic$Meiosis < 12))), ]
    ti <- ei[sample(nrow(ei), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 12 & meiotic$Meiosis < 13))), ]
    tj <- ej[sample(nrow(ej), length(subset(meiotic$Meiosis, meiotic$Meiosis >= 13 & meiotic$Meiosis < 14))), ]
}


for(i in 1:2){
    repeat{
        eq <- equal
        temp <- eq[sample(nrow(eq), 50), ]
        blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
        a <- a+1
        print(a)
        if(blabs$p.value >= 0.05) break
    }
    a=0
    b <- b+1
    will <- wilcox.test(temp$alpha, meiotic$alpha)
    bootz[b] <- will$p.value
}




