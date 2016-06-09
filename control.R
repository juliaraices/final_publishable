# June 2016
# Júlia Raíces

droso <- read.table("final.output", header=T)

equal <- subset(droso, droso$Group=="Equal")

meiotic <- subset(droso, droso$Group=="Meiotic")

post <- subset(droso, droso$Group=="PostMeiotic")

haploid <- subset(droso, droso$Group=="Meiotic" | droso$Group=="PostMeiotic" | droso$Group=="MeioticPostmeiotic")

bootz <- rep("NA", 100)

for(i in 1:2){
    repeat{
        temp <- sample(equal, length(meiotic))
        blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
        if(blabs$p.value >= 0.05) break
    }
    will <- wilcox.test(temp$alpha, meiotic$alpha)
    bootz <- will$p.value
}

bootz

