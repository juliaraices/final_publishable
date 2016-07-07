# July 2016
# Julia Raices

# creating new controsl as the old ones may not be so good

# reading the file and creating interest groups (as in final.R)
droso <- read.table("final.output", header=T)

equal <- subset(droso, droso$Group=="Equal")
eq.new <- subset(equal, equal$Age=="new")
eq.old <- subset(equal, equal$Age=="old")

meiotic <- subset(droso, droso$Group=="Meiotic")
mei.new <- subset(meiotic, meiotic$Age=="new")
mei.old <- subset(meiotic, meiotic$Age=="old")

post <- subset(droso, droso$Group=="PostMeiotic")

haploid <- subset(droso, droso$Group=="Meiotic" | droso$Group=="PostMeiotic" | droso$Group=="MeioticPostmeiotic")

bootz <- rep("NA", 100)

a = 0
b = 0

summary(meiotic$Meiosis)

y <- rbind(z, w)

# trying to make a bootstrap in wich the subset of equal new genes has the same expression as the meiotic new genes (in this case) and check if the alpha is the same or not.

for(i in 1:100){
    repeat{
        eq <- eq.new
        temp <- eq[sample(nrow(eq), 200), ]
        blabs <- wilcox.test(temp$Meiosis, mei.new$Meiosis)
        a <- a+1
        print(a)
        if(blabs$p.value >= 0.05) break
    }
    a=0
    b <- b+1
    will <- wilcox.test(temp$alpha, mei.new$alpha)
    bootz[b] <- will$p.value
}

# this takes a hell of long time to run... may be I did sth wrong...


