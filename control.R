# June 2016
# Júlia Raíces

droso <- read.table("final.output", header=T)

equal <- subset(droso, droso$Group=="Equal")

meiotic <- subset(droso, droso$Group=="Meiotic")

post <- subset(droso, droso$Group=="PostMeiotic")

haploid <- subset(droso, droso$Group=="Meiotic" | droso$Group=="PostMeiotic" | droso$Group=="MeioticPostmeiotic")

bootz <- rep("NA", 100)

<<<<<<< HEAD
for(i in 1:2){
=======
# 1st try
a = 0
for(i in 1:100){
>>>>>>> 518dd22773e1c81c88ec925843b9592bc37f7d4c
    repeat{
        eq <- equal
        temp <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
        blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
        a <- a+1
        print(a)
        if(blabs$p.value >= 0.05) break
    }
    will <- wilcox.test(temp$alpha, meiotic$alpha)
    bootz <- will$p.value
}

<<<<<<< HEAD
bootz
=======


# 2nd try
b = 0
for(i in 1:100){
    while(TRUE){
        eq <- equal
        temp <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
        blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
        b <- b+1
        print(b)
        if(blabs$p.value >= 0.05) break
    }
    will <- wilcox.test(temp$alpha, meiotic$alpha)
    bootz <- will$p.value
}

# 3nd try
c = 0
random.sample <- function(interest) {
    eq <- equal
    x <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
    t <- wilcox.test(x$Meiosis, meiotic$Meiosis)
    t$p.value >= 0.05
    c <- c + 1
    print(c)
    if(t$p.value >= 0.05) return(x)
    else Recall(x)# run the function again
}

for(i in 1:100){
    y <- random.sample(equal)
    will <- wilcox.test(y$alpha, meiotic$alpha)
    bootz <- will$p.value
}




# 4nd try
d = 0
sucess <- FALSE
for(i in 1:100){
    while(!sucess){
        eq <- equal
        temp <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
        blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
        d <- d+1
        print(d)
        sucess <- blabs$p.value >= 0.05
    }
    will <- wilcox.test(temp$alpha, meiotic$alpha)
    bootz <- will$p.value
}



>>>>>>> 518dd22773e1c81c88ec925843b9592bc37f7d4c

