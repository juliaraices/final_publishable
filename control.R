# June 2016
# Júlia Raíces

# creating new controsl as the old ones may not be so good

# reading the file and creating interest groups (as in final.R)

droso <- read.table("final.output", header=T)

equal <- subset(droso, droso$Group=="Equal")

meiotic <- subset(droso, droso$Group=="Meiotic")

post <- subset(droso, droso$Group=="PostMeiotic")

haploid <- subset(droso, droso$Group=="Meiotic" | droso$Group=="PostMeiotic" | droso$Group=="MeioticPostmeiotic")

bootz <- rep("NA", 100)

# trying to make a bootstrap in wich the subset of equal genes has the same expression as the meiotic genes (in this case) and check if the alpha is the same or not.
# 1st try
a = 0
b = 0
for(i in 1:100){
    repeat{
        eq <- equal
        temp <- eq[sample(nrow(eq), 500, replace = TRUE), ]
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

bootz

sink("out_cross.txt")
print("Output do programa cross.R\n\n\n")
print(bootz)

sink()

# 
# # 2nd try
# b = 0
# for(i in 1:100){
#     while(TRUE){
#         eq <- equal
#         temp <- eq[sample(nrow(eq), 2), ]
#         blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
#         b <- b+1
#         print(b)
#         if(blabs$p.value >= 0.05) break
#     }
#     will <- wilcox.test(temp$alpha, meiotic$alpha)
#     bootz <- will$p.value
# }
# 
# # 3nd try
# c = 0
# random.sample <- function(interest) {
#     eq <- equal
#     x <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
#     t <- wilcox.test(x$Meiosis, meiotic$Meiosis)
#     t$p.value >= 0.05
#     c <- c + 1
#     print(c)
#     if(t$p.value >= 0.05) return(x)
#     else Recall(x)# run the function again
# }
# 
# for(i in 1:100){
#     y <- random.sample(equal)
#     will <- wilcox.test(y$alpha, meiotic$alpha)
#     bootz <- will$p.value
# }
# 
# 
# 
# 
# # 4nd try
# d = 0
# sucess <- FALSE
# for(i in 1:100){
#     while(!sucess){
#         eq <- equal
#         temp <- eq[sample(nrow(eq), length(meiotic$Meiosis)), ]
#         blabs <- wilcox.test(temp$Meiosis, meiotic$Meiosis)
#         d <- d+1
#         print(d)
#         sucess <- blabs$p.value >= 0.05
#     }
#     will <- wilcox.test(temp$alpha, meiotic$alpha)
#     bootz <- will$p.value
# }
# 
# 
# 
