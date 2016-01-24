# January 2016
# Final R program to make:
# - make controls
# - make graphs 1, 2, 3, 4, suplementals et al

al <- read.table("final.output", header=T)
al$age <- factor(al$age, levels=c("old","new"))