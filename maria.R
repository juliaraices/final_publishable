tot <- read.table("final.output", header = T, sep = "\t")
tots <- subset(tot, tot$dn != "NA")
tot.a <- subset(tot, tot$XorA=="A")
tot.x <- subset(tot, tot$XorA=="X")

new.a <- subset(tot.a, tot.a$age=="new")
old.a <- subset(tot.a, tot.a$age=="old")

new.x <- subset(tot.x, tot.x$age=="new")
old.x <- subset(tot.x, tot.x$age=="old")

news <- subset(tot, tot$age=="new")
olds <- subset(tot, tot$age=="old")

mit <- subset(tot, tot$Group=="Mitotic")
mit.a <- subset(mit, mit$XorA=="A")
mit.new.a <- subset(mit.a, mit.a$age=="new")
mit.old.a <- subset(mit.a, mit.a$age=="old")
mit.x <- subset(mit, mit$XorA=="X")
mit.new.x <- subset(mit.x, mit.x$age=="new")
mit.old.x <- subset(mit.x, mit.x$age=="old")

mei <- subset(tot, tot$Group=="Meiotic")
mei.a <- subset(mei, mei$XorA=="A")
mei.new.a <- subset(mei.a, mei.a$age=="new")
mei.old.a <- subset(mei.a, mei.a$age=="old")
mei.x <- subset(mei, mei$XorA=="X")
mei.new.x <- subset(mei.x, mei.x$age=="new")
mei.old.x <- subset(mei.x, mei.x$age=="old")

pmei <- subset(tot, tot$Group=="PostMeiotic")
pmei.a <- subset(pmei, pmei$XorA=="A")
pmei.new.a <- subset(pmei.a, pmei.a$age=="new")
pmei.old.a <- subset(pmei.a, pmei.a$age=="old")
pmei.x <- subset(pmei, pmei$XorA=="X")
pmei.new.x <- subset(pmei.x, pmei.x$age=="new")
pmei.old.x <- subset(pmei.x, pmei.x$age=="old")

length(tot$id)

length(news$id)
length(subset(tot$id, tot$age=="new"))

length(olds$id)
length(subset(tot$id, tot$age=="old"))

length(tots$id)
length(subset(tot$id, tot$dn!="NA"))

### Expressão ####

summary(new.a$Mitosis)
summary(subset(tot.a$Mitosis, tot.a$age=="new"))
summary(subset(tot$Mitosis, tot$age=="new" & tot$XorA=="A"))

summary(old.a$Mitosis)
summary(subset(tot.a$Mitosis, tot.a$age=="old"))
summary(subset(tot$Mitosis, tot$age=="old" & tot$XorA=="A"))

length(mit.new.a$Group)
length(subset(new.a$id, new.a$Group=="Mitotic"))
length(subset(tot.a$Mitosis, tot.a$age=="new" & tot.a$Group=="Mitotic"))
length(subset(tot$Mitosis, tot$age=="new" & tot$XorA=="A" & tot$Group=="Mitotic"))

length(mit.old.a$Mitosis)
length(subset(old.a$id, old.a$Group=="Mitotic"))
length(subset(tot.a$Mitosis, tot.a$age=="old" & tot.a$Group=="Mitotic"))
length(subset(tot$Mitosis, tot$age=="old" & tot$XorA=="A" & tot$Group=="Mitotic"))


summary(new.a$Meiosis)
summary(subset(tot.a$Meiosis, tot.a$age=="new"))
summary(subset(tot$Meiosis, tot$age=="new" & tot$XorA=="A"))

summary(old.a$Meiosis)
summary(subset(tot.a$Meiosis, tot.a$age=="old"))
summary(subset(tot$Meiosis, tot$age=="old" & tot$XorA=="A"))

length(mei.new.a$Meiosis)
length(subset(tot.a$Meiosis, tot.a$age=="new" & tot.a$Group=="Meiotic"))
length(subset(tot$Meiosis, tot$age=="new" & tot$XorA=="A"& tot$Group=="Meiotic"))

length(mei.old.a$Meiosis)
length(subset(tot.a$Meiosis, tot.a$age=="old" & tot.a$Group=="Meiotic"))
length(subset(tot$Meiosis, tot$age=="old" & tot$XorA=="A" & tot$Group=="Meiotic"))


summary(new.a$PostMeiosis)
summary(subset(tot.a$PostMeiosis, tot.a$age=="new"))
summary(subset(tot$PostMeiosis, tot$age=="new" & tot$XorA=="A"))

summary(old.a$PostMeiosis)
summary(subset(tot.a$PostMeiosis, tot.a$age=="old"))
summary(subset(tot$PostMeiosis, tot$age=="old" & tot$XorA=="A"))

length(pmei.new.a$PostMeiosis)
length(subset(tot.a$PostMeiosis, tot.a$age=="new" & tot.a$Group=="PostMeiotic"))
length(subset(tot$PostMeiosis, tot$age=="new" & tot$XorA=="A" & tot$Group=="PostMeiotic"))

length(pmei.old.a$PostMeiosis)
length(subset(tot.a$PostMeiosis, tot.a$age=="old" & tot.a$Group=="PostMeiotic"))
length(subset(tot$PostMeiosis, tot$age=="old" & tot$XorA=="A"  & tot$Group=="PostMeiotic"))



### X and A ###
length(mit.x$id)
length(subset(mit$id, mit$XorA=="X"))
length(subset(tot$id, tot$XorA=="X" & tot$Group=="Mitotic"))

length(mei.x$id)
length(subset(mei$id, mei$XorA=="X"))
length(subset(tot$id, tot$XorA=="X" & tot$Group=="Meiotic"))

length(pmei.x$id)
length(subset(pmei$id, pmei$XorA=="X"))
length(subset(tot$id, tot$XorA=="X" & tot$Group=="PostMeiotic"))


length(mit.a$id)
length(subset(mit$id, mit$XorA=="A"))
length(subset(tot$id, tot$XorA=="A" & tot$Group=="Mitotic"))

length(mei.a$id)
length(subset(mei$id, mei$XorA=="A"))
length(subset(tot$id, tot$XorA=="A" & tot$Group=="Meiotic"))

length(pmei.a$id)
length(subset(pmei$id, pmei$XorA=="A"))
length(subset(tot$id, tot$XorA=="A" & tot$Group=="PostMeiotic"))


length(mit.new.x$id)
length(subset(mit.x$id, mit.x$age=="new"))
length(subset(mit$id, mit$age=="new" & mit$XorA=="X"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="X" & tot$Group=="Mitotic"))

length(mei.new.x$id)
length(subset(mei.x$id, mei.x$age=="new"))
length(subset(mei$id, mei$age=="new" & mei$XorA=="X"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="X" & tot$Group=="Meiotic"))

length(pmei.new.x$id)
length(subset(pmei.x$id, pmei.x$age=="new"))
length(subset(pmei$id, pmei$age=="new" & pmei$XorA=="X"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="X" & tot$Group=="PostMeiotic"))


length(mit.old.x$id)
length(subset(mit.x$id, mit.x$age=="old"))
length(subset(mit$id, mit$age=="old" & mit$XorA=="X"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="X" & tot$Group=="Mitotic"))

length(mei.old.x$id)
length(subset(mei.x$id, mei.x$age=="old"))
length(subset(mei$id, mei$age=="old" & mei$XorA=="X"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="X" & tot$Group=="Meiotic"))

length(pmei.old.x$id)
length(subset(pmei.x$id, pmei.x$age=="old"))
length(subset(pmei$id, pmei$age=="old" & pmei$XorA=="X"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="X" & tot$Group=="PostMeiotic"))


length(mit.new.a$id)
length(subset(mit.a$id, mit.a$age=="new"))
length(subset(mit$id, mit$age=="new" & mit$XorA=="A"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="A" & tot$Group=="Mitotic"))

length(mei.new.a$id)
length(subset(mei.a$id, mei.a$age=="new"))
length(subset(mei$id, mei$age=="new" & mei$XorA=="A"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="A" & tot$Group=="Meiotic"))

length(pmei.new.a$id)
length(subset(pmei.a$id, pmei.a$age=="new"))
length(subset(pmei$id, pmei$age=="new" & pmei$XorA=="A"))
length(subset(tot$id, tot$age=="new" & tot$XorA=="A" & tot$Group=="PostMeiotic"))


length(mit.old.a$id)
length(subset(mit.a$id, mit.a$age=="old"))
length(subset(mit$id, mit$age=="old" & mit$XorA=="A"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="A" & tot$Group=="Mitotic"))

length(mei.old.a$id)
length(subset(mei.a$id, mei.a$age=="old"))
length(subset(mei$id, mei$age=="old" & mei$XorA=="A"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="A" & tot$Group=="Meiotic"))

length(pmei.old.a$id)
length(subset(pmei.a$id, pmei.a$age=="old"))
length(subset(pmei$id, pmei$age=="old" & pmei$XorA=="A"))
length(subset(tot$id, tot$age=="old" & tot$XorA=="A" & tot$Group=="PostMeiotic"))




### dN/dS ###
summary(new.a$dnds)
summary(subset(news$dnds, news$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="new"))
summary(subset(tot.a$dnds, tot.a$age=="new"))

summary(old.a$dnds)
summary(subset(olds$dnds, olds$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="old"))
summary(subset(tot.a$dnds, tot.a$age=="old"))


summary(mit.new.a$dnds)
summary(subset(mit.a$dnds, mit.a$age=="new"))
summary(subset(mit$dnds, mit$age=="new" & mit$XorA=="A"))
summary(subset(new.a$dnds, new.a$Group=="Mitotic"))
summary(subset(news$dnds, news$Group=="Mitotic" & news$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="new" & tot$Group=="Mitotic"))
summary(subset(tot.a$dnds, tot.a$age=="new" & tot.a$Group=="Mitotic"))

summary(mit.old.a$dnds)
summary(subset(mit.a$dnds, mit.a$age=="old"))
summary(subset(mit$dnds, mit$age=="old" & mit$XorA=="A"))
summary(subset(old.a$dnds, old.a$Group=="Mitotic"))
summary(subset(olds$dnds, olds$Group=="Mitotic" & olds$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="old" & tot$Group=="Mitotic"))
summary(subset(tot.a$dnds, tot.a$age=="old" & tot.a$Group=="Mitotic"))


summary(mei.new.a$dnds)
summary(subset(mei.a$dnds, mei.a$age=="new"))
summary(subset(mei$dnds, mei$age=="new" & mei$XorA=="A"))
summary(subset(new.a$dnds, new.a$Group=="Meiotic"))
summary(subset(news$dnds, news$Group=="Meiotic" & news$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="new" & tot$Group=="Meiotic"))
summary(subset(tot.a$dnds, tot.a$age=="new" & tot.a$Group=="Meiotic"))

summary(mei.old.a$dnds)
summary(subset(mei.a$dnds, mei.a$age=="old"))
summary(subset(mei$dnds, mei$age=="old" & mei$XorA=="A"))
summary(subset(old.a$dnds, old.a$Group=="Meiotic"))
summary(subset(olds$dnds, olds$Group=="Meiotic" & olds$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="old" & tot$Group=="Meiotic"))
summary(subset(tot.a$dnds, tot.a$age=="old" & tot.a$Group=="Meiotic"))


summary(pmei.new.a$dnds)
summary(subset(pmei.a$dnds, pmei.a$age=="new"))
summary(subset(pmei$dnds, pmei$age=="new" & pmei$XorA=="A"))
summary(subset(new.a$dnds, new.a$Group=="PostMeiotic"))
summary(subset(news$dnds, news$Group=="PostMeiotic" & news$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="new" & tot$Group=="PostMeiotic"))
summary(subset(tot.a$dnds, tot.a$age=="new" & tot.a$Group=="PostMeiotic"))

summary(pmei.old.a$dnds)
summary(subset(pmei.a$dnds, pmei.a$age=="old"))
summary(subset(pmei$dnds, pmei$age=="old" & pmei$XorA=="A"))
summary(subset(old.a$dnds, old.a$Group=="PostMeiotic"))
summary(subset(olds$dnds, olds$Group=="PostMeiotic" & olds$XorA=="A"))
summary(subset(tot$dnds, tot$XorA=="A" & tot$age=="old" & tot$Group=="PostMeiotic"))
summary(subset(tot.a$dnds, tot.a$age=="old" & tot.a$Group=="PostMeiotic"))

