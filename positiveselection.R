# Julia Raices, April, 2017

crossed <- read.table("output/final.output", header = T)

divas <- read.table("input/melsubgroup_analysis_results_flydivas_v1.2", header=T, sep = "\t")

mydivas <- merge(crossed, divas, by.x="id", by.y = "cg")
new.mydivas <- subset(mydivas, mydivas$age=='new')
old.mydivas <- subset(mydivas, mydivas$age=='old')

new.old <- matrix(data=c(length(subset(new.mydivas$id, new.mydivas$pos78=='does' | new.mydivas$pos_12=='does' | new.mydivas$pos_88a=='does')), (length(new.mydivas$id) - length(subset(new.mydivas$id, new.mydivas$pos78=='does' | new.mydivas$pos_12=='does' | new.mydivas$pos_88a=='does'))), length(subset(old.mydivas$id, old.mydivas$pos78=='does' | old.mydivas$pos_12=='does' | old.mydivas$pos_88a=='does')), (length(old.mydivas$id) - length(subset(old.mydivas$id, old.mydivas$pos78=='does' | old.mydivas$pos_12=='does' | old.mydivas$pos_88a=='does')))), ncol=2)

fisher.test(new.old) # between new and old genes from all classes there is no difference in the number of genes with and without positive selection markers (p = 0.05762)

# no V, impossible or equal
interest <- subset(mydivas, mydivas$Group!="TheV" | mydivas$Group!="Impossible" | mydivas!="Equal")
interest.old <- subset(interest, interest$age=='old')
interest.new <- subset(interest, interest$age=='new')

new.old.interest <- matrix(data=c(length(subset(interest.new$id, interest.new$pos78=='does' | interest.new$pos_12=='does' | interest.new$pos_88a=='does')), (length(interest.new$id) - length(subset(interest.new$id, interest.new$pos78=='does' | interest.new$pos_12=='does' | interest.new$pos_88a=='does'))), length(subset(interest.old$id, interest.old$pos78=='does' | interest.old$pos_12=='does' | interest.old$pos_88a=='does')), (length(interest.old$id) - length(subset(interest.old$id, interest.old$pos78=='does' | interest.old$pos_12=='does' | interest.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.interest) # between new and old genes from interest groups (all but equal, impossible e theV), there is no difference in the number of genes with positive selection markers (p = 0.05762)

# only mitotic and haploid groups (interest groups)
mithap <- subset(mydivas, mydivas$Group=="Mitotic" | mydivas$HaploidGroup=="haploid_group")
mithap.new <- subset(mithap, mithap$age=='new')
mithap.old <- subset(mithap, mithap$age=='old')

new.old.mithap <- matrix(data=c(length(subset(mithap.new$id, mithap.new$pos78=='does' | mithap.new$pos_12=='does' | mithap.new$pos_88a=='does')), (length(mithap.new$id) - length(subset(mithap.new$id, mithap.new$pos78=='does' | mithap.new$pos_12=='does' | mithap.new$pos_88a=='does'))), length(subset(mithap.old$id, mithap.old$pos78=='does' | mithap.old$pos_12=='does' | mithap.old$pos_88a=='does')), (length(mithap.old$id) - length(subset(mithap.old$id, mithap.old$pos78=='does' | mithap.old$pos_12=='does' | mithap.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.mithap) # between new and old genes from mitotic+haploid groups there is no difference in the number of genes with and without positive selection markers (p = 0.163)

# considering only haploid group genes:
hap <- subset(mydivas, mydivas$HaploidGroup=='haploid_group')
hap.new <- subset(hap, hap$age=='new')
hap.old <- subset(hap, hap$age=='old')

new.old.hap <- matrix(data=c(length(subset(hap.new$id, hap.new$pos78=='does' | hap.new$pos_12=='does' | hap.new$pos_88a=='does')), (length(hap.new$id) - length(subset(hap.new$id, hap.new$pos78=='does' | hap.new$pos_12=='does' | hap.new$pos_88a=='does'))), length(subset(hap.old$id, hap.old$pos78=='does' | hap.old$pos_12=='does' | hap.old$pos_88a=='does')), (length(hap.old$id) - length(subset(hap.old$id, hap.old$pos78=='does' | hap.old$pos_12=='does' | hap.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.hap) # no significative difference in genes with or without positive selection between new and old genes from haploid group (p = 0.6663)

# considering only mitotic group genes:
mit <- subset(mydivas, mydivas$Group=='Mitotic')
mit.new <- subset(mit, mit$age=='new')
mit.old <- subset(mit, mit$age=='old')

new.old.mit <- matrix(data=c(length(subset(mit.new$id, mit.new$pos78=='does' | mit.new$pos_12=='does' | mit.new$pos_88a=='does')), (length(mit.new$id) - length(subset(mit.new$id, mit.new$pos78=='does' | mit.new$pos_12=='does' | mit.new$pos_88a=='does'))), length(subset(mit.old$id, mit.old$pos78=='does' | mit.old$pos_12=='does' | mit.old$pos_88a=='does')), (length(mit.old$id) - length(subset(mit.old$id, mit.old$pos78=='does' | mit.old$pos_12=='does' | mit.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.mit) # for mitotic genes, there is no difference in the number of genes with positive selection markers for new and old genes [small sample size] (p = 0.4199)

# testing between haploid group and mitotic group genes
mit.hap.test <- matrix(data=c(length(subset(hap$id, hap$pos78=='does' | hap$pos_12=='does' | hap$pos_88a=='does')), (length(hap$id) - length(subset(hap$id, hap$pos78=='does' | hap$pos_12=='does' | hap$pos_88a=='does'))), length(subset(mit$id, mit$pos78=='does' | mit$pos_12=='does' | mit$pos_88a=='does')), (length(mit$id) - length(subset(mit$id, mit$pos78=='does' | mit$pos_12=='does' | mit$pos_88a=='does')))), ncol=2)

fisher.test(mit.hap.test) # there is a significative difference in the number of genes with positive slection markers when comparing haploid group and mitotic genes (p = 3.042e-06), haploid group has more positively selected genes.

old.mit.hap.test <- matrix(data=c(length(subset(hap.old$id, hap.old$pos78=='does' | hap.old$pos_12=='does' | hap.old$pos_88a=='does')), (length(hap.old$id) - length(subset(hap.old$id, hap.old$pos78=='does' | hap.old$pos_12=='does' | hap.old$pos_88a=='does'))), length(subset(mit.old$id, mit.old$pos78=='does' | mit.old$pos_12=='does' | mit.old$pos_88a=='does')), (length(mit.old$id) - length(subset(mit.old$id, mit.old$pos78=='does' | mit.old$pos_12=='does' | mit.old$pos_88a=='does')))), ncol=2)

fisher.test(old.mit.hap.test) # there is a significative difference in the number of genes with and without positive selection markers for old genes from mitotic and haploid groups (p = 7.766e-06)


######################################################################################################################################

#### considering only autosomal genes:
a.mydivas <- subset(mydivas, mydivas$XorA=="A")
new.a.mydivas <- subset(a.mydivas, a.mydivas$age=='new')
old.a.mydivas <- subset(a.mydivas, a.mydivas$age=='old')

new.old.a <- matrix(data=c(length(subset(new.a.mydivas$id, new.a.mydivas$pos78=='does' | new.a.mydivas$pos_12=='does' | new.a.mydivas$pos_88a=='does')), (length(new.a.mydivas$id) - length(subset(new.a.mydivas$id, new.a.mydivas$pos78=='does' | new.a.mydivas$pos_12=='does' | new.a.mydivas$pos_88a=='does'))), length(subset(old.a.mydivas$id, old.a.mydivas$pos78=='does' | old.a.mydivas$pos_12=='does' | old.a.mydivas$pos_88a=='does')), (length(old.a.mydivas$id) - length(subset(old.a.mydivas$id, old.a.mydivas$pos78=='does' | old.a.mydivas$pos_12=='does' | old.a.mydivas$pos_88a=='does')))), ncol=2)

fisher.test(new.old.a) # between new and old autossomal genes from all classes there is no greater number of genes with positive selection in any age group (p = 0.2635)

# no V, impossible or equal
interest.a <- subset(a.mydivas, a.mydivas$Group!="TheV" | a.mydivas$Group!="Impossible" | a.mydivas!="Equal")
interest.a.old <- subset(interest.a, interest.a$age=='old')
interest.a.new <- subset(interest.a, interest.a$age=='new')

new.old.a.interest <- matrix(data=c(length(subset(interest.a.new$id, interest.a.new$pos78=='does' | interest.a.new$pos_12=='does' | interest.a.new$pos_88a=='does')), (length(interest.a.new$id) - length(subset(interest.a.new$id, interest.a.new$pos78=='does' | interest.a.new$pos_12=='does' | interest.a.new$pos_88a=='does'))), length(subset(interest.a.old$id, interest.a.old$pos78=='does' | interest.a.old$pos_12=='does' | interest.a.old$pos_88a=='does')), (length(interest.a.old$id) - length(subset(interest.a.old$id, interest.a.old$pos78=='does' | interest.a.old$pos_12=='does' | interest.a.old$pos_88a=='does')))), ncol=2)
new.old.a.interest.p <- prop.table(new.old.a.interest, 2)*100

fisher.test(new.old.a.interest) # between new and old genes from interest groups (all but equal, impossible e theV), there is no difference in the number of genes with positive selection markers (p = 0.2635)


# only mitotic and haploid groups (interest groups)
mithap.a <- subset(a.mydivas, a.mydivas$Group=="Mitotic" | a.mydivas$HaploidGroup=="haploid_group")
mithap.a.new <- subset(mithap.a, mithap.a$age=='new')
mithap.a.old <- subset(mithap.a, mithap.a$age=='old')

new.old.a.mithap <- matrix(data=c(length(subset(mithap.a.new$id, mithap.a.new$pos78=='does' | mithap.a.new$pos_12=='does' | mithap.a.new$pos_88a=='does')), (length(mithap.a.new$id) - length(subset(mithap.a.new$id, mithap.a.new$pos78=='does' | mithap.a.new$pos_12=='does' | mithap.a.new$pos_88a=='does'))), length(subset(mithap.a.old$id, mithap.a.old$pos78=='does' | mithap.a.old$pos_12=='does' | mithap.a.old$pos_88a=='does')), (length(mithap.a.old$id) - length(subset(mithap.a.old$id, mithap.a.old$pos78=='does' | mithap.a.old$pos_12=='does' | mithap.a.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.a.mithap) # no age group is enriched with positively selected genes for autosomal genes from mitotic+haploid group (p = 0.6087)

# considering only haploid group genes:
hap.a <- subset(a.mydivas, a.mydivas$HaploidGroup=='haploid_group')
hap.a.new <- subset(hap.a, hap.a$age=='new')
hap.a.old <- subset(hap.a, hap.a$age=='old')

new.old.a.hap <- matrix(data=c(length(subset(hap.a.new$id, hap.a.new$pos78=='does' | hap.a.new$pos_12=='does' | hap.a.new$pos_88a=='does')), (length(hap.a.new$id) - length(subset(hap.a.new$id, hap.a.new$pos78=='does' | hap.a.new$pos_12=='does' | hap.a.new$pos_88a=='does'))), length(subset(hap.a.old$id, hap.a.old$pos78=='does' | hap.a.old$pos_12=='does' | hap.a.old$pos_88a=='does')), (length(hap.a.old$id) - length(subset(hap.a.old$id, hap.a.old$pos78=='does' | hap.a.old$pos_12=='does' | hap.a.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.a.hap) # no significative difference in autosomal genes with or without positive selection between groups [sample sizes are very small] (p = 1)

# considering only mitotic group genes:
mit.a <- subset(a.mydivas, a.mydivas$Group=='Mitotic')
mit.a.new <- subset(mit.a, mit.a$age=='new')
mit.a.old <- subset(mit.a, mit.a$age=='old')

new.old.a.mit <- matrix(data=c(length(subset(mit.a.new$id, mit.a.new$pos78=='does' | mit.a.new$pos_12=='does' | mit.a.new$pos_88a=='does')), (length(mit.a.new$id) - length(subset(mit.a.new$id, mit.a.new$pos78=='does' | mit.a.new$pos_12=='does' | mit.a.new$pos_88a=='does'))), length(subset(mit.a.old$id, mit.a.old$pos78=='does' | mit.a.old$pos_12=='does' | mit.a.old$pos_88a=='does')), (length(mit.a.old$id) - length(subset(mit.a.old$id, mit.a.old$pos78=='does' | mit.a.old$pos_12=='does' | mit.a.old$pos_88a=='does')))), ncol=2)

fisher.test(new.old.a.mit) # no significative difference in autossomal genes with or without positive selection between new and old genes from mitotic group [sample sizes are very small, zero for new genes with positive selection] (p = 1)

# testing between haploid group and mitotic group genes
mit.hap.a.test <- matrix(data=c(length(subset(hap.a$id, hap.a$pos78=='does' | hap.a$pos_12=='does' | hap.a$pos_88a=='does')), (length(hap.a$id) - length(subset(hap.a$id, hap.a$pos78=='does' | hap.a$pos_12=='does' | hap.a$pos_88a=='does'))), length(subset(mit.a$id, mit.a$pos78=='does' | mit.a$pos_12=='does' | mit.a$pos_88a=='does')), (length(mit.a$id) - length(subset(mit.a$id, mit.a$pos78=='does' | mit.a$pos_12=='does' | mit.a$pos_88a=='does')))), ncol=2)
mit.hap.a.test.p <- prop.table(mit.hap.a.test, 2)*100

fisher.test(mit.hap.a.test) # there are more autosomal genes with signatures of poositive selection in the haploid group then in the mitotic group (p = 2.795e-06)

old.mit.hap.a.test <- matrix(data=c(length(subset(hap.a.old$id, hap.a.old$pos78=='does' | hap.a.old$pos_12=='does' | hap.a.old$pos_88a=='does')), (length(hap.a.old$id) - length(subset(hap.a.old$id, hap.a.old$pos78=='does' | hap.a.old$pos_12=='does' | hap.a.old$pos_88a=='does'))), length(subset(mit.a.old$id, mit.a.old$pos78=='does' | mit.a.old$pos_12=='does' | mit.a.old$pos_88a=='does')), (length(mit.a.old$id) - length(subset(mit.a.old$id, mit.a.old$pos78=='does' | mit.a.old$pos_12=='does' | mit.a.old$pos_88a=='does')))), ncol=2)
old.mit.hap.a.test.p <- prop.table(old.mit.hap.a.test, 2)*100

fisher.test(old.mit.hap.a.test) # and it holds for old autosomal genes. since new genes usually have more signatures of positive selection, even when w remove them, we still have a higher amount of haploid genes with signatures of positive selection (p = 4.397e-06)

