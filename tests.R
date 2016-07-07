# All statistical tests made to see if expression among groups is different or not.
## Groups used are those made in final.R
wilcox.test(postmeiotic$PostMeiosis, control.1$PostMeiosis) # 0.1164
wilcox.test(postmeiotic.a$PostMeiosis, control.1.a$PostMeiosis) # 0.07363
wilcox.test(postmeiotic.x$PostMeiosis, control.1.x$PostMeiosis) # 0.952

wilcox.test(meiotic$Meiosis, control.2$Meiosis) #0.3925
wilcox.test(meiotic.a$Meiosis, control.2.a$Meiosis) #0.2091
wilcox.test(meiotic.x$Meiosis, control.2.x$Meiosis) #0.3873

wilcox.test(meioticpostmeiotic$PostMeiosis, control.3$PostMeiosis) #0.3019
wilcox.test(meioticpostmeiotic.a$PostMeiosis, control.3.a$PostMeiosis) #0.3658
wilcox.test(meioticpostmeiotic.x$PostMeiosis, control.3.x$PostMeiosis) #0.8889

wilcox.test(meioticpostmeiotic$Meiosis, control.4$Meiosis) #0.2533
wilcox.test(meioticpostmeiotic.a$Meiosis, control.4.a$Meiosis) #0.3504
wilcox.test(meioticpostmeiotic.x$Meiosis, control.4.x$Meiosis) #0.6667

wilcox.test(meiotic.meioticpostmeiotic$Meiosis, control.5$Meiosis) #0.0876
wilcox.test(meiotic.meioticpostmeiotic.a$Meiosis, control.5.a$Meiosis) #0.07714
wilcox.test(meiotic.meioticpostmeiotic.x$Meiosis, control.5.x$Meiosis) #0.9656

wilcox.test(meioticpostmeiotic.postmeiotic$PostMeiosis, control.6$PostMeiosis) #0.203
wilcox.test(meioticpostmeiotic.postmeiotic.a$PostMeiosis, control.6.a$PostMeiosis) #0.1372
wilcox.test(meioticpostmeiotic.postmeiotic.x$PostMeiosis, control.6.x$PostMeiosis) #0.9336

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$Meiosis, control.7$Meiosis) #0.1004
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$Meiosis, control.7.a$Meiosis) #0.09668
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$Meiosis, control.7.x$Meiosis) #0.7562

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$PostMeiosis, control.8$PostMeiosis) #0.125
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$PostMeiosis, control.8.a$PostMeiosis) #0.08351
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$PostMeiosis, control.8.x$PostMeiosis) #0.9239

wilcox.test(mitotic$Mitosis, control.mit$Mitosis) #0.3909
wilcox.test(mitotic.a$Mitosis, control.mit.a$Mitosis) #0.2627
wilcox.test(mitotic.x$Mitosis, control.mit.x$Mitosis) #1


################################ dNdS #####################################
wilcox.test(postmeiotic$dnds~postmeiotic$age) #0.8392
wilcox.test(control.1$dnds~control.1$age) # 0.1562
wilcox.test(postmeiotic.a$dnds~postmeiotic.a$age) #0.845
wilcox.test(control.1.a$dnds~control.1.a$age) #0.1481
wilcox.test(postmeiotic.x$dnds~postmeiotic.x$age) #only old genes
wilcox.test(control.1.x$dnds~control.1.x$age) #only old genes

wilcox.test(meiotic$dnds~meiotic$age) # 0.001096 ###OK###
wilcox.test(control.2$dnds~control.2$age) # only old genes
wilcox.test(meiotic.a$dnds~meiotic.a$age) # 0.06695
wilcox.test(control.2.a$dnds~control.2.a$age)  # only old genes
wilcox.test(meiotic.x$dnds~meiotic.x$age) #0.01046 ###OK###
wilcox.test(control.2.x$dnds~control.2.x$age)  # only old genes

wilcox.test(meioticpostmeiotic$dnds~meioticpostmeiotic$age) #0.02385 ###OK###
wilcox.test(control.3$dnds~control.3$age)  # only old genes
wilcox.test(control.4$dnds~control.4$age) # only old genes
wilcox.test(meioticpostmeiotic.a$dnds~meioticpostmeiotic.a$age) #0.08442
wilcox.test(control.3.a$dnds~control.3.a$age) # only old genes
wilcox.test(control.4.a$dnds~control.4.a$age) # only old genes
wilcox.test(meioticpostmeiotic.x$dnds~meioticpostmeiotic.x$age) #0.2857
wilcox.test(control.3.x$dnds~control.3.x$age) # only old genes
wilcox.test(control.4.x$dnds~control.4.x$age) # only old genes

wilcox.test(meiotic.meioticpostmeiotic$dnds~meiotic.meioticpostmeiotic$age) # 7.115e-05 ###OK###
wilcox.test(control.5$dnds~control.5$age) # only old genes
wilcox.test(meiotic.meioticpostmeiotic.a$dnds~meiotic.meioticpostmeiotic.a$age) # 0.01028 ###OK###
wilcox.test(control.5.a$dnds~control.5.a) # only old genes
wilcox.test(meiotic.meioticpostmeiotic.x$dnds~meiotic.meioticpostmeiotic.x$age) # 0.002956 ###OK###
wilcox.test(control.5.x$dnds~control.5.x$age) # only old genes

wilcox.test(meioticpostmeiotic.postmeiotic$dnds~meioticpostmeiotic.postmeiotic$age) #0.01002 ###OK###
wilcox.test(control.6$dnds~control.6$age) #0.1622
wilcox.test(meioticpostmeiotic.postmeiotic.a$dnds~meioticpostmeiotic.postmeiotic.a$age) #0.03406 ###OK###
wilcox.test(control.6.a$dnds~control.6.a$age) #0.16
wilcox.test(meioticpostmeiotic.postmeiotic.x$dnds~meioticpostmeiotic.postmeiotic.x$age) #0.1
wilcox.test(control.6.x$dnds~control.6.x$age)  # only old genes

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$dnds~meiotic.meioticpostmeiotic.postmeiotic$age) # 1.058e-05 ###OK###
wilcox.test(control.7$dnds~control.7$age) # 0.5694
wilcox.test(control.8$dnds~control.8$age) # 0.151
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$dnds~meiotic.meioticpostmeiotic.postmeiotic.a$age) # 0.002438 ###OK###
wilcox.test(control.7.a$dnds~control.7.a$age) # 0.7313
wilcox.test(control.8.a$dnds~control.8.a$age) # 0.1373
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$dnds~meiotic.meioticpostmeiotic.postmeiotic.x$age) # 0.0007177 ###OK###
wilcox.test(control.7.x$dnds~control.7.x$age) # only old genes
wilcox.test(control.8.x$dnds~control.8.x$age) # only old genes

wilcox.test(mitotic$dnds~mitotic$age) # 0.07325
wilcox.test(control.mit$dnds~control.mit$age) # only old genes
wilcox.test(mitotic.a$dnds~mitotic.a$age) #0.1014
wilcox.test(control.mit.a$dnds~control.mit.a$age) # only old genes
wilcox.test(mitotic.x$dnds~mitotic.x$age) # only old genes
wilcox.test(control.mit.x$dnds~control.mit.x$age) # only old genes

#entre

wilcox.test(postmeiotic$dnds, control.1$dnds) # 0.04361 ###OK###
wilcox.test(postmeiotic.a$dnds, control.1.a$dnds) # 0.007231 ###OK###
wilcox.test(postmeiotic.x$dnds, control.1.x$dnds) # 0.4321

wilcox.test(meiotic$dnds, control.2$dnds) # 0.9794
wilcox.test(meiotic.a$dnds, control.2.a$dnds) # 0.6056
wilcox.test(meiotic.x$dnds, control.2.x$dnds) # 0.7287

wilcox.test(meioticpostmeiotic$dnds, control.3$dnds) # 0.9562
wilcox.test(meioticpostmeiotic.a$dnds, control.3.a$dnds) #0.7186
wilcox.test(meioticpostmeiotic.x$dnds, control.3.x$dnds) #0.5

wilcox.test(meioticpostmeiotic$dnds, control.4$dnds) #0.9562
wilcox.test(meioticpostmeiotic.a$dnds, control.4.a$dnds) #0.7186
wilcox.test(meioticpostmeiotic.x$dnds, control.4.x$dnds) #0.5

wilcox.test(meiotic.meioticpostmeiotic$dnds, control.5$dnds) #0.6395
wilcox.test(meiotic.meioticpostmeiotic.a$dnds, control.5.a$dnds) #0.3152
wilcox.test(meiotic.meioticpostmeiotic.x$dnds, control.5.x$dnds) #0.583

wilcox.test(meioticpostmeiotic.postmeiotic$dnds, control.6$dnds) #0.2048
wilcox.test(meioticpostmeiotic.postmeiotic.a$dnds, control.6.a$dnds) #0.04701 ###OK###
wilcox.test(meioticpostmeiotic.postmeiotic.x$dnds, control.6.x$dnds) # 0.3729

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$dnds, control.7$dnds) #0.5282
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.7.a$dnds) #0.2228
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$dnds, control.7.x$dnds) #0.4485

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$dnds, control.8$dnds) # 0.4814
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$dnds, control.8.a$dnds) #0.2084
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$dnds, control.8.x$dnds) #0.5338

wilcox.test(mitotic$dnds, control.mit$dnds) #0.009922 ###OK###
wilcox.test(mitotic.a$dnds, control.mit.a$dnds) #0.004819 ###OK###
wilcox.test(mitotic.x$dnds, control.mit.x$dnds) #0.6834


#################################### Alpha ##################################
wilcox.test(postmeiotic$alpha~postmeiotic$age) #0.1779
wilcox.test(control.1$alpha~control.1$age) #0.1043
wilcox.test(postmeiotic.a$alpha~postmeiotic.a$age) #0.1667
wilcox.test(control.1.a$alpha~control.1.a$age) #0.1084
wilcox.test(postmeiotic.x$alpha~postmeiotic.x$age) # only old genes
wilcox.test(control.1.x$alpha~control.1.x$age) # only old genes

wilcox.test(meiotic$alpha~meiotic$age) # 0.0144 ###OK###
wilcox.test(control.2$alpha~control.2$age) # only old genes
wilcox.test(meiotic.a$alpha~meiotic.a$age) # 0.1896
wilcox.test(control.2.a$alpha~control.2.a$age) # only old genes
wilcox.test(meiotic.x$alpha~meiotic.x$age) #0.01905 ###OK###
wilcox.test(control.2.x$alpha~control.2.x$age) # only old genes

wilcox.test(meioticpostmeiotic$alpha~meioticpostmeiotic$age) #0.07019
wilcox.test(control.3$alpha~control.3$age) # only old genes
wilcox.test(control.4$alpha~control.4$age) # only old genes
wilcox.test(meioticpostmeiotic.a$alpha~meioticpostmeiotic.a$age) #0.1487
wilcox.test(control.3.a$alpha~control.3.a$age) # only old genes
wilcox.test(control.4.a$alpha~control.4.a$age) # only old genes
wilcox.test(meioticpostmeiotic.x$alpha~meioticpostmeiotic.x$age) #0.3333
wilcox.test(control.3.x$alpha~control.3.x$age) # only old genes
wilcox.test(control.4.x$alpha~control.4.x$age) # only old genes

wilcox.test(meiotic.meioticpostmeiotic$alpha~meiotic.meioticpostmeiotic$age) #0.001861 ###OK###
wilcox.test(control.5$alpha~control.5$age) # only old genes
wilcox.test(meiotic.meioticpostmeiotic.a$alpha~meiotic.meioticpostmeiotic.a$age) #0.0403 ###OK###
wilcox.test(control.5.a$alpha~control.5.a) # only old genes
wilcox.test(meiotic.meioticpostmeiotic.x$alpha~meiotic.meioticpostmeiotic.x$age) #0.009119 ###OK###
wilcox.test(control.5.x$alpha~control.5.x$age) # only old genes

wilcox.test(meioticpostmeiotic.postmeiotic$alpha~meioticpostmeiotic.postmeiotic$age) #0.003004 ###OK###
wilcox.test(control.6$alpha~control.6$age) #0.1058
wilcox.test(meioticpostmeiotic.postmeiotic.a$alpha~meioticpostmeiotic.postmeiotic.a$age) #0.009161 ###OK###
wilcox.test(control.6.a$alpha~control.6.a$age) #0.1105
wilcox.test(meioticpostmeiotic.postmeiotic.x$alpha~meioticpostmeiotic.postmeiotic.x$age) #0.1311
wilcox.test(control.6.x$alpha~control.6.x$age) # only old genes

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$alpha~meiotic.meioticpostmeiotic.postmeiotic$age) # 4.018e-05 ###OK###
wilcox.test(control.7$alpha~control.7$age) #0.7635
wilcox.test(control.8$alpha~control.8$age) #0.103
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha~meiotic.meioticpostmeiotic.postmeiotic.a$age) #0.002432 ###OK###
wilcox.test(control.7.a$alpha~control.7.a$age) #0.879
wilcox.test(control.8.a$alpha~control.8.a$age) #0.1074
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$alpha~meiotic.meioticpostmeiotic.postmeiotic.x$age) #0.003856 ###OK###
wilcox.test(control.7.x$alpha~control.7.x$age) # only old genes
wilcox.test(control.8.x$alpha~control.8.x$age) # only old genes

wilcox.test(mitotic$alpha~mitotic$age) #0.124
wilcox.test(control.mit$alpha~control.mit$age) # only old genes
wilcox.test(mitotic.a$alpha~mitotic.a$age) #0.1434
wilcox.test(control.mit.a$alpha~control.mit.a$age) # only old genes
wilcox.test(mitotic.x$alpha~mitotic.x$age) # only old genes
wilcox.test(control.mit.x$alpha~control.mit.x$age) # only old genes

#entre

wilcox.test(postmeiotic$alpha, control.1$alpha) #0.3123
wilcox.test(postmeiotic.a$alpha, control.1.a$alpha) #0.1356
wilcox.test(postmeiotic.x$alpha, control.1.x$alpha) #0.4359

wilcox.test(meiotic$alpha, control.2$alpha) #1
wilcox.test(meiotic.a$alpha, control.2.a$alpha) #0.6403
ilcox.test(meiotic.x$alpha, control.2.x$alpha) #0.3037

wilcox.test(meioticpostmeiotic$alpha, control.3$alpha) #0.7996
wilcox.test(meioticpostmeiotic.a$alpha, control.3.a$alpha) #0.6565
wilcox.test(meioticpostmeiotic.x$alpha, control.3.x$alpha) #0.8571

wilcox.test(meioticpostmeiotic$alpha, control.4$alpha) #0.7996
wilcox.test(meioticpostmeiotic.a$alpha, control.4.a$alpha) #0.6565
wilcox.test(meioticpostmeiotic.x$alpha, control.4.x$alpha) #0.8571

wilcox.test(meiotic.meioticpostmeiotic$alpha, control.5$alpha) #0.3873
wilcox.test(meiotic.meioticpostmeiotic.a$alpha, control.5.a$alpha) #0.2454
wilcox.test(meiotic.meioticpostmeiotic.x$alpha, control.5.x$alpha) #0.4395

wilcox.test(meioticpostmeiotic.postmeiotic$alpha, control.6$alpha) #0.6989
wilcox.test(meioticpostmeiotic.postmeiotic.a$alpha, control.6.a$alpha) #0.3602
wilcox.test(meioticpostmeiotic.postmeiotic.x$alpha, control.6.x$alpha) #0.4102

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$alpha, control.7$alpha) #0.9388
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.7.a$alpha) #0.5869
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$alpha, control.7.x$alpha) #0.3807

wilcox.test(meiotic.meioticpostmeiotic.postmeiotic$alpha, control.8$alpha) #0.9256
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.a$alpha, control.8.a$alpha) #0.6904
wilcox.test(meiotic.meioticpostmeiotic.postmeiotic.x$alpha, control.8.x$alpha) #0.3981

wilcox.test(mitotic$alpha, control.mit$alpha) #0.112
wilcox.test(mitotic.a$alpha, control.mit.a$alpha) #0.08239
wilcox.test(mitotic.x$alpha, control.mit.x$alpha) #0.9461
