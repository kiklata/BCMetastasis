mat = propor[c(-1,-8,-10,-11,-13,-16,-17,-18,-19,-20,-21,-23,-24,-25,-29,-30,-32,-33,-44),]

site.lymph = c("atac","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312",
               "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321",
               "oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2")
site.bone = c( "BoM1","BoM2")
site.brain = c("cell.pid1","cell.pid2", "cell.pid3","GSE143423","GSE202501",
               "science.GSM4555888", "science.GSM4555889","science.GSM4555891" )
site.primary = c("primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                 "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                 "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                 "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                 "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                 "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                 "primary.CID4530N" ,  "primary.CID4535")
mat$sample = if_else(rownames(mat) %in% site.lymph,'lymph',
                     if_else(rownames(mat) %in% site.bone,'bone',
                             if_else(rownames(mat) %in% site.brain,'brain',
                                     if_else(rownames(mat) %in% site.primary,'primary','na'))))

  
library(corrplot)

mycol = rev(COL2('RdBu',200))

# ALL
mat.calc = mat[,1:10]

mat.calc=apply(mat.calc,2,function(x) as.numeric(as.character(x)))

mat.cor = cor(mat.calc,method = 'pearson')

corrplot(mat.cor,method = 'square',type = 'lower',order = 'original',
         tl.pos = 'l',tl.cex = .6,cl.pos = 'b',tl.col = 'black',#addCoef.col = 'grey',
         col = mycol,col.lim = c(-1,1),cl.length = 5,diag = T)


# primary
mat.primary = filter(mat,sample == 'primary')

mat.calc = mat.primary[,1:10]

mat.calc=apply(mat.calc,2,function(x) as.numeric(as.character(x)))

mat.cor = cor(mat.calc,method = 'pearson')

corrplot(mat.cor,method = 'square',type = 'lower',order = 'original',
         tl.pos = 'l',tl.cex = .6,cl.pos = 'b',tl.col = 'black',#addCoef.col = 'grey',
         col = mycol,col.lim = c(-1,1),cl.length = 5,diag = T)


# brain
mat.brain = filter(mat,sample == 'brain')
mat.calc = mat.brain[,1:10]

mat.calc=apply(mat.calc,2,function(x) as.numeric(as.character(x)))

mat.cor = cor(mat.calc,method = 'pearson')

corrplot(mat.cor,method = 'square',type = 'lower',order = 'original',
         tl.pos = 'l',tl.cex = .6,cl.pos = 'b',tl.col = 'black',#addCoef.col = 'grey',
         col = mycol,col.lim = c(-1,1),cl.length = 5,diag = T)


# lymph
mat.lymph = filter(mat,sample == 'lymph')
mat.calc = mat.lymph[,1:10]

mat.calc=apply(mat.calc,2,function(x) as.numeric(as.character(x)))

mat.cor = cor(mat.calc,method = 'pearson')

corrplot(mat.cor,method = 'square',type = 'lower',order = 'original',
         tl.pos = 'l',tl.cex = .6,cl.pos = 'b',tl.col = 'black',#addCoef.col = 'grey',
         col = mycol,col.lim = c(-1,1),cl.length = 5,diag = T)

# bone
# insufficient data to calc
#mat.bone = filter(mat,sample == 'bone')
#mat.calc = mat.bone[,1:10]

#mat.calc=apply(mat.calc,2,function(x) as.numeric(as.character(x)))

#mat.cor = cor(mat.calc,method = 'pearson')

#corrplot(mat.cor,method = 'square',type = 'lower',order = 'original',
#         tl.pos = 'l',tl.cex = .6,cl.pos = 'b',tl.col = 'black',#addCoef.col = 'grey',
#         col = mycol,col.lim = c(-1,1),cl.length = 5,diag = T)

mat.all = mat
mat.all[,1:10]=apply(mat.all[,1:10],2,function(x) as.numeric(as.character(x)))

mat.all$site = mat.all$sample
mat.all$sample = rownames(mat.all)
rownames(mat.all) = NULL

mycomparsion = list(c('lymph','primary'),c('brain','primary'),c('bone','primary'))
# P1
mat.all.p1 = filter(mat.all,P1<0.15)
p1 = ggplot(mat.all.p1, aes(x = site, y = P1,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p2
mat.all.p2 = filter(mat.all,P2<0.15)
p2 = ggplot(mat.all.p2, aes(x = site, y = P2,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p3
mat.all.p3 = filter(mat.all,P3<0.1)
p3 = ggplot(mat.all.p3, aes(x = site, y = P3,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p4
hist(mat.all$P4)
mat.all.p4 = filter(mat.all,P4<0.8)
p4 = ggplot(mat.all.p4, aes(x = site, y = P4,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p5
hist(mat.all$P5)
mat.all.p5 = filter(mat.all,P5<0.5)
p5 = ggplot(mat.all.p5, aes(x = site, y = P5,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p6
hist(mat.all$P6)
mat.all.p6 = filter(mat.all,P6<0.5)
p6 = ggplot(mat.all.p6, aes(x = site, y = P6,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p7
hist(mat.all$P7)
mat.all.p7 = filter(mat.all,P7<0.15)
p7 = ggplot(mat.all.p7, aes(x = site, y = P7,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p8
hist(mat.all$P8)
mat.all.p8 = filter(mat.all,P8<0.5)
p8 = ggplot(mat.all.p8, aes(x = site, y = P8,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p9
hist(mat.all$P9)
mat.all.p9 = filter(mat.all,P9<0.5)
p9 = ggplot(mat.all.p9, aes(x = site, y = P9,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  
#p10
hist(mat.all$P10)
mat.all.p10 = filter(mat.all,P10<0.5)
p10 = ggplot(mat.all.p10, aes(x = site, y = P10,color = site))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + stat_compare_means(method = 't.test',comparisons = mycomparsion)  

ggsave(filename = '~/p1.png',p1,width = 6,height = 4)
ggsave(filename = '~/p2.png',p2,width = 6,height = 4)
ggsave(filename = '~/p3.png',p3,width = 6,height = 4)
ggsave(filename = '~/p4.png',p4,width = 6,height = 4)
ggsave(filename = '~/p5.png',p5,width = 6,height = 4)
ggsave(filename = '~/p6.png',p6,width = 6,height = 4)
ggsave(filename = '~/p7.png',p7,width = 6,height = 4)
ggsave(filename = '~/p8.png',p8,width = 6,height = 4)
ggsave(filename = '~/p9.png',p9,width = 6,height = 4)
ggsave(filename = '~/p10.png',p10,width = 6,height = 4)
