# multi-anno stack bar plot -----------------------------------------------
library(reshape2)
library(ggnewscale)
library(ggplot2)
library(RColorBrewer)

df = sampleimmune

# df reGroup--------------------------------------------------------------
df$sample = factor(df$sample,levels = c("atac","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312",
                                        "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321",
                                        "oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2",
                                        "BoM1","BoM2",
                                        "cell.pid1","cell.pid2", "cell.pid3","GSE143423","GSE202501",
                                        "science.GSM4555888", "science.GSM4555889","science.GSM4555891",             
                                        "primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                                        "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                                        "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                                        "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                                        "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                                        "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                                        "primary.CID4530N" ,  "primary.CID4535"))

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
study.bone.annal = c("BoM1","BoM2")
study.brain.cell = c("cell.pid1","cell.pid2", "cell.pid3")
study.brain.ctm = c("GSE202501")
study.brain.gse143423 = c("GSE143423")
study.brain.science = c("science.GSM4555888", "science.GSM4555889","science.GSM4555891" )
study.lymph.emboj = c( "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312")
study.lymph.hmg = c('atac')
study.lymph.oncogenesis = c("oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2")
study.primary.ng = c("primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                     "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                     "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                     "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                     "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                     "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                     "primary.CID4530N" ,  "primary.CID4535")

pam.luminal = c('oncogenesis.B2','atac',study.bone.annal,study.lymph.emboj,'cell.pid1','cell.pid3','primary.CID3586','primary.CID3941',
                "primary.CID3948" , "primary.CID3963" , "primary.CID4040","primary.CID4066" ,"primary.CID4067", "primary.CID4290A",
                "primary.CID4398","primary.CID4461" ,   "primary.CID4463","primary.CID4471", "primary.CID4530N" ,  "primary.CID4535")
pam.her = c('oncogenesis.C2','oncogenesis.E2',"primary.CID3838" , "primary.CID3921","primary.CID45171")
pam.tnbc = c("oncogenesis.D2","oncogenesis.D3",'cell.pid2','GSE143423','primary.CID3946','primary.CID44041','primary.CID4465',
             "primary.CID4495" ,"primary.CID44971" ,  "primary.CID44991","primary.CID4513" , "primary.CID4515",'primary.CID4523')
pam.notknow = c('GSE202501',"science.GSM4555888", "science.GSM4555889","science.GSM4555891")

df$sample =df$sample[order(df$sample)]

df$site = if_else(df$sample %in% site.lymph,'Lymph node',
                  if_else(df$sample %in% site.bone,'Bone',
                          if_else(df$sample %in% site.brain,'Brain',
                                  if_else(df$sample %in% site.primary,'Primary site','na'))))

df$study = if_else(df$sample %in% study.bone.annal,'Annals',
                   if_else(df$sample %in% study.brain.cell,'Cell',
                           if_else(df$sample %in% study.brain.ctm,'CTM',
                                   if_else(df$sample %in% study.brain.gse143423,'GSE143423',
                                           if_else(df$sample %in% study.brain.science,'Science',
                                                   if_else(df$sample %in% study.lymph.emboj,'EMBOJ',
                                                           if_else(df$sample %in% study.lymph.hmg,'HMG',
                                                                   if_else(df$sample %in% study.lymph.oncogenesis,'Oncogenesis',
                                                                           if_else(df$sample %in% study.primary.ng,'NG','na')
                                                                   ))))))))

df$pamsub = if_else(df$sample %in% pam.luminal,'Luminal',
                    if_else(df$sample %in% pam.her,'HER2',
                            if_else(df$sample %in% pam.tnbc,'TNBC',
                                    if_else(df$sample %in% pam.notknow,'NA','na'))))


# number ---------------------------------------------------------------

da <- melt(df[,c(1:8,12:14)],id.vars =  c('sample','site','study','pamsub'))
da$site = factor(da$site,levels = c('Lymph node','Bone','Brain','Primary site'))
da$variable = gsub(pattern = "\\.",' ',da$variable)


# rename variable 
da$variable = if_else(da$variable %in% c('Monocyte Macrophage'),'Monocyte/Macrophage',da$variable)
da$variable = factor(da$variable,levels = c('T cell','B cell','DC','ILC','Monocyte/Macrophage',
                                            'Myelocyte','Early developed or Proliferative cell'))

da$sample = factor(da$sample,levels = c("atac","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312",
                                          "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321",
                                          "oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2",
                                          "BoM1","BoM2",
                                          "cell.pid1","cell.pid2", "cell.pid3","GSE143423","GSE202501",
                                          "science.GSM4555888", "science.GSM4555889","science.GSM4555891",             
                                          "primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                                          "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                                          "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                                          "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                                          "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                                          "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                                          "primary.CID4530N" ,  "primary.CID4535"))


p1 <- ggplot(data = da) +
  geom_col(aes(x = sample,y = value,fill = variable),
             position = position_stack(),#stat = 'identity',
           width = 0.9) +
  scale_fill_manual(values = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC","#053061"),
                    name = 'Cell Type') +
  #scale_y_continuous(labels = scales::percent_format()) +
  theme_gray(base_size = 18) +
  xlab('') + ylab('Number') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        #plot.margin = margin(t = 3,b = 3,r = 1,unit = 'cm')
        )


# site anno

anno1 <- data.frame(name1 = c('Lymph node','Bone','Brain','Primary site'),
                    xmin1 = c(1,14,16,24) - 0.45,
                    xmax1 = c(13,15,23,49) + 0.45)

anno1$name1 <- factor(anno1$name1,levels = anno1$name1)


maxy = 10000

p2 = p1 + new_scale('fill') +geom_rect(data = anno1,
                    aes(xmin = xmin1,xmax = xmax1,
                        ymin = -0.075*maxy,ymax = -0.07*maxy),fill = 'black',show.legend = F,)


# study anno

anno2 = data.frame(name2 = c('HMG','EMBOJ','Oncogenesis','Annals','Cell','GSE143423','CTM','Science','NG'),
                   xmin2 = c(1,2,9,14,16,19,20,21,24) - 0.45,
                   xmax2 = c(1,8,13,15,18,19,20,23,49) + 0.45)


anno2$name2 <- factor(anno2$name2,levels = anno2$name2)


p3 = p2 + new_scale('fill') +geom_rect(data = anno2,
                                       aes(xmin = xmin2,xmax = xmax2,
                                           ymin = -0.06*maxy,ymax = -0.04*maxy,fill = name2))+
  scale_fill_manual(values = brewer.pal(9,'Paired'),name = 'Study')

# pam50 

anno3 = data.frame(name3 = da$pamsub,
                   xmin3 = c(1:49) - 0.45,
                   xmax3 = c(1:49) + 0.45)

p4 = p3 + new_scale('fill') +geom_rect(data = anno3,
                                       aes(xmin = xmin3,xmax = xmax3,
                                           ymin = -0.03*maxy,ymax = -0.01*maxy,fill = name3))+
  scale_fill_manual(values = brewer.pal(4,'Dark2'),name = 'PAM50',limits = c('Luminal','HER2','TNBC','NA'))

p5 = p4 +
  geom_text(aes(x = 6,y = -0.1*maxy,label = 'Lymph Node'),size = 5) +
  geom_text(aes(x = 14.5,y = -0.1*maxy,label = 'Bone'),size = 5) +
  geom_text(aes(x = 20,y = -0.1*maxy,label = 'Brain'),size = 5) +
  geom_text(aes(x = 36,y = -0.1*maxy,label = 'Primary site'),size = 5)
count = p5
ggsave('conut.pdf',p5,width = 15,height = 10)



# proportion --------------------------------------------------------------


dp = df
dp$T.cell = dp$T.cell/dp$notMal
dp$B.cell = dp$B.cell/dp$notMal
dp$ILC = dp$ILC/dp$notMal
dp$DC = dp$DC/dp$notMal
dp$Monocyte.Macrophage = dp$Monocyte.Macrophage/dp$notMal
dp$Myelocyte = dp$Myelocyte/dp$notMal
dp$Early.developed.or.Proliferative.cell = dp$Early.developed.or.Proliferative.cell/dp$notMal
da <- melt(dp[,c(1:8,12:14)],id.vars =  c('sample','site','study','pamsub'))



# rename variable 
da$variable = if_else(da$variable %in% c('Monocyte Macrophage'),'Monocyte/Macrophage',da$variable)
da$variable = factor(da$variable,levels = c('T cell','B cell','DC','ILC','Monocyte/Macrophage',
                                            'Myelocyte','Early developed or Proliferative cell'))

da$sample = factor(da$sample,levels = c("atac","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312",
                                        "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321",
                                        "oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2",
                                        "BoM1","BoM2",
                                        "cell.pid1","cell.pid2", "cell.pid3","GSE143423","GSE202501",
                                        "science.GSM4555888", "science.GSM4555889","science.GSM4555891",             
                                        "primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                                        "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                                        "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                                        "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                                        "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                                        "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                                        "primary.CID4530N" ,  "primary.CID4535"))


p1 <- ggplot(data = da) +
  geom_col(aes(x = sample,y = value,fill = variable),
           position = position_stack(),#stat = 'identity',
           width = 0.9) +
  scale_fill_manual(values = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC","#053061"),
                    name = 'Cell Type') +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_gray(base_size = 18) +
  xlab('') + ylab('Number') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y = element_blank(),
        panel.background = element_blank(),
        #plot.margin = margin(t = 3,b = 3,r = 1,unit = 'cm')
  )


# site anno

anno1 <- data.frame(name1 = c('Lymph node','Bone','Brain','Primary site'),
                    xmin1 = c(1,14,16,24) - 0.45,
                    xmax1 = c(13,15,23,49) + 0.45)

anno1$name1 <- factor(anno1$name1,levels = anno1$name1)

maxy = 1

p2 = p1 + new_scale('fill') +geom_rect(data = anno1,
                                       aes(xmin = xmin1,xmax = xmax1,
                                           ymin = -0.075*maxy,ymax = -0.07*maxy),fill = 'black',show.legend = F,)


# study anno

anno2 = data.frame(name2 = c('HMG','EMBOJ','Oncogenesis','Annals','Cell','GSE143423','CTM','Science','NG'),
                   xmin2 = c(1,2,9,14,16,19,20,21,24) - 0.45,
                   xmax2 = c(1,8,13,15,18,19,20,23,49) + 0.45)


anno2$name2 <- factor(anno2$name2,levels = anno2$name2)


p3 = p2 + new_scale('fill') +geom_rect(data = anno2,
                                       aes(xmin = xmin2,xmax = xmax2,
                                           ymin = -0.06*maxy,ymax = -0.04*maxy,fill = name2))+
  scale_fill_manual(values = brewer.pal(9,'Paired'),name = 'Study')

# pam50 

anno3 = data.frame(name3 = da$pamsub,
                   xmin3 = c(1:49) - 0.45,
                   xmax3 = c(1:49) + 0.45)

p4 = p3 + new_scale('fill') +geom_rect(data = anno3,
                                       aes(xmin = xmin3,xmax = xmax3,
                                           ymin = -0.03*maxy,ymax = -0.01*maxy,fill = name3))+
  scale_fill_manual(values = brewer.pal(4,'Dark2'),name = 'PAM50',limits = c('Luminal','HER2','TNBC','NA'))

p5 = p4 +
  geom_text(aes(x = 6,y = -0.1*maxy,label = 'Lymph Node'),size = 5) +
  geom_text(aes(x = 14.5,y = -0.1*maxy,label = 'Bone'),size = 5) +
  geom_text(aes(x = 20,y = -0.1*maxy,label = 'Brain'),size = 5) +
  geom_text(aes(x = 36,y = -0.1*maxy,label = 'Primary site'),size = 5)
propor = p5





















