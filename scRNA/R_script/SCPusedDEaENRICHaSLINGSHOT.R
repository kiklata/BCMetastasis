tcell_count <- readRDS("~/Project/MultiOmics/data/snRNA/Object/summary/Immune/tcell_count.rds")
tcell_cluster <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/tcell_cluster.tsv", row.names=1)
tcell_count = AddMetaData(tcell_count,tcell_cluster)

all.count =tcell_count@assays$RNA@counts

mt.name = grep ("^MT-", rownames(all.count),value = T)
rp.name = grep("^RP[L|S]",rownames(all.count),value = T)
hb.name = grep("^HB[^(P)]",rownames(all.count),value = T)

new.count = all.count[!rownames(all.count) %in% hb.name,]

new_tcell = CreateSeuratObject(new.count, meta.data = tcell_count@meta.data)

new_tcell = NormalizeData(new_tcell,assay = 'RNA')

# DE -----------------------------------------------------------------
library(SCP)

cd4t_sub <- RunDEtest(srt = new_tcell %>% subset(.,Tcell_minor == 'CD4-Treg-FOXP3'), group_by = "SampleTimepoint",group1 = 'S2', 
                      fc.threshold = 1, only.pos = FALSE)
cd8t_sub <- RunDEtest(srt = new_tcell %>% subset(.,Tcell_minor == 'CD8-Trm-ZNF683'), group_by = "SampleTimepoint", group1 = 'S2', 
                      fc.threshold = 1, only.pos = FALSE)
voltheme = theme(
  plot.title = element_text(hjust = 0.5),
  strip.background = element_blank(),
  strip.text.x = element_blank()
  )

VolcanoPlot(srt = cd4t_sub)+labs(title = 'CD4-Treg-FOXP3')+voltheme
  
VolcanoPlot(srt = cd8t_sub)+labs(title = 'CD8-Trm-ZNF683')+voltheme



# myeloid -----------------------------------------------------------------
