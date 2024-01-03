# bashsed '/^all_references/d' <  infercnv.cell_groupings > trimmed_infercnv.cell_groupings 

##this is to match the chr-regions to the chr-arms
hg38 = read.delim("~/Reference/cytoBand.txt.gz", header = F)
cnv = read.table("HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",header = T) ##from inferCNV resulting files
#cytoband = hg19
cytoband = hg38
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL
names(cytoband) = c("chr", "start", "end", "arm")
cytoband$chr_arm = paste0(cytoband$chr, cytoband$arm)
cytoband$chromosome = paste0("chr", cytoband$chr)
cytoband$chr = cytoband$chromosome
cytoband = cytoband[,1:5]
library(dplyr)
cnv$chr = as.character(cnv$chr)
cytoband$chr = as.character(cytoband$chr)
cnv$mid = (as.numeric(cnv$start) + as.numeric(cnv$end))/2
new = inner_join(cnv,cytoband, by=c("chr"),relationship = "many-to-many")%>%
mutate(arm = ifelse(mid >= as.numeric(start.y) & mid <= as.numeric(end.y), chr_arm, NA))%>%
  group_by(chr)%>%
  dplyr::filter(!is.na(arm)|n()==1)
final_cnv = new[,c(1:6, 10)]
names(final_cnv) = c("cell_group_name", "cnv_name", "state", "chr", "start", "end", "arm")
final_cnv$event = ifelse(final_cnv$state >=4, "gain", "loss")
final_cnv$large_event = paste(final_cnv$arm, final_cnv$event, sep = " ")
final_cnv = final_cnv[!duplicated(final_cnv[,c('cell_group_name', 'large_event')]),]
large_events = unique(final_cnv$large_event)
final_cnv$group = gsub("sample.", "", final_cnv$cell_group_name)
##change the data frame to wide to summarize the event
library(reshape2)
wide = reshape2::dcast(final_cnv, large_event ~ group, value.var = "event", fun.aggregate = length)
##change the table to 0/1
rownames(wide) = wide$large_event
wide = wide[,-1]
wide = as.matrix(ifelse(wide >= 1, 1, 0))
wide = as.data.frame(wide)
wide = wide[colnames(wide) %>% grep('observat',.,value = T)]
wide$total = rowSums(wide)
new = wide[order(-wide$total),]
new$"1.1.1" = ifelse((new$`1.1.1.1` + new$`1.1.1.2`) >=2, 1, 0)
new$"1.1.2" = ifelse((new$`1.1.2.1` + new$`1.1.2.2`) >=2, 1, 0)
new$"1.1" = ifelse((new$`1.1.1` + new$`1.1.2`) >=2, 1, 0)
new$"1.2.1" = ifelse((new$`1.2.1.1` + new$`1.2.1.2`) >= 2, 1, 0)
new$"1.2.2" = ifelse((new$`1.2.2.1` + new$`1.2.2.2`) >= 2, 1, 0)
new$"1.2" = ifelse((new$`1.2.1` + new$`1.2.2`) >=2, 1, 0)
new$"1" = ifelse((new$`1.1` + new$`1.2`) >=2, 1, 0)
new = new[,c(14,11,9,1,2,10,3,4,13,12,5,6,7,8)]
write.csv(new, "events.csv")
##read in the cell grouping data from UPhyloplot2 to calculate the cell number in each branch
cells = read.table("12_HMM_preds.cell_groupings", header = T)
sum = as.data.frame(table(cells$cell_group_name))
write.csv(sum, "cell_group_frequency.csv")
##calculate the percentage of each CNV event
cell_group = read.csv("12_HMM_preds.cell_groupings.csv", header = F)
group = as.data.frame(te[,1])
rownames(group) = te[,2]
per = new[,-14]
percentage = group[,1][col(per)]
percentage$sample = rowSums(percentage)
write.csv(percentage, "sample_events_percentage.csv")