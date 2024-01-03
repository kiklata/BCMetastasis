## nichenetr

library(nichenetr)

ligand_target_matrix <- readRDS("~/Reference/nichenetr/ligand_target_matrix_nsga2r_final.rds")
lr_network <- readRDS("~/Reference/nichenetr/lr_network_human_21122021.rds")
weighted_networks <- readRDS("~/Reference/nichenetr/weighted_networks_nsga2r_final.rds")

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

cellranger_filter_count <- readRDS("~/Project/MultiOmics/data/snRNA/Object/summary/cellranger_filter_count.rds")
anno <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/cellranger_filter_anno.tsv", row.names=1)

seu = AddMetaData(cellranger_filter_count,anno)
seu = seu %>% NormalizeData(.)

select_cellmajor = c('CAF','Myeloid','T cell')
input_seu = seu %>% subset(., celltype_major_order %in% select_cellmajor)
Idents(input_seu) = input_seu$celltype_minor_order

receiver = c('CD8-Tex-ITM2C')
expressed_genes_receiver = get_expressed_genes(receiver, input_seu, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

sender_celltypes = c("Macrophage-TREM2")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, input_seu, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

celltype_sign = FindAllMarkers(input_seu)

geneset_oi = celltype_sign %>% group_by(.,cluster) %>% top_n(.,50,avg_log2FC) %>%  filter(cluster %in% c('Macrophage-TREM2','CD8-Tex-ITM2C') & p_val <= 0.05 ) %>% pull(gene) %>% unique()
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 20) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes",
                      color = "purple",legend_position = "bottom",
                      x_axis_position = "top",legend_title = "Regulatory potential")  +
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

