marker_genes_dict = list(
  "KC-Basal" = c("KRT15", "COL17A1"),
  "KC-Suprabasal" = c("KRT10","KRT1","CALN1"),
  "KC-Spinous" = c("EPS8L1","SPINK5","IL18"),
  "KC-Proliferating" = c("DIAPH3", "ASPM"),
  "KC-Cycling" = c("KRT6A","S100A8"),
  "Mast cell" = c("HDC","CPA3"),
  "Melanocyte" = c("TRPM1", "DCT"),
  "Endo-Vas" = c("SELE","ADGRL4"), 
  "Endo-Lymph" = c('VWF','PECAM1',"PKHD1L1", "CCL21"), 
  "PVL" = c("ACTA2"),
  "Fibroblast" = c("PDGFRA","COL1A1", "COL3A1"), 
  "Schwann" = c("NRXN1", "CADM2"),
  "Eccrine gland" = c("EDAR"),
  "DC" = c("CD1C", "CCL22", "CSF2RA"),  
  "Macrophage" = c("CD163","MS4A7","SPP1","FOLR2"),
  "LC" = c("PRKCB","CD207"),
  "NK" = c("GNLY","CD96","KLRC1"),
  "T-CD8" = c("GZMA","GZMK","IFNG","SCML4"),
  "T-CD4" = c('PTPRC','CD2','CD3G',"ITK","ICOS","CD40LG")
)

dotplot_cl_order = c('T-CD4','T-CD8','NK','Macrophage','DC','LC','Mast cell',
                     'KC-Basal','KC-Suprabasal','KC-Spinous','KC-Proliferating','KC-Cycling',
                     'Eccrine gland','Melanocyte',
                     'Endo','Endo-Lymph','Endo-Vas','PVL','Fibroblast',
                     'Schwann')

kera_marker_genes_dict = list(
  "KC-Basal-COL17A1" = c("COL17A1", "KRT15",'SOX6'),  
  "KC-Basal-ITGA6" = c("ITGA6",'S100A2','TNC'),
  "KC-Suprabasal-DSC3" = c("DSC3","SLC24A3"),
  "KC-Suprabasal-KRT6A" = c('KRT6A'),
  "KC-Suprabasal-LYPD3" = c('LYPD3','KRT6B','KRT17'),
  "KC-Suprabasal-PLD1" = c('PLD1','KRT10','KRT1'),
  "KC-Spinous-AZGP1" = c('AZGP1'),
  "KC-Spinous-SPINK5" = c('SPINK5','STK10'),
  "KC-Proliferating-DIAPH3" = c('DIAPH3','ASPM','CENPP'),
  "KC-Cycling" = c('CELF2','PITPNC1')
)

mye_marker_genes_dict = list(
  "DC-ETV6" = c("ETV6", "CSF2RA"),  
  "DC-WDFY4" = c("WDFY4",'TOX'),
  "DC-WNT5B" = c("WNT5B","ADAM12"),
  "LC-DIAPH3" = c('DIAPH3','ASPM','KIF18B'),
  "LC-PRKCB" = c('CDH20','PRKCB'),
  "Macrophage-ARHGAP15" = c('DCC','ARHGAP15','CCDC88C','APOE','TREM2'),
  "Macrophage-CHIT1" = c('CHIT1','SPP1'),
  "Macrophage-F13A1" = c('F13A1','CD163','SLC9A9','FOLR2'),
  "Macrophage-RORA" = c('RORA','EGFR'),
  "Macrophage-TCF4" = c('TCF4')
)

fib_marker_genes_dict = list(
  "Fib-ANXA1" = c("ANXA1", "CCDC39",'MFAP5','DCN'),  
  "Fib-ASAP1" = c("ASAP1",'RUNX1','IL1R1'),
  "Fib-BNC2" = c("BNC2","EBF2",'SLC22A3','PDZRN4'),
  "Fib-COL23A1" = c('COL23A1','PLCB1','F13A1','COL21A1'),
  "Fib-IGFBP7" = c('IGFBP7','C7','ITGA8'),
  "Fib-ITPR2" = c('ITPR2','MAGI1','SOX5'),
  "Fib-TCF4" = c('TCF4','MKX')
)

endo_marker_genes_dict = list(
  "Endo-Lymph-DOCK5" = c("DOCK5", "PKHD1L1",'SEMA3A','COLEC12'),  
  "Endo-Vas-ARL15" = c("ARL15",'BTNL9'),
  "Endo-Vas-CD44" = c("CD44","ANK3"),
  "Endo-Vas-FLG" = c('FLG','FLG2','CALML5'),
  "Endo-Vas-MCTP1" = c('MCTP1','ZNF385D','ACKR1','SELE'),
  "Endo-Vas-PPARG" = c('PPARG','FABP4','MCC')
)

t_marker_genes_dict = list(
  "T-CD4" = c("IL7R","THEMIS",'INPP4B','LEF1'),
  "T-CD8" = c("ITGA4","GZMA","GZMK","IFNG"),
  "T-CD4-FOXP3" = c("FOXP3","CTLA4",'TOX','IL2RA'),
  "NK" = c('GNLY','KLRD1','NKG7')
)

mye_dotplot_cl_order = c('Macrophage-TCF4','Macrophage-RORA','Macrophage-F13A1','Macrophage-CHIT1','Macrophage-ARHGAP15',
                         'DC-ETV6','DC-WDFY4','DC-WNT5B',
                         'LC-DIAPH3','LC-PRKCB')

kera_dotplot_cl_order = c('KC-Basal-COL17A1','KC-Basal-ITGA6',
                          'KC-Suprabasal-DSC3','KC-Suprabasal-KRT6A','KC-Suprabasal-LYPD3','KC-Suprabasal-PLD1',
                          'KC-Spinous-AZGP1','KC-Spinous-SPINK5', 
                          'KC-Proliferating-DIAPH3','KC-Cycling')

fib_dotplot_cl_order = c('Fib-ANXA1','Fib-ASAP1','Fib-BNC2','Fib-COL23A1','Fib-IGFBP7','Fib-ITPR2','Fib-TCF4')

endo_dotplot_cl_order = c('Endo-Lymph-DOCK5','Endo-Vas-ARL15','Endo-Vas-CD44','Endo-Vas-FLG','Endo-Vas-MCTP1','Endo-Vas-PPARG')

t_dotplot_cl_order = c('T-CD4','T-CD4-FOXP3','T-CD8','NK')
