
inhouse_path = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse'

bom_site = 'BoM'
bom_sample = c('BoM1','BoM3')
inhouse_bom_path = file.path(inhouse_path,bom_site,bom_sample,'cNMF')

brm_site = 'BrM'
brm_sample = c('BrM1','BrM2','BrM3')
inhouse_brm_path = file.path(inhouse_path,brm_site,brm_sample,'cNMF')

liverm_site = 'LiverM'
liverm_sample = c('LiverM1','LiverM2','LiverM3','LiverM4','LiverM5','LiverM6')
inhouse_liverm_path = file.path(inhouse_path,liverm_site,liverm_sample,'cNMF')

lungm_site = 'LungM'
lungm_sample = c('LungM1')
inhouse_lungm_path = file.path(inhouse_path,lungm_site,lungm_sample,'cNMF')

public_path = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public'

p_bom_sample = c('BoM1','BoM2')
p_bom_path = file.path(public_path,bom_site,p_bom_sample,'cNMF')

p_brm_sample = c('cell/BrM1','cell/BrM2','cell/BrM3','GSE143423','GSE202501','science/BrM1','science/BrM2','science/BrM3')
p_brm_path = file.path(public_path,brm_site,p_brm_sample,'cNMF')

p_lymphm_sample = c('EMBOJ/LymphM1','EMBOJ/LymphM2','EMBOJ/LymphM3','EMBOJ/LymphM4','EMBOJ/LymphM5','EMBOJ/LymphM6','EMBOJ/LymphM7',
                    'NC/LymphM1','NC/LymphM5','NC/LymphM6','NC/LymphM7','NC/LymphM8',
                    'oncogenesis/LymphM1','oncogenesis/LymphM2','oncogenesis/LymphM3','oncogenesis/LymphM4','oncogenesis/LymphM5')
p_lymphm_path = file.path(public_path,'LymphM',p_lymphm_sample,'cNMF')

p_primary_sample = c('CID3921', 'CID3941', 'CID3948', 'CID3963', 'CID4066', 'CID4067', 'CID4290A',
                     'CID4461', 'CID4463', 'CID4465', 'CID4471', 'CID4495', 'CID44971', 'CID44991',
                     'CID4513', 'CID4515', 'CID45171', 'CID4523', 'CID4530N', 'CID4535')
p_primary_path = file.path(public_path,'Primary',p_primary_sample,'cNMF')

c(inhouse_bom_path,inhouse_brm_path,inhouse_liverm_path,inhouse_lungm_path,p_bom_path,p_brm_path,p_lymphm_path,p_primary_path)
