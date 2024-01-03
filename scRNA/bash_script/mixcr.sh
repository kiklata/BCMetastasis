organ="Lymph"
study="atac"
fastq1="~/scRNA_BC_metastases/Lymph/atac/SRR12695061_1.fastq.gz"
fastq2="~/scRNA_BC_metastases/Lymph/atac/SRR12695061_2.fastq.gz"

cd scRNA_BC_metastases/$organ/$study
mkdir mixcr
cd mixcr

mixcr align \
    -p rna-seq -s hs \
    --tag-pattern '^(CELL:N{16})(UMI:N{10}) \ ^(R1:*)' \
    -OvParameters.geneFeatureToAlign=VTranscript \
    -OallowPartialAlignments=true \
    -OallowNoCDR3PartAlignments=true \
    --report align.report \
    ~/scRNA_BC_metastases/Lymph/atac/SRR12695061_1.fastq.gz \
    ~/scRNA_BC_metastases/Lymph/atac/SRR12695061_2.fastq.gz \
    alignments.vdjca

mixcr correctAndSortTags --report barcode.correction.report alignments.vdjca alignments.corrected.vdjca
mixcr assemblePartial alignments.corrected.vdjca alignments.partial.assemble.vdjca
mixcr assemble  -a alignments.partial.assemble.vdjca clones.clna
mixcr assembleContigs clones.clna clones.clns

mixcr exportClones --split-by-tag CELL -tag CELL -uniqueTagCount UMI -p full clones.clns clones.clns.txt


