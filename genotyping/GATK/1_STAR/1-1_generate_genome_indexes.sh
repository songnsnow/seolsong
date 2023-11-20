#1. generate genome indexes--------------------------------------------------------
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /PATH/indices \
--genomeFastaFiles /reference/Homo_sapiens_assembly38.fasta \
--sjdbGTFfile /reference/gencode.v32.primary_assembly.annotation.gtf \
--sjdbOverhang 100
#---------------------------------------------------------------------------------
