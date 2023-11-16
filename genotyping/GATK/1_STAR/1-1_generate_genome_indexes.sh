#2. generate genome indexes--------------------------------------------------------
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /data/project/RCC_HWS/SS/bulk/STAR/indices \
--genomeFastaFiles /home/songnsnow/download/STAR_run/reference/Homo_sapiens_assembly38.fasta \
--sjdbGTFfile /home/songnsnow/download/STAR_run/reference/gencode.v32.primary_assembly.annotation.gtf \
--sjdbOverhang 100
#---------------------------------------------------------------------------------
