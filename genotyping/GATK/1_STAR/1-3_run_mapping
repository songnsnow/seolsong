#3. run mapping-------------------------------------------------------------------
STAR --runThreadN 20 \
--genomeDir /data/project/RCC_HWS/SS/bulk/STAR/indices \
--readFilesIn /data/project/RCC_HWS/mRNA-bulk/CSB_1.fastq.gz /data/project/RCC_HWS/mRNA-bulk/CSB_2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /data/project/RCC_HWS/SS/bulk/STAR/outs/CSB/CSB \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 6 \
--twopassMode Basic \
--sjdbOverhang 100
<!-- --limitBAMsortRAM
--limitOutSJcollapsed -->
#patients: CSB, CSK, KJT, KSH, KYE, LHJ, MYS, YJH
#----------------------------------------------------------------------------------

