#2. run mapping-------------------------------------------------------------------
STAR --runThreadN 20 \
--genomeDir /PATH/indices \
--readFilesIn mRNA-bulk/Name_1.fastq.gz /mRNA-bulk/Name_2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /PATH/outs/Name \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 6 \
--twopassMode Basic \
--sjdbOverhang 100
<!-- --limitBAMsortRAM
--limitOutSJcollapsed -->
#----------------------------------------------------------------------------------

