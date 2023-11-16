#MergeBamAlignment -- need dict. file of ref fasta
picard MergeBamAlignment \
    ALIGNED=/data/project/RCC_HWS/SS/bulk/STAR/outs/CSB/CSBAligned.sortedByCoord.out.bam \
    UNMAPPED=/data/project/RCC_HWS/SS/bulk/STAR/outs/CSB/CSBUnaligned.bam \
    O=/data/project/RCC_HWS/SS/bulk/STAR/outs/CSB/CSBMerged.bam \
    R=/home/songnsnow/download/STAR_run/reference/Homo_sapiens_assembly38.fasta \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    VALIDATION_STRINGENCY=SILENT