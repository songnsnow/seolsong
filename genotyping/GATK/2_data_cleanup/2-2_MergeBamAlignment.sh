#MergeBamAlignment -- need dict. file of ref fasta
picard MergeBamAlignment \
    ALIGNED=/PATH/outs/Aligned.sortedByCoord.out.bam \
    UNMAPPED=/PATH/outs/Unaligned.bam \
    O=/PATH/outs/Merged.bam \
    R=/reference/Homo_sapiens_assembly38.fasta \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    VALIDATION_STRINGENCY=SILENT