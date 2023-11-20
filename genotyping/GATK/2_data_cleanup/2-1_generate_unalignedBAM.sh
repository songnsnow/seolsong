#make unalignedBAM for mergebam -> revertsam
picard RevertSam \
    I=/PATH/outs/Aligned.sortedByCoord.out.bam \
    O=/PATH/outs/Unaligned.bam
