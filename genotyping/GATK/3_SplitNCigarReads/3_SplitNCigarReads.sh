#SplitNCigarReads------------------------------------------------------------------
gatk SplitNCigarReads \
    -R /reference/Homo_sapiens_assembly38.fasta \
    -I /PATH/outs/Mark.bam \
    -O /PATH/outs/Split.bam
#----------------------------------------------------------------------------------
