#MarkDuplicates
picard MarkDuplicates \
    I=/PATH/outs/Merged.bam \
    O=/PATH/outs/Mark.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    METRICS_FILE=CSBMark.metrics