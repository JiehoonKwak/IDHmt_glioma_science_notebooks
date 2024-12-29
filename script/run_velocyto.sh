#!/bin/bash
transcriptome_path="/home/jiehoonk/mnt/annotations/refdata-gex-GRCm39-2024-A"
gft_path="/home/jiehoonk/mnt/annotations/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

for sample in ../data/tnl/*; do
    if [ -d "$sample" ]; then
        sample=$(basename $sample)
        echo "Running cellranger for $sample"
        cellranger count --id=${sample} \
            --transcriptome=$transcriptome_path \
            --fastqs=../data/tnl/${sample}/fastq \
            --output-dir=../data/tnl/${sample}/${sample}_cellranger_count \
            --create-bam=true \
            --nosecondary
            # --cell-annotation-model=auto \

        echo "Running velocyto for $sample"
        velocyto run10x ../data/tnl/${sample}/${sample}_cellranger_count $gft_path
    fi
done