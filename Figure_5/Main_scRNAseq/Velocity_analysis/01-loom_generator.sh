#!/usr/bin/env bash


gtf=/project/hnisz_seq_data/Velocyto/gencode.vM10.annotation.gtf

repeat=/project/hnisz_bioinf/Genomes/mm10/mm10_rmsk.gtf

RNA756=/project/hnisz_seq_data/Velocyto/soupX_AM-RNA_756_RNA_757_filt/features.tsv.gz
RNA738=/project/hnisz_seq_data/Velocyto/soupX_AM-RNA_738_filt/features.tsv.gz
RNA739=/project/hnisz_seq_data/Velocyto/soupX_AM-RNA-739_filt/features.tsv.gz
RNA740=/project/hnisz_seq_data/Velocyto/soupX_AM-RNA-740_filt/features.tsv.gz
RNA755=/project/hnisz_seq_data/Velocyto/soupX_AM-RNA-755_filt/features.tsv.gz

bam756=/project/hnisz_seq_data/scRNAseq/Project_MatHep/RNA-756_RNA-757/outs/possorted_genome_bam.bam
bam738=/project/hnisz_seq_data/scRNAseq/Project_MatHep/RNA-738/outs/possorted_genome_bam.bam
bam739=/project/hnisz_seq_data/scRNAseq/Project_MatHep/RNA-739/outs/possorted_genome_bam.bam
bam740=/project/hnisz_seq_data/scRNAseq/Project_MatHep/RNA-740/outs/possorted_genome_bam.bam
bam755=/project/hnisz_seq_data/scRNAseq/Project_MatHep/RNA-755/outs/possorted_genome_bam.bam

velocyto run -@ 20 -b $RNA756 -o ./loomRNA756 -m $gtf $bam756 $repeat
velocyto run -@ 20 -b $RNA738 -o ./loomRNA738 -m $gtf $bam738 $repeat
velocyto run -@ 20 -b $RNA739 -o ./loomRNA739 -m $gtf $bam739 $repeat
velocyto run -@ 20 -b $RNA740 -o ./loomRNA740 -m $gtf $bam740 $repeat
velocyto run -@ 20 -b $RNA755 -o ./loomRNA755 -m $gtf $bam755 $repeat
