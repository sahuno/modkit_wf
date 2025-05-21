import os, glob

configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/samples_merged_bedmethyl.yaml"
SAMPLES = config["samples"]

NORMAL_BED="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/mergedBedmethyl/DMSO/merged_DMSO.bed.gz"

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/modkit_dmr_pair.smk --cores alll --forcerun -np
REF = "/data1/greenbab/database/mm10/mm10.fa"
minCov = 5
THREADS = 8

# 0) path to your ModKit Singularity image
MODKIT_IMG = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"

# 1) load treatments from sample.yaml


rule all:
    input:
        expand("results/modkit_dmr_pair_unphased/{sample}/{sample}_unphased_raw_dmr.bed", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_raw_segmentation.bed", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr.bed.gz", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr.bed.gz.tbi", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_segmentation.bed.gz", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_segmentation.bed.gz.tbi", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed.gz", sample=SAMPLES.keys()),
        expand("results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed.gz.tbi", sample=SAMPLES.keys())


rule modkit_dmr_unphased:
    """
    Call DMRs comparing unphased normal vs tumor 5mC BED outputs,
    then emit raw and post‐processed beds, bgzip/tabix index, and filtered DMRs.
    """
    input:
        tumor_bed  = lambda wildcards: config["samples"][wildcards.sample],
        normal_bed = NORMAL_BED
    output:
        # keep raw outputs
        raw_dmr      = "results/modkit_dmr_pair_unphased/{sample}/{sample}_unphased_raw_dmr.bed",
        raw_seg      = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_raw_segmentation.bed",
        # compressed/indexed DMR
        dmr_gz       = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr.bed.gz",
        dmr_tbi      = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr.bed.gz.tbi",
        # compressed/indexed segmentation
        seg_gz       = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_segmentation.bed.gz",
        seg_tbi      = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_segmentation.bed.gz.tbi",
        # filtered DMRs (uncompressed + compressed/indexed)
        dmr_diff     = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed",
        dmr_diff_gz  = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed.gz",
        dmr_diff_tbi = "results/modkit_dmr_unphased/{sample}/{sample}_unphased_TumorNormal_dmr_diff.bed.gz.tbi"
    params:
        ref           = REF,
        base          = "C",
        mincov        = minCov,
        min_dmr_sites = 3
    threads: THREADS
    singularity: MODKIT_IMG
    log:
        "logs/modkit_dmr_pair_unphased/{sample}.log"
    shell:
        r"""
        mkdir -p results/modkit_dmr_unphased/{wildcards.sample}.log

        # 1) Raw DMR + segmentation
        modkit dmr pair \
          -a {input.normal_bed} \
          -b {input.tumor_bed} \
          -o {output.raw_dmr} \
          --segment {output.raw_seg} \
          --ref {params.ref} \
          --base {params.base} \
          --min-valid-coverage {params.mincov} \
          --threads {threads} \
          --log-filepath {log} \
          --careful \
          --header

        # 2) Compress & index raw DMR
        sort -k1,1 -k2,2n {output.raw_dmr} \
          | bgzip -c > {output.dmr_gz}
        tabix -p bed {output.dmr_gz}

        # 3) Compress & index raw segmentation
        sort -k1,1 -k2,2n {output.raw_seg} \
          | bgzip -c > {output.seg_gz}
        tabix -p bed {output.seg_gz}

        # 4) Filter for “different” DMRs
        awk -F '\t' 'NR>1 && $4=="different" && $6>={params.min_dmr_sites} && ($14>=0.5||$14<=-0.5) && ($15*$16>0)' \
          {output.raw_seg} > {output.dmr_diff}

        # 5) Compress & index filtered DMRs
        sort -k1,1 -k2,2n {output.dmr_diff} \
          | bgzip -c > {output.dmr_diff_gz}
        tabix -p bed {output.dmr_diff_gz}
        """


