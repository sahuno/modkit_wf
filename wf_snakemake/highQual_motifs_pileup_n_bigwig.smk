# Load multiple config files
# configfile: ["/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_markDup.yaml",
#     "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml",
#     "/data1/greenbab/database/db_config.yaml"]

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/highQual_motifs_pileup_n_bigwig.smk --cores 32 --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/slurmMinimal  --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --quiet --forceall -np



configfile: "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_markDup.yaml"
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml"
configfile: "/data1/greenbab/database/db_config.yaml"


# Workflow parameters
outdir = config.get("modpileup_outdir", "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/unphased_DNAme")
fasta  = config["genomes_references"]["ref_hg38"]
sizes  = config["genome_sizes"]["gsize_hg38"]
region = config.get("region", "")
minCov = config.get("minCov", 10)

# Sample -> BAM mapping
modbam_dict = config["samples"]
samples     = list(modbam_dict.keys())

rule all:
    input:
        # Generate pileup.gz and index for each sample
        expand(f"{outdir}/results/{{sample}}/{{sample}}.pileup_minCov{minCov}.bed.gz", sample=samples),
        expand(f"{outdir}/results/{{sample}}/{{sample}}.pileup_minCov{minCov}.bed.gz.tbi", sample=samples),
        # BigWig tracks: a, hm, h, m
        expand(f"{outdir}/results/{{sample}}/{{sample}}.pileup_{{track}}_minCov{minCov}.bw", sample=samples, track=["a","hm","h","m"])

rule pileup_and_tracks_minCov:
    input:
        modbam=lambda wc: config["samples"][wc.sample]
    output:
        pileup   = f"{outdir}/results/{{sample}}/{{sample}}.pileup_minCov{minCov}.bed.gz",
        tabix    = f"{outdir}/results/{{sample}}/{{sample}}.pileup_minCov{minCov}.bed.gz.tbi",
        a_bw     = f"{outdir}/results/{{sample}}/{{sample}}.pileup_a_minCov{minCov}.bw",
        hm_bw    = f"{outdir}/results/{{sample}}/{{sample}}.pileup_hm_minCov{minCov}.bw",
        h_bw     = f"{outdir}/results/{{sample}}/{{sample}}.pileup_h_minCov{minCov}.bw",
        m_bw     = f"{outdir}/results/{{sample}}/{{sample}}.pileup_m_minCov{minCov}.bw"
    params:
        region_flag=lambda wc: f"--region {region}" if region else "",
        mincov=minCov,
        fasta=fasta,
        sizes=sizes,
        modkit_path="~/.cargo/bin/modkit",
        bgzip_path="/data1/greenbab/users/ahunos/apps/htslib-1.20/bgzip",
        tabix_path="/data1/greenbab/users/ahunos/apps/htslib-1.20/tabix"
    singularity:
        "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"
    log:
        a_log = f"{outdir}/results/{{sample}}/{{sample}}.bedmethyl_a_minCov{minCov}.log",
        hm_log= f"{outdir}/results/{{sample}}/{{sample}}.bedmethyl_hm_minCov{minCov}.log",
        h_log = f"{outdir}/results/{{sample}}/{{sample}}.bedmethyl_h_minCov{minCov}.log",
        m_log = f"{outdir}/results/{{sample}}/{{sample}}.bedmethyl_m_minCov{minCov}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.pileup})
        modkit pileup --threads 6 --ref {params.fasta} {input.modbam} - \
          | sort -k1,1 -k2,2n \
          | awk -v min_cov={params.mincov} '$5 >= min_cov' \
          | tee >(modkit bedmethyl tobigwig - {output.a_bw} --mod-code a -g {params.sizes} --log {log.a_log} --negative-strand-values --suppress-progress) \
          | tee >(modkit bm          tobigwig - {output.hm_bw} --mod-codes h,m -g {params.sizes} --log {log.hm_log} --negative-strand-values --suppress-progress) \
          | tee >(modkit bm          tobigwig - {output.h_bw} --mod-codes h   -g {params.sizes} --log {log.h_log}  --negative-strand-values --suppress-progress) \
          | tee >(modkit bm          tobigwig - {output.m_bw} --mod-codes m   -g {params.sizes} --log {log.m_log}  --negative-strand-values --suppress-progress) \
          | sort -k1,1 -k2,2n | bgzip -c > {output.pileup}
        tabix -p bed {output.pileup}
        """
# {params.region_flag}
#--region "chr19:70343-524369"
