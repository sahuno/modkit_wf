# // Description: Snakemake module for generating modkit summary statistics
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/modkit_summary.smk \
# --cores 32 --directory . \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/slurmMinimal \
# --jobs 10 --cores all --keep-going --forceall --quiet -np

# --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3"
# --configfile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/config.yaml \
# --configfile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_markDup.yaml \
configfile: "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/config.yaml"
configfile: "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_markDup.yaml"

SAMPLES = list(config["samples"].keys())
WKDIR = WKDIR = config.get("workdir", ".")

chr_names = []
with open(config["hg38_chrSizes"]) as f:
    for line in f:
        chr_name = line.strip().split()[0]
        chr_names.append(chr_name)

chr_names=chr_names[0:2]
# print(chr_names)


# uncomment this to run the module as a standalone
rule all:
    input:
        expand(f"{WKDIR}/modkit_summary/{{sample}}/{{sample}}_{{chr}}_summary.tsv", sample=SAMPLES, chr=chr_names, allow_missing=True)

rule modkit_summary:
    input:
        bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        modSummary=WKDIR + "/modkit_summary/{sample}/{sample}_{chr}_summary.tsv"
        #log=WKDIR + "/modkit_summary/{sample}/modkit_summary.log"
    threads: config["threads"]
    log:
        WKDIR + "/modkit_summary/{sample}/{sample}.{chr}.modkit_summary.log"
    # singularity: "/data1/greenbab/users/ahunos/apps/containers/onttools_0.0.2.sif"
    params:
        outdir = lambda wc: f"{WKDIR}/modkit_summary/{wc.sample}",
        region_arg = lambda wc: f"--region {wc.chr}",
        # region_arg = lambda wc: (
        #     f"--region {config['regions']}"
        #     if config.get("regions") else ""
        # ),
        extra_args = lambda wc: config.get("modkit_summary_extraArgs", ""),
        fullSummary = lambda wc: config.get("modFullSummary", "")
    resources:
        mem_mb = 256000
        # mem_mb_per_cpu: 512000
    shell:
        """
        mkdir -p {params.outdir}
            ~/.cargo/bin/modkit summary --tsv \
            {params.extra_args} \
            {params.region_arg} \
            {params.fullSummary} \
            --threads 2 \
            --interval-size 500000 \
            {input.bam} \
            --log-filepath {log} > {output.modSummary}
        """ 