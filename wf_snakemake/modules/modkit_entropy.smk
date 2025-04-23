# rm -r .snakemake 


configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
WKDIR = config["workdir"]

# uncomment this to run the module as a standalone
rule all:
    input:
        expand(f"{WKDIR}/modkit_entropy/{{sample}}/{{sample}}_regions.bed", sample=SAMPLES),
        expand(f"{WKDIR}/modkit_entropy/{{sample}}/{{sample}}_windows.bedgraph", sample=SAMPLES)


rule modkit_entropy:
    input:
        bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        bed=WKDIR + "/modkit_entropy/{sample}/{sample}_regions.bed",
        bedgraph=WKDIR + "/modkit_entropy/{sample}/{sample}_windows.bedgraph"
        #log=WKDIR + "/modkit_entropy/{sample}/modkit_entropy.log"
    threads: config["threads"]
    log:
        WKDIR + "/modkit_entropy/{sample}/{sample}.modkit_entropy.log"
    params:
        outdir = lambda wc: f"{WKDIR}/modkit_entropy/{wc.sample}",
        ModMotifs = config["modkit_entropyMotifs"],
        ref = config["ref"],
        region_arg = lambda wc: (
            f"--regions {config['regions']}"
            if config.get("regions") else ""
        ),
        extra_args = lambda wc: config.get("modkit_entropy_extraArgs", "")
        # combined_args = lambda wc: f"{params.extra_args(wc)}"  # New line to combine extra args
    shell:
        """
        mkdir -p {params.outdir}
        modkit entropy \
            --in-bam {input.bam} \
            -o {params.outdir} \
            --prefix {wildcards.sample} \
            {params.region_arg} \
            {params.ModMotifs} \
            {params.extra_args} \
            --ref {params.ref} \
            --threads {threads} \
            --log-filepath {log}
        """ 

#--cpg \


# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/modkit_entropy.smk --cores 32 --directory . --configfile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml -np

# --config workdir=./ --use-conda --conda-prefix /path/to/conda_envs --conda-frontend mamba


#NOTE: If `--regions` then `--out-bed` should be directory. files written will be
# 1. region.bed
# 2. windows.bedgraph
# if `--regions` is not provided, then `--out-bed` should be a .bed file.
