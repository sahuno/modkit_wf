# rm -rf mPileup_HapModStrandSpecific .snakemake

SAMPLES = list(config["samples"].keys())
WKDIR = config["workdir"]

motif = "CG"
strands = ["negative", "positive"]
haplotypes = ["1", "2", "ungrouped"]
mods = ["h", "m"]

# set_species = hg38
# set_species = config["species"]

# uncomment this to run the module as a standalone
#should be followed by tabix_mpileup

rule all:
    input:
        #expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}_{{haplotypes}}_{{mods}}_CG0_{{strands}}.bedgraph", sample=SAMPLES, mods=mods, haplotypes=haplotypes, strands=strands)
        expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}.done", sample=SAMPLES)
        # expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}_{{haplotypes}}.bed", sample=SAMPLES)
        #expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}_{{haplotypes}}_{{mods}}_CG0_{{strands}}.sorted.bedgraph", sample=SAMPLES, mods=mods, haplotypes=haplotypes, strands=strands)



rule mPileup_HapModStrandSpecific:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        # unsortedBedgraphs=WKDIR + "/mPileup_HapModStrandSpecific/{sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.bedgraph",
        #sortedBedgraphs=WKDIR +"/mPileup_HapModStrandSpecific/{sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.sorted.bedgraph",
        done=WKDIR + "/mPileup_HapModStrandSpecific/{sample}/{sample}.done"
    params:
        # reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        reference_genome=config["ref"],
        modkit_threads=16,
        modkit_prob_percentiles=0.1,
        partition_tag="HP",
        region_arg = lambda wc: (
            f"--region {config['pileupRegion']}"
            if config.get("pileupRegion") else ""
        ),
        outtype_arg = lambda wc: (
            f"--bedgraph"
            if config.get("pileupBedgraph") else ""
        ),
        outdir = lambda wc: f"{WKDIR}/mPileup_HapModStrandSpecific/{wc.sample}",
        log=lambda wc: f"{WKDIR}/mPileup_HapModStrandSpecific/{wc.sample}/{wc.sample}.log"
    # log:
    #     log=WKDIR + "/log/mPileup_HapModStrandSpecific/{wc.sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.log"
    # #  threads: 16
    shell:
        """
        mkdir -p {params.outdir}

        ~/.cargo/bin/modkit pileup --threads {params.modkit_threads} \
        {params.outtype_arg} \
        --prefix {wildcards.sample} \
        --motif CG 0 \
        --ref {params.reference_genome} \
        --partition-tag {params.partition_tag} \
        {params.region_arg} \
        --log-filepath {params.log} \
        {input} \
        {params.outdir} && touch {output.done}
        """

# 

# rule sortmPileup_HapModStrandSpecific:
#     input:
#         unsortedBedgraphs=WKDIR + "/mPileup_HapModStrandSpecific/{sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.bedgraph"
#     output:
#         sortedBedgraphs=WKDIR +"/mPileup_HapModStrandSpecific/{sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.sorted.bedgraph",
#         #done=WKDIR + "/mPileup_HapModStrandSpecific/{sample}/{sample}.done"
#     params:
#         # reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
#         outdir = lambda wc: f"{WKDIR}/mPileup_HapModStrandSpecific/{wc.sample}",
#         log=lambda wc: f"{WKDIR}/mPileup_HapModStrandSpecific/{wc.sample}/{wc.sample}.sort.log"
#     # log:
#     #     log=WKDIR + "/log/mPileup_HapModStrandSpecific/{wc.sample}/{sample}_{haplotypes}_{mods}_CG0_{strands}.log"
#     # #  threads: 16
#     shell:
#         """
#         mkdir -p {params.outdir}
#         sort -k1,1 -k2,2n {input.unsortedBedgraphs} > {output.sortedBedgraphs} 2> {params.log}
#         """


    # sort -k1,1 -k2,2n {output.unsortedBedgraphs} > {output.sortedBedgraphs}

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/mPileup_HapModStrandSpecific.smk --cores 32 --directory . --configfile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal -np

