# rm -rf mPileup_HapModStrandSpecific .snakemake

SAMPLES = list(config["samples"].keys())
WKDIR = config["workdir"]

motif = "CG"
strands = ["negative", "positive"]
haplotypes = ["1", "2", "ungrouped"]
mods = ["h", "m"]

# uncomment this to run the module as a standalone
rule all:
    input:
        #expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}_{{haplotypes}}_{{mods}}_CG0_{{strands}}.bedgraph", sample=SAMPLES, mods=mods, haplotypes=haplotypes, strands=strands)
        expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/{{sample}}.done", sample=SAMPLES),
        expand(f"{WKDIR}/mPileup_HapModStrandSpecific/{{sample}}/sorted/{{sample}}_{{haplotypes}}_{{mods}}_CG0_{{strands}}.sorted.bedgraph", sample=SAMPLES, mods=mods, haplotypes=haplotypes, strands=strands)



rule mPileup_HapModStrandSpecific:
    input:
        unsortedBedgraphs = unsorted_bed = lambda wc: os.path.join(WKDIR,"mPileup_HapModStrandSpecific",wc.sample, f"{wc.sample}_{wc.haplotypes}_{wc.mods}_CG0_{wc.strands}.bedgraph")
    output:
        sortedBedgraphs=WKDIR +"/mPileup_HapModStrandSpecific/{sample}/sorted/{sample}_{haplotypes}_{mods}_CG0_{strands}.sorted.bedgraph"
    params:
        outdir = lambda wc: f"{WKDIR}/mPileup_HapModStrandSpecific/{wc.sample}/sorted"
    shell:
        """
        mkdir -p {params.outdir}
        sort -k1,1 -k2,2n {input.unsortedBedgraphs} > {output.sortedBedgraphs}
        """
    
    # sort -k1,1 -k2,2n {output.unsortedBedgraphs} > {output.sortedBedgraphs}
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/sort_mPileup.smk --cores 32 --directory . --configfile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal -np
# unsortedBedgraphs = unsorted_bed = lambda wc: os.path.join(WKDIR,"mPileup_HapModStrandSpecific",wc.sample, f"{wc.sample}_{wc.haplotypes}_{wc.mods}_CG0_{wc.strands}.bedgraph")

