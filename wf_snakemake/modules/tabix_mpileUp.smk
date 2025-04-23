# rm -rf tabixSorted .snakemake

SAMPLES = list(config["samples"].keys())
WKDIR = config["workdir"]
tabixSortedWkdir = os.path.join(WKDIR, "mPileup_HapModStrandSpecific") # this change this to run multiple times anywhere

motif = "CG"
strands = ["negative", "positive"]
haplotypes = ["1", "2", "ungrouped"]
# mods = ["h", "m"]

# uncomment this to run the module as a standalone
rule all:
    input:
        expand(f"{tabixSortedWkdir}/{{sample}}/tabixSorted/{{sample}}_{{haplotypes}}.sorted.bed", sample=SAMPLES, haplotypes=haplotypes),
        expand(f"{tabixSortedWkdir}/{{sample}}/tabixSorted/{{sample}}_{{haplotypes}}.sorted.bed.gz", sample=SAMPLES, haplotypes=haplotypes)

rule mPileup_HapModStrandSpecific:
    input:
        unsorted_bed = lambda wc: os.path.join(tabixSortedWkdir,wc.sample, f"{wc.sample}_{wc.haplotypes}.bed")
    output:
        sortedBed=tabixSortedWkdir +"/{sample}/tabixSorted/{sample}_{haplotypes}.sorted.bed",
        tabixSorted=tabixSortedWkdir +"/{sample}/tabixSorted/{sample}_{haplotypes}.sorted.bed.gz"
    params:
        outdir = lambda wc: f"{tabixSortedWkdir}/{wc.sample}/tabixSorted"
    shell:
        """
        mkdir -p {params.outdir}
        /data1/greenbab/users/ahunos/apps/htslib-1.20/bgzip -k {input.unsorted_bed} --output {output.sortedBed}
        /data1/greenbab/users/ahunos/apps/htslib-1.20/tabix -p bed {output.sortedBed}
        """

    
    # sort -k1,1 -k2,2n {output.unsortedBedgraphs} > {output.sortedBedgraphs}
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/modules/tabix_mpileUp.smk --cores 32 --directory . --configfile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal -np



# unsortedBedgraphs = unsorted_bed = lambda wc: os.path.join(WKDIR,"mPileup_HapModStrandSpecific",wc.sample, f"{wc.sample}_{wc.haplotypes}_{wc.mods}_CG0_{wc.strands}.bedgraph")

