# Load configuration
configfile: "config.yaml"

# Set working directory to execution dir
import os

SAMPLES = list(config["samples"].keys())# workdir: os.getcwd()
WKDIR = config["workdir"]

rule all:
    input:
        # entropy outputs:
        expand(
            os.path.join(WKDIR, "modkit_entropy/{sample}/{sample}.entropy.bedgraph"),
            sample=SAMPLES
        )

# Include the modkit entropy module
include: "modules/modkit_entropy.smk"

#optional modules
# for m in ["modkit_entropy", "another_module", "final_report"]:
#     include: f"modules/{m}.smk"

# snakemake --snakefile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/main.smk --cores 32 --configfile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/modkit_wf/wf_snakemake/config.yaml -np

# snakemake \
#   --snakefile main.smk \
#   --cores 32 \
#   --use-conda \
#   --conda-frontend mamba \
#   --conda-prefix ~/envs

