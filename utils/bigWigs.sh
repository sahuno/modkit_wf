# 1) Ensure you have an index
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

# 2) Prepend "chr" to every header line
sed 's/^>/>chr/' Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  > GRCh38_primary_withChrPrefix.fa

# 3) Build a new index
samtools faidx GRCh38_primary_withChrPrefix.fa

# 4) (Optional) extract only chr20 now if you like:
samtools faidx GRCh38_primary_withChrPrefix.fa chr20 > GRCh38_chr20_withChrPrefix.fa
samtools faidx GRCh38_chr20_withChrPrefix.fa

# 5) Rerun entropy
modkit entropy \
  --in-bam "${mod_bam}" \
  --ref GRCh38_primary_withChrPrefix.fa \
  --cpg \
  -o output_entropy_bedgraph \
  --threads 32


outdir=${PWD}
# you can run pileup and stream the output through the new command '-' means output to stdout
modkit pileup ${modbam} - \
  | tee >(modkit bedmethyl tobigwig - ${outdir}/subsample_pileup_a.bw --mod-code a -g ${sizes} --log ${outdir}/bm_a.log --negative-strand-values --suppress-progress) \
  | tee >(modkit bm tobigwig - ${outdir}/subsample_pileup_hm.bw --mod-codes h,m -g ${sizes} --log ${outdir}/bm_hm.log --negative-strand-values --suppress-progress) \
  > ${pileup}

ref=${HOME}/sta/apps/methyl_wk/resources/genomes/hg38/ensembl/GRCh38_chr20_withChrPrefix.fa
mod_bam=${HOME}/sta/DNAme/softwares/modkit/tests/resources/HG002_small.ch20._other.sorted.bam

# you can use a file also, if you have conflicting bases (e.g. 6mA and 5mC) - you'll get an error.
modkit pileup --ref ${ref} --cpg ${mod_bam} pileup.bed



# you can run pileup and stream the output through the new command '-' means output to stdout
modkit pileup ${modbam} - \
  # --mod-code a will make a track of 6mA
  | tee >(modkit bedmethyl tobigwig - ${outdir}/subsample_pileup_a.bw --mod-code a -g ${sizes} --log ${outdir}/bm_a.log --negative-strand-values --suppress-progress) \
  # you can use "bm" or "bedmethyl"
  # --mod-code (or --mod-codes h,m) will combine the 5hmC and 5mC counts into a track
  | tee >(modkit bm tobigwig - ${outdir}/subsample_pileup_hm.bw --mod-codes h,m -g ${sizes} --log ${outdir}/bm_hm.log --negative-strand-values --suppress-progress) \
  # finally pipe out to the pileup
  > ${pileup}
  
#output bigwisgs for vidualization
modkit bm tobigwig ${pileup} ${outbw} --mod-code m -g ${sizes} --negative-strand-values --log-filepath ${outdir}/bm_m.log
modkit bm tobigwig ${pileup} ${outbw} --mod-code m,h -g ${sizes} --negative-strand-values --log-filepath ${outdir}/bm_mh.log
modkit bm tobigwig ${pileup} ${outbw} --mod-code a -g ${sizes} --negative-strand-values --log-filepath ${outdir}/bm_a.log
