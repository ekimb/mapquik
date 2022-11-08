# pbsim (version 1) always adds some errors, so it can't simulate exactly perfect reads, but it can simulate with variable error rates
# note: pbsim (version 1) only has a 1% precision for error rate (can't simulate 0.5% error rate for example)

ref=chm13.genome.fa
reads=nearperfect-chm13.10X.24kb

pbsim \
       $ref \
       --model_qc  example/model_qc_clr \
       --accuracy-mean 0.99\
       --accuracy-sd 0\
       --depth 10\
       --prefix $reads\
       --length-mean 24000 #hifi

samtools faidx $ref
paftools.js pbsim2fq $ref.fai "$reads"_*.maf > $reads.fa
rm -f "$reads"_*.maf "$reads"_*.ref "$reads"_*.fastq
