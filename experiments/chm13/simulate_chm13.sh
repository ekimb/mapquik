# pbsim (version 1) always adds some errors, so it can't simulate exactly perfect reads, but it can simulate with variable error rates
# note: pbsim (version 1) only has a 1% precision for error rate (can't simulate 0.5% error rate for example)

ref=chm13.genome.fa
reads=nearperfect-chm13.10X

pbsim \
       $ref \
       --model_qc  ~/model_qc_clr \
       --accuracy-mean 0.99\
       --accuracy-sd 0\
       --depth 10\
       --prefix $reads\
       --length-mean 12000 #hifi

#paftools.js pbsim2fq <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]
paftools.js pbsim2fq $ref.fai "$reads"_000*.maf > $reads.fa

rm -rf  "$reads"_000*
