# pbsim (version 1) always adds some errors, so it can't simulate exactly perfect reads, but it can simulate with variable error rates
# note: pbsim (version 1) only has a 1% precision for error rate (can't simulate 0.5% error rate for example)

ref=ecoli.genome.fa
reads=nearperfect-ecoli.10X

pbsim \
       $ref \
       --model_qc model_qc_clr \
       --accuracy-mean 0.99\
       --accuracy-sd 0\
       --depth 10\
       --prefix $reads\
       --length-mean 24000 #hifi

paftools.js pbsim2fq $ref.fai "$reads"_00*.maf > $reads.fa

rm -rf  "$reads"_00*

