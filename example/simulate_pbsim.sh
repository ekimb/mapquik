# pbsim (version 1) always adds some errors, so it can't simulate exactly perfect reads, but it can simulate with variable error rates
# note: pbsim (version 1) only has a 1% precision for error rate (can't simulate 0.5% error rate for example)

/pasteur/sonic/homes/rchikhi/tools/PBSIM-PacBio-Simulator/src/pbsim \
       ecoli.genome.100k.fa \
       --model_qc  /pasteur/sonic/homes/rchikhi/tools/PBSIM-PacBio-Simulator/data/model_qc_clr \
       --accuracy-mean 0.99\
       --accuracy-sd 0\
       --depth 10\
       --prefix nearperfect-ecoli.10X.100k

#paftools.js pbsim2fq <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]
paftools.js pbsim2fq ecoli.genome.100k.fa.fai nearperfect-ecoli.10X.100k_0001.maf > nearperfect-ecoli.10X.100k.fa

rm -rf  nearperfect-ecoli.10X.100k_0001*

