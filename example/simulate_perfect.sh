# simulate perfect reads using pbsim

/pasteur/sonic/homes/rchikhi/tools/PBSIM-PacBio-Simulator/src/pbsim \
       ecoli.genome.100k.fa \
       --model_qc  /pasteur/sonic/homes/rchikhi/tools/PBSIM-PacBio-Simulator/data/model_qc_ccs \
       --accuracy-mean 1\
       --depth 10\
       --prefix perfect-ecoli.10X.100k

#paftools.js pbsim2fq <ref.fa.fai> <pbsim1.maf> [[pbsim2.maf] ...]
paftools.js pbsim2fq ecoli.genome.100k.fa.fai perfect-ecoli.10X.100k_0001.maf > perfect-ecoli.10X.100k.fa

rm -rf  perfect-ecoli.10X.100k_0001*

