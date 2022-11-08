ref=maize.genome.fa
reads=simulated-maize.30X

pbsim \
       $ref \
       --model_qc  example/model_qc_clr \
       --accuracy-mean 0.99\
       --accuracy-sd 0\
       --depth 30\
       --prefix $reads\
       --length-mean 24000 #hifi

samtools faidx $ref
paftools.js pbsim2fq $ref.fai "$reads"_*.maf > $reads.fa
rm -f "$reads"_*.maf "$reads"_*.ref "$reads"_*.fastq