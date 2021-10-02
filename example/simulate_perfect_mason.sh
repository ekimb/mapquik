# simulate perfect reads using mason

# crashes..  length(intervals) == 1u was: 0 != 1
# and dunno how to fix

mason_simulator  \
               --fragment-size-model normal \
               --fragment-min-size 5000 \
               --fragment-max-size 5100 \
               --fragment-mean-size 5050 \
               --fragment-size-std-dev 2 \
               --seq-technology sanger \
               --sanger-prob-mismatch-begin 0 \
               --sanger-prob-mismatch-end 0 \
               --sanger-prob-insertion-begin 0 \
               --sanger-prob-insertion-end 0 \
               --sanger-prob-deletion-begin 0 \
               --sanger-prob-deletion-end 0 \
               --sanger-read-length-error 0 \
               --sanger-read-length-min 5000  \
               --sanger-read-length-mean 5001  \
               --sanger-read-length-max 5002  \
       -ir ecoli.genome.100k.fa -n 113 -oa perfect-ecoli.10X.100k.sam -o osef.fastq

rm -f osef.fastq

paftools.js mason2fq perfect-ecoli.10X.100k.sam |seqtk seq -A > perfect-ecoli.10X.100k.fa

