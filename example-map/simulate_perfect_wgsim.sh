# simulate perfect reads using wgsim and converts to mapeval format 

wgsim \
       ecoli.genome.100k.fa \
       perfect-ecoli.10X.100k.fq /dev/null \
       -e  0  \
       -N 113 \
       -1 10000 \
       -2 30 \
       -r 0 \
       -R 0\
       -X 0 \
       -h

# get strands
minimap2 ecoli.genome.100k.fa perfect-ecoli.10X.100k.fq > minimap2.paf

python wgsim_to_mapeval.py  perfect-ecoli.10X.100k.fq minimap2.paf |seqtk seq -A > perfect-ecoli.10X.100k.fa
rm -f perfect-ecoli.10X.100k.fq

