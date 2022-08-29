K=$1
L=$2
D=$3
G=$4
echo "------------- hifimap -------------"

/usr/bin/time cargo run --release -- nearperfect-chm13.10X.fa --reference chm13.genome.fa -k $K -l $L -d $D -g $G -p hifimap-$K-$L-$D-$G --threads 48 > hifimap-$K-$L-$D-$G.out
paftools.js mapeval hifimap-$K-$L-$D-$G.paf

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq nearperfect-chm13.10X.fa hifimap-$K-$L-$D-$G.unmapped.out > hifimap-$K-$L-$D-$G.unmapped.fq

echo "------------- minimap2 -------------"

/usr/bin/time minimap2 -ax asm20 -t 48 chm13.genome.fa hifimap-$K-$L-$D-$G.unmapped.fq > minimap2-$K-$L-$D-$G.unmapped.paf
paftools.js mapeval minimap2-$K-$L-$D-$G.unmapped.paf

