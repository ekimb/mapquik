K=$1
L=$2
D=$3
F=$4
echo "------------- hifimap -------------"

/usr/bin/time cargo run --release -- nearperfect-chm13.10X.fa --reference chm13.genome.fa -k $K -l $L -d $D -f $F -p hifimap-$K-$L-$D-$F --threads 11 > hifimap-$K-$L-$D-$F.out
paftools.js mapeval hifimap-$K-$L-$D-$F.paf

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq nearperfect-chm13.10X.fa hifimap-$K-$L-$D-$F.unmapped.out > hifimap-$K-$L-$D-$F.unmapped.fq

echo "------------- minimap2 -------------"

/usr/bin/time minimap2-2.24 -x map-hifi -t 11 chm13.genome.fa hifimap-$K-$L-$D-$F.unmapped.fq > minimap2-$K-$L-$D-$F.unmapped.paf
paftools.js mapeval minimap2-$K-$L-$D-$F.unmapped.paf

