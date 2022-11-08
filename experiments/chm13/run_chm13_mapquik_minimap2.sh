K=$1
L=$2
D=$3
echo "------------- mapquik -------------"

/usr/bin/time cargo run --release -- nearperfect-chm13.10X.24kb.fa --reference chm13.genome.fa -k $K -l $L -d $D -p mapquik-$K-$L-$D --threads 11 > mapquik-$K-$L-$D.out
paftools.js mapeval mapquik-$K-$L-$D.paf
tail -4 mapquik-$K-$L-$D.out
UN=$(cat mapquik-$K-$L-$D.unmapped.out | wc -l)
echo "$UN reads unmapped."

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq nearperfect-chm13.10X.24kb.fa mapquik-$K-$L-$D.unmapped.out > mapquik-$K-$L-$D.unmapped.fq

echo "------------- minimap2 -------------"

/usr/bin/time ./minimap2-2.24 -x map-hifi -t 11 chm13.genome.fa mapquik-$K-$L-$D.unmapped.fq > minimap2-$K-$L-$D.unmapped.paf
paftools.js mapeval minimap2-$K-$L-$D.unmapped.paf

