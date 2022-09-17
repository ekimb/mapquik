K=$1
L=$2
D=$3
echo "------------- hifimap -------------"

/usr/bin/time cargo run --release -- nearperfect-chm13.10X.24kb.fa --reference chm13.genome.fa -k $K -l $L -d $D -p hifimap-$K-$L-$D --threads 11 > hifimap-$K-$L-$D.out
paftools.js mapeval hifimap-$K-$L-$D.paf
tail -4 hifimap-$K-$L-$D.out
UN=$(cat hifimap-$K-$L-$D.unmapped.out | wc -l)
echo "$UN reads unmapped."

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq nearperfect-chm13.10X.24kb.fa hifimap-$K-$L-$D.unmapped.out > hifimap-$K-$L-$D.unmapped.fq

echo "------------- minimap2 -------------"

/usr/bin/time ./minimap2-2.24 -x map-hifi -t 11 chm13.genome.fa hifimap-$K-$L-$D.unmapped.fq > minimap2-$K-$L-$D.unmapped.paf
paftools.js mapeval minimap2-$K-$L-$D.unmapped.paf

