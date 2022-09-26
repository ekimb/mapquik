K=$1
L=$2
D=$3
K2=$4
L2=$5
D2=$6
echo "hifimap -------------"
/usr/bin/time cargo run --release -- ~/hifimap/experiments/chm13/nearperfect-chm13.10X.24kb.fa --reference ~/hifimap/experiments/chm13/chm13.genome.fa -k $K -l $L -d $D -p hifimap-$K-$L-$D --threads 11 > hifimap-$K-$L-$D.out

paftools.js mapeval hifimap-$K-$L-$D.paf

tail -4 hifimap-$K-$L-$D.out

UN=$(cat hifimap-$K-$L-$D.unmapped.out | wc -l)
echo "$UN reads unmapped."

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq ~/hifimap/experiments/chm13/nearperfect-chm13.10X.24kb.fa hifimap-$K-$L-$D.unmapped.out > hifimap-$K-$L-$D-$K2-$L2-$D2.fa

/usr/bin/time cargo run --release -- hifimap-$K-$L-$D-$K2-$L2-$D2.fa --reference ~/hifimap/experiments/chm13/chm13.genome.fa -k $K2 -l $L2 -d $D2 -p hifimap-$K-$L-$D-$K2-$L2-$D2 --threads 11 > hifimap-$K-$L-$D-$K2-$L2-$D2.out
paftools.js mapeval hifimap-$K-$L-$D-$K2-$L2-$D2.paf
tail -4 hifimap-$K-$L-$D-$K2-$L2-$D2.out
UN=$(cat hifimap-$K-$L-$D-$K2-$L2-$D2.unmapped.out | wc -l)
echo "$UN reads unmapped."

