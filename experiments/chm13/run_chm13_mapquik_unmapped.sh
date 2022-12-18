K=$1
L=$2
D=$3
K2=$4
L2=$5
D2=$6
echo "mapquik -------------"
/usr/bin/time cargo run --release -- ~/mapquik/experiments/chm13/nearperfect-chm13.10X.24kb.fa --reference ~/mapquik/experiments/chm13/chm13.genome.fa -k $K -l $L -d $D -p mapquik-$K-$L-$D --threads 11 > mapquik-$K-$L-$D.out

paftools.js mapeval mapquik-$K-$L-$D.paf

tail -4 mapquik-$K-$L-$D.out

UN=$(cat mapquik-$K-$L-$D.unmapped.out | wc -l)
echo "$UN reads unmapped."

echo "------------- generating FASTA from unmapped reads -------------"

seqtk subseq ~/mapquik/experiments/chm13/nearperfect-chm13.10X.24kb.fa mapquik-$K-$L-$D.unmapped.out > mapquik-$K-$L-$D-$K2-$L2-$D2.fa

/usr/bin/time cargo run --release -- mapquik-$K-$L-$D-$K2-$L2-$D2.fa --reference ~/mapquik/experiments/chm13/chm13.genome.fa -k $K2 -l $L2 -d $D2 -p mapquik-$K-$L-$D-$K2-$L2-$D2 --threads 11 > mapquik-$K-$L-$D-$K2-$L2-$D2.out
paftools.js mapeval mapquik-$K-$L-$D-$K2-$L2-$D2.paf
tail -4 mapquik-$K-$L-$D-$K2-$L2-$D2.out
UN=$(cat mapquik-$K-$L-$D-$K2-$L2-$D2.unmapped.out | wc -l)
echo "$UN reads unmapped."

