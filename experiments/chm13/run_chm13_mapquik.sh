K=$1
L=$2
D=$3
C=$4
S=$5
G=$6
echo "------------- mapquik -------------"

/usr/bin/time cargo run --release -- simulated-chm13v2.0-10X.fa --reference chm13v2.0.oneline.fa -k $K -l $L -d $D -c $C -s $S -g $G -p mapquik-$K-$L-$D --threads 10 > mapquik-$K-$L-$D.out
paftools.js mapeval mapquik-$K-$L-$D.paf
tail -4 mapquik-$K-$L-$D.out
