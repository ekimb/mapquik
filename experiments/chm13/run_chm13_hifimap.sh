K=$1
L=$2
D=$3
F=$4
echo "hifimap -------------"


/usr/bin/time cargo run --release -- nearperfect-chm13.10X.10k.fa --reference chm13.genome.fa -k $K -l $L -d $D -f $F -p hifimap-$K-$L-$D-$F --threads 11 > hifimap-$K-$L-$D-$F.out
paftools.js mapeval hifimap-$K-$L-$D-$F.paf
