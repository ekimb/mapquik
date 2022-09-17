K=$1
L=$2
D=$3
echo "------------- hifimap -------------"

/usr/bin/time cargo run --release -- nearperfect-chm13.10X.24kb.fa --reference chm13.genome.fa -k $K -l $L -d $D -p hifimap-$K-$L-$D --threads 11 > hifimap-$K-$L-$D.out
paftools.js mapeval hifimap-$K-$L-$D.paf
tail -4 hifimap-$K-$L-$D.out
UN=$(cat hifimap-$K-$L-$D.unmapped.out | wc -l)
echo "$UN reads unmapped."

