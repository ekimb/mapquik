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

echo "------------- winnowmap2 -------------"

meryl count k=15 threads=11 output chm13.db chm13.genome.fa
meryl print greater-than distinct=0.9998 chm13.db > chm13.repetitive.k15.txt
/usr/bin/time winnowmap -W chm13.repetitive.k15.txt -t 11 -x map-pb  chm13.genome.fa mapquik-$K-$L-$D.unmapped.fq > winnowmap2-$K-$L-$D.unmapped.paf
paftools.js mapeval winnowmap2-$K-$L-$D.unmapped.paf

