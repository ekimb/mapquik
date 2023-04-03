k=$1
for l in {10..31}
do
echo "----- Running for $k, $l -----"

/usr/bin/time cargo run --release -- ~/seqs/reads/simulated-chm13v2.0-10X.fa --reference ~/seqs/refs/chm13v2.0.oneline.fa -k $k -l $l -c 4 -s 11 -g 2000 -p mapquik-$k-$l --threads 10 > mapquik-$k-$l.out 

# get first line in format [total Q60] [wrong Q60]

paftools.js mapeval mapquik-$k-$l.paf | cut -f 3,4 > mapquik-$k-$l.stats
tail mapquik-$k-$l.out
cat mapquik-$k-$l.stats
done
