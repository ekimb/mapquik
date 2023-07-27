k=$1
l=$2
for d in {1..9}
do
echo "----- Running for $k, $l, 0.00$d -----"

/usr/bin/time cargo run --release -- ~/seqs/reads/simulated-chm13v2.0-10X.fa --reference ~/seqs/refs/chm13v2.0.oneline.fa -k $k -l $l -d 0.00$d -c 4 -s 11 -g 2000 -p mapquik-$k-$l-0.00$d --threads 10 > mapquik-$k-$l-0.00$d.out 

# get first line in format [total Q60] [wrong Q60]

paftools.js mapeval mapquik-$k-$l-0.00$d.paf | cut -f 3,4 > mapquik-$k-$l-0.00$d.stats
tail mapquik-$k-$l-0.00$d.out
cat mapquik-$k-$l-0.00$d.stats
done
for d in {0..9}
do
echo "----- Running for $k, $l, 0.01$d -----"

/usr/bin/time cargo run --release -- ~/seqs/reads/simulated-chm13v2.0-10X.fa --reference ~/seqs/refs/chm13v2.0.oneline.fa -k $k -l $l -d 0.01$d -c 4 -s 11 -g 2000 -p mapquik-$k-$l-0.01$d --threads 10 > mapquik-$k-$l-0.01$d.out 

# get first line in format [total Q60] [wrong Q60]

paftools.js mapeval mapquik-$k-$l-0.01$d.paf | cut -f 3,4 > mapquik-$k-$l-0.01$d.stats
tail mapquik-$k-$l-0.01$d.out
cat mapquik-$k-$l-0.01$d.stats
done