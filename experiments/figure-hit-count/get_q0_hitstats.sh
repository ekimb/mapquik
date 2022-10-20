k=$1
echo "----- Running for $k -----"

/usr/bin/time cargo run --release -- ../chm13/simulated-chm13v2.0-10X.fa --reference ../chm13/chm13v2.0.oneline.fa -k $k -c 0 -s 0 -p hifimap-default-$k --threads 10 > hifimap-default-$k.out 

# get first line in format [total Q60] [wrong Q60]

paftools.js mapeval -Q 0 hifimap-default-$k.paf | cut -f 2,8 > hifimap-default-$k-Q0-incorrect.stats

