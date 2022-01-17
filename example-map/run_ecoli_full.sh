echo "hifimap -------------"


/usr/bin/time cargo run --release -- nearperfect-ecoli.10X.fa --reference ecoli.genome.fa -k 6 -d 0.1 -l 22 -p hifimap --threads 11 --fmin 2
paftools.js mapeval hifimap.paf

echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -x map-hifi ecoli.genome.fa nearperfect-ecoli.10X.fa > minimap2.paf
paftools.js mapeval minimap2.paf
