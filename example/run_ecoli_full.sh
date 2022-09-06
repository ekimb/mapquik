echo "hifimap -------------"


/usr/bin/time cargo run --release -- nearperfect-ecoli.10X.fa --reference ecoli.genome.fa -k 8 -d 0.01 -l 16 -p hifimap -g 100 --threads 11
paftools.js mapeval hifimap.paf

echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -ax asm20 -t 8 ecoli.genome.fa nearperfect-ecoli.10X.fa > minimap2.paf
paftools.js mapeval minimap2.paf
