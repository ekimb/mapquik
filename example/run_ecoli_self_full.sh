echo "hifimap -------------"


/usr/bin/time cargo run --release -- nearperfect-ecoli.10X.fa --reference nearperfect-ecoli.10X.fa -k 6 -d 0.1 -l 22 -p hifimap.self --threads 8
paftools.js mapeval hifimap.self.paf

echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -X -ax map-pb -t 8 nearperfect-ecoli.10X.fa nearperfect-ecoli.10X.fa > minimap2.self.paf
paftools.js mapeval minimap2.self.paf
