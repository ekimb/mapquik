echo "hifimap -------------"


/usr/bin/time cargo run --release -- ../example-map/nearperfect-ecoli.10X.fa --reference ../example-map/nearperfect-ecoli.10X.fa -k 6 -d 0.1 -l 22 -p hifimap.self --threads 11 --fmin 2

echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -X -x map-hifi ../example-map/nearperfect-ecoli.10X.fa ../example-map/nearperfect-ecoli.10X.fa > minimap2.self.paf
