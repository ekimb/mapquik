echo "minimap2 -------------"

# minimap2
/usr/bin/time ./minimap2-2.24 -x map-hifi -t 11 chm13.genome.fa nearperfect-chm13.10X.24kb.fa > minimap2.chm13.paf
paftools.js mapeval minimap2.chm13.paf
