echo "minimap2 -------------"

# minimap2
/usr/bin/time minimap2 -ax asm20 -t 8 chm13.genome.fa nearperfect-chm13.10X.fa > minimap2.chm13.paf
paftools.js mapeval minimap2.paf
