dsk=~/tools/dsk/buildk2048/bin/dsk

#for k in 11 13 15 17 19 21 25 31 41 51 61 81 101 151 201 251 301 501 1001 2001
for k in 1001 2001
do
    mkdir ~/scratch/tmp-dsk
       \time $dsk -kmer-size $k -abundance-min 1 -histo 1 -file chm13v2.0.fa -out-tmp ~/scratch/tmp-dsk/
    rm -Rf ~/scratch/tmp-dsk
       mv chm13v2.0.histo chm13v2.0-$k.histo
       rm -f chm13v2.0.h5
done
