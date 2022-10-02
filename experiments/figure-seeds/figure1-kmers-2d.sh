dsk=~/tools/dsk/buildk2048/bin/dsk

#for k in 11 13 15 17 19 21 25 31 41 51 61 81 101 151 201 251 301 501 1001 2001
for k in 1001 2001
do
    mkdir ~/scratch/tmp-dsk
       \time $dsk -kmer-size $k -abundance-min 1 -out-tmp ~/scratch/tmp-dsk/ -histo2D 1 -file chm13v2.0.fa,../HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed -max-memory 1000000 -nb-cores 10
    rm -Rf ~/scratch/tmp-dsk
       mv chm13v2.0.histo2D chm13v2.0-$k.histo2D
       rm -f chm13v2.0.h5
done
