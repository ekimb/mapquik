rustmdbg=~/rust-mdbg/target/release/rust-mdbg

#cat chm13v2.0.fa | seqkit seq -u |seqtk seq -A > chm13v2.0.oneline.fa

l=31
d=0.007

for k in 2 4 6 8 10 12 15 20
do
       \time $rustmdbg --no-basespace --reference --minabund 1 --density $d -l $l -k $k  chm13v2.0.oneline.fa 
       awk '/^S/ {print $5}' graph-k"$k"-d0.007-l31.gfa | cut -c6- |sort -n  | uniq -c |awk '{print $2" "$1}'   >  kminmers-k"$k".hist
       rm -f *.sequences *.ec_data *.gfa
done
