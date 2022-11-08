l=31
d=0.01

for k in 2 4 6 8 10 12 15 20
do
        \time cargo run --manifest-path ~/tools/mapquik/Cargo.toml -- --reference ../chm13v2.0.oneline.fa --density $d -l $l -k $k ../HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed 
done

for f in `ls *.read_stats`; do awk '{print $2}' $f > $f.abbrv; done
for f in `ls *.read_stats`; do  sort $f.abbrv|uniq -c >$f.abbrv.hist; done

