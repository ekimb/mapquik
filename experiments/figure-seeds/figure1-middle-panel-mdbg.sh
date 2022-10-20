l=31
d=0.01

for k in 2 4 6 8 10 12 15 20
do
        #\time cargo run --manifest-path ~/rust-mdbg/Cargo.toml -- ../chm13v2.0.oneline.fa --density $d -l $l -k $k  --read-stats ../HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed  --minabund 1
        #mv ../HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed.read_stats HG002-kminmer-abundance-in-CHM13v2-k$k.read_stats
        #rm -f graph-k*
        #python figure1-middle-panel-mdbg.py  HG002-kminmer-abundance-in-CHM13v2-k$k.reads >  HG002-kminmer-abundance-in-CHM13v2-k$k.read_stats &
        shuf -n 50000 HG002-kminmer-abundance-in-CHM13v2-k$k.read_stats > HG002-kminmer-abundance-in-CHM13v2-k$k.50k.read_stats
done
wait

