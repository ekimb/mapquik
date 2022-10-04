
for k in 11 15 20 25 50 100 200 500
do
#\time ~/tools/Jellyfish/bin/jellyfish  count -m $k -s 3000000000 -t 10 -o mer_counts_k$k.jf ../chm13v2.0.fa
#LD_LIBRARY_PATH=/pasteur/appa/homes/rchikhi/tools/Jellyfish/lib/ \time ~/tools/Jellyfish/query_per_sequence/query_per_sequence mer_counts_k$k.jf ../HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed |head -n 100000 > jellyfish-k$k.reads
#rm -f mer_counts_k$k.jf 
python figure1-middle-panel-jellyfish.py  jellyfish-k$k.reads >  jellyfish-k$k.reads_stats &
done


wait
