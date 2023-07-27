export PATH=$PATH:/home/baris/mapquik/target/release/

VCF=~/seqs/hg002-hg38.vcf
REF=~/seqs/refs/hg38.fa
SEQ=~/seqs/reads/dc.24kb.fastq
MQPRE=mapquik-hg002-hg38
MMPRE=minimap2-hg002-hg38

# run mapquik
\time mapquik $SEQ --reference $REF --prefix $MQPRE --threads 10
# extract Q0 reads
awk '$12 == 0 {print}' $MQPRE.paf | sort | uniq | gawk -F'\t' -vOFS='\t' '{ print $6, $8, $9, $1 ":" $3 "-" $4, $12, $5 }' > $MQPRE.Q0.bed
# intersect
bedtools intersect -a $MQPRE.Q0.bed -b $VCF > $MQPRE.int

# run minimap2
\time minimap2 -t 10 -x map-hifi $REF $SEQ > $MMPRE.paf
# extract Q0 reads
awk '$12 == 0 {print}' $MMPRE.paf | sort | uniq | gawk -F'\t' -vOFS='\t' '{ print $6, $8, $9, $1 ":" $3 "-" $4, $12, $5 }' > $MMPRE.Q0.bed
# intersect
bedtools intersect -a $MMPRE.Q0.bed -b $VCF > $MMPRE.int