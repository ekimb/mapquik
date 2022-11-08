
# pasteur specific paths
export PATH=$PATH:/pasteur/appa/homes/rchikhi/tools/BLEND/bin/
export PATH=$PATH:/pasteur/appa/homes/rchikhi/tools/mapquik/target/release/


# ---------------- simulated 10X


# deprecated, we don't use hard masking anymore as regular reference works even better now thanks to Baris' heuristic of keeping unique kminmers only
#LD_LIBRARY_PATH=. \time mapquik simulated-chm13v2.0-10X.fa  --reference ../human-genome/chm13v2.0.hardmasked.oneline.fa -k 5 -l 31 -d 0.007 -f 1 -p mapquik-5-31-0.007-1-hardmasked --threads 10
LD_LIBRARY_PATH=. \time mapquik simulated-chm13v2.0-10X.fa  --reference chm13v2.0.oneline.fa -p mapquik-7-31-0.01-1 --threads 10



\time blend -t 10 -x map-hifi chm13v2.0.fa simulated-chm13v2.0-10X.fa > blend-sim10X.paf

#[M::main] CMD: /pasteur/appa/homes/rchikhi/tools/BLEND/bin/blend -t 10 -x map-hifi chm13v2.0.fa simulated-chm13v2.0-10X.fa
# [M::main] Real time: 323.651 sec; CPU: 2558.931 sec; Peak RSS: 5.893 GB
# 2499.09user 60.27system 5:24.12elapsed 789%CPU (0avgtext+0avgdata 6179508maxresident)k
# 66412152inputs+487832outputs (1major+13906184minor)pagefaults 0swaps

\time minimap2 -x map-hifi -t 10 chm13v2.0.fa simulated-chm13v2.0-10X.fa > minimap-sim10X.paf


#[M::main] CMD: minimap2 -x map-hifi -t 10 chm13v2.0.fa simulated-chm13v2.0-10X.fa
#[M::main] Real time: 1323.936 sec; CPU: 12524.641 sec; Peak RSS: 10.458 GB
#12395.82user 129.61system 22:04.75elapsed 945%CPU (0avgtext+0avgdata 10965632maxresident)k
#0inputs+494512outputs (0major+72188819minor)pagefaults 0swaps


\time ~/tools/mm2-fast/mm2-fast -x map-hifi -t 10 chm13v2.0.fa simulated-chm13v2.0-10X.fa > mm2-fast-sim10X.paf

#[M::main] Real time: 1253.776 sec; CPU: 12494.347 sec; Peak RSS: 10.492 GB
#minimizer-lookup: 0 dp: 0 rmq: 0 rmq_t1: 0 rmq_t2: 0 rmq_t3: 0 rmq_t4: 0 alignment: 0 0
#12343.30user 151.85system 21:54.58elapsed 950%CPU (0avgtext+0avgdata 11001160maxresident)k


winnowmap -t 10 -W repetitive_k15.txt -x map-pb chm13v2.0.fa simulated-chm13v2.0-10X.fa > winnowmap-sim10X.paf

#[M::main] Version: 2.03, pthreads=10, omp_threads=3
#[M::main] CMD: winnowmap -t 10 -W repetitive_k15.txt -x map-pb chm13v2.0.fa simulated-chm13v2.0-10X.fa
#[M::main] Real time: 17168.046 sec; CPU: 167135.983 sec; Peak RSS: 10.696 GB


#  ------------------- DeepConsensus 30X



LD_LIBRARY_PATH=. \time mapquik HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq  --reference chm13v2.0.oneline.fa -k 7 -l 31 -d 0.01 -p mapquik-7-31-0.01-1 --threads 10

#Total execution time: 229.499133928s
#Maximum RSS: 12.61264GB
#1050.20user 78.71system 3:49.81elapsed 491%CPU (0avgtext+0avgdata 13225312maxresident)k
#165785376inputs+357920outputs (12major+6943058minor)pagefaults 0swaps



minimap2 -x map-hifi -t 10 ../human-genome/chm13v2.0.hardmasked.oneline.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq
#[M::main] Real time: 2419.934 sec; CPU: 23239.105 sec; Peak RSS: 6.499 GB
#22922.75user 316.69system 40:20.34elapsed 960%CPU (0avgtext+0avgdata 6814592maxresident)k
#373132696inputs+3591808outputs (1major+108024891minor)pagefaults 0swaps



minimap2 -ax map-hifi -t 10 ../human-genome/chm13v2.0.hardmasked.oneline.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq
#[M::main] Real time: 14933.478 sec; CPU: 147166.823 sec; Peak RSS: 17.563 GB
#146257.68user 909.35system 4:08:54elapsed 985%CPU (0avgtext+0avgdata 18415968maxresident)k
#372118208inputs+464324352outputs (1major+378197886minor)pagefaults 0swaps


minimap2 -x map-hifi -t 10 ../human-genome/chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq
#[M::main] Real time: 4375.615 sec; CPU: 42675.539 sec; Peak RSS: 10.541 GB
#42286.92user 388.88system 1:12:55elapsed 975%CPU (0avgtext+0avgdata 11052904maxresident)k
#368208696inputs+1651744outputs (1major+207294704minor)pagefaults 0swaps


minimap2 -ax map-hifi -t 10 ../human-genome/chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq
#[M::main] Real time: 14282.141 sec; CPU: 140844.753 sec; Peak RSS: 20.087 GB
#139634.98user 1209.92system 3:58:02elapsed 986%CPU (0avgtext+0avgdata 21062716maxresident)k
#368209040inputs+377585848outputs (2major+529446327minor)pagefaults 0swaps


\time blend -t 10 -x map-hifi chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed
#[M::main] Real time: 947.034 sec; CPU: 8246.362 sec; Peak RSS: 5.562 GB
#8159.24user 87.43system 15:47.37elapsed 870%CPU (0avgtext+0avgdata 5832068maxresident)k
#57663408inputs+1494392outputs (0major+44394728minor)pagefaults 0swaps

# very long
winnowmap -W repetitive_k15.txt -x map-pb chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed > winnowmap.paf
#[M::main] CMD: winnowmap -t 10 -W repetitive_k15.txt -x map-pb chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed
#[M::main] Real time: 130862.657 sec; CPU: 1293943.525 sec; Peak RSS: 8.579 GB


\time ~/tools/mm2-fast/mm2-fast -x map-hifi -t 10 chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed > mm2-fast.paf

#[M::main] Version: 2.24-r1122
#[M::main] CMD: /pasteur/appa/homes/rchikhi/tools/mm2-fast/mm2-fast -x map-hifi -t 10 chm13v2.0.fa HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.fixed
#[M::main] Real time: 4315.603 sec; CPU: 42542.181 sec; Peak RSS: 10.544 GB
#minimizer-lookup: 0 dp: 0 rmq: 0 rmq_t1: 0 rmq_t2: 0 rmq_t3: 0 rmq_t4: 0 alignment: 0 0
#42061.77user 481.00system 1:12:56elapsed 972%CPU (0avgtext+0avgdata 11055884maxresident)k


# getting number of Q60 reads from all mappers

 awk '$12 == 60 {print}' pafs/file.paf |awk '{print $1}' |sort|uniq|wc -l


# -------------------- the simulated 10X centromeres analysis


#preparation of the nocensat reads file
cat simulated-chm13v2.0-10X.fa|grep ">" |awk -F'!' '{print $2" "$3" "$4}' |sort -n> simulated-chm13v2.0-10X.bed.tmp
sed 's/ /\t/g' simulated-chm13v2.0-10X.bed.tmp |bedtools sort -i /dev/stdin > simulated-chm13v2.0-10X.bed
rm -f simulated-chm13v2.0-10X.bed.tmp
bedtools intersect -b chm13v2.0_censat_v2.0.bed -a simulated-chm13v2.0-10X.bed -wa -v |awk '{print $1"!"$2"!"$3}' > simulated-chm13v2.0-10X.nocensat.txt
grep --no-group-separator -A1 -wFf simulated-chm13v2.0-10X.censat.txt simulated-chm13v2.0-10X.fa > simulated-chm13v2.0-10X.nocensat.fa


#gathering reads not mapped at Q60
awk '$12 == 60' mapquik.paf|awk '{print $1}'> mapquik.mappedQ60.txt
wc -l mapquik.mappedQ60.txt
-> 1448212 mapquik.mappedQ60.txt
grep -vFwf mapquik.mappedQ60.txt simulated-chm13v2.0-10X.txt > mapquik.unmappedQ60.txt

grep -Fwf mapquik.unmappedQ60.txt  simulated-chm13v2.0-10X.nocensat.fa|wc -l
-> 2803
$ wc -l mapquik.unmappedQ60.txt
-> 42198 mapquik.unmappedQ60.txt

